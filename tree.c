#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<limits.h>

#if MACOSX
#include<fcntl.h>
#include<unistd.h>
#endif

#include<sys/random.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "this.h"

#define bitspword (uint32)(8*sizeof(uint64))
#define bitspint  (uint32)(8*sizeof(uint32))

#define INVALID(cond,args,...) \
{ \
    if((cond)) \
    { \
        fprintf(stderr,args,##__VA_ARGS__); \
        exit(1); \
    } \
}

#define INVALID_ARG(cond,name,arg) \
    if((cond)) \
    { \
        fprintf(stderr,"Invalid value for " #name ": %s\n",arg); \
        exit(1); \
    }

#define FULLORPART(_i) \
    ( \
      (_i) == nwords - 1?residual:bitspword \
    ) 

#if defined(DIAG) || defined(STEPWISE)
#define ASSERT_PRINT(stmt,bits) \
{ \
    assert((stmt) || !bitprint(bits,1)); \
}

#define ASSERT_GOOD_GENOTYPE(bits) \
{ \
    uint32 _i,_j; \
    uint64 _word; \
    for(_i = 0 ; _i < nwords ; _i++) \
    { \
        _word = (bits)[_i]; \
        for(_j = 0 ; _j < FULLORPART(_i); _j += allele_size, _word >>= allele_size) \
        ASSERT_PRINT( (_word & allele_mask) < nalleles,bits); \
    } \
    ASSERT_PRINT(_word <= 1,bits); \
} 
#endif
#define ENVELOPE_DEGREE (6)

extern void parse_rates(char *s);
extern void parse_contrib(char *s);

typedef struct _bitstr
{
    uint64 *bits;
    /* weight is not initialized until inserted into tree */
    uint32 weight;
} bitstr ;

/* Globals set by cmdline arguments */
uint32 nloci;
double discount;
double shift_rate;
double shift_size;
double *mutation_rate;
uint32 *mutation_contrib;
gsl_ran_discrete_t **mutant_tables;
double sex_mutation_rate;
uint32 sex_change;

uint32 nindiv;
uint32 nwords;
uint32 residual;
uint32 nalleles;
uint32 allele_size;
uint64 allele_mask;

uint32 maximum_weight;
double env;

/* Global for location of sexuals */
int found_sex = -1;

/* High water marks */
#ifdef DIAG
uint32 hwm_cache;
uint64 mutation_events;
uint32 hwm_tree_size;
uint32 hwm_tree_depth;
#endif

/* Forward declarations */
void mutate(bitstr bs);
int bitprint(uint64 *bits,int nl);

/* Bitstring operations */
uint32 weight(uint64 *bits)
{
    uint32 i,j,res = 0;
    uint64 s;
    for(i = 0; i < nwords; i++)
        for(j = 0, s = bits[i]; j < FULLORPART(i); j += allele_size, s >>= allele_size)
            res += mutation_contrib[s & allele_mask];
    return res;
}

int cmp(bitstr b1, bitstr b2)
{
    int i;
    /* Bitstrings are ordered word-wise little endian */
    for(i = nwords - 1 ; i >= 0 ; i--)
    {
        if(b1.bits[i] == b2.bits[i])
            continue;
        return b1.bits[i] < b2.bits[i]?1:-1;
    }
    return 0;
}

void xor(uint64 *p, uint64 *q, uint64 *res)
{
    uint32 i;
    for(i=0;i<nwords;i++)
        res[i] = p[i]^q[i];
}

void xor_me(uint64 *p, uint64 *q)
{
    uint32 i;
    for(i=0;i<nwords;i++)
        p[i] ^= q[i];
}

/* Single data structure */
struct node
{
    bitstr bs;
    uint32 n;
    /* Nodes are either in the cache or in the tree,
     * never both, so unionize */
    union
    {
        struct node *left;
        /* for cache */
        struct node *next;
    };
    struct node *right;
};

int check_sex_bit(bitstr bs)
{
    return !!(bs.bits[nwords-1] & (1UL << residual));
}

void set_sex_bit(bitstr bs)
{
    bs.bits[nwords-1] |= (1UL << residual);
}

struct nodecache
{
    struct node* head;
#if defined(STEPWISE) || defined(DIAG)
    size_t size;
#endif
} ncache = {NULL
#if defined(STEPWISE) || defined(DIAG)
    ,0} ;
#else
    };
#endif

struct _tree
{
    struct node *root;
#if defined(STEPWISE) || defined(DIAG)
    size_t size;
#endif
} tree = {NULL
#if defined(STEPWISE) || defined(DIAG)
    ,0};
#else
};
#endif

void putnode(struct node *n)
{
    n->next = ncache.head;
    ncache.head = n;
#if defined(STEPWISE) || defined(DIAG)
    ncache.size++;
#ifdef DIAG
    hwm_cache = (hwm_cache < ncache.size)?ncache.size:hwm_cache;
#endif
#endif
}

/* Node and tree methods */

/* bitstring weights are computed always and only here */
struct node *getnode(bitstr bs)
{
    struct node *n;
    if(ncache.head)
    {
        n = ncache.head;
        ncache.head = n->next;
#if defined(STEPWISE) || defined(DIAG)
        assert(ncache.size);
        ncache.size--;
#endif
    }
    else
    {
        n = malloc(sizeof(*n));
        n->bs.bits = malloc(sizeof(uint64)*nwords);
    }
    memcpy(n->bs.bits,bs.bits,nwords*sizeof(uint64));
    n->bs.weight = weight(n->bs.bits);
#if defined(DIAG) || defined(STEPWISE)
    ASSERT_PRINT(n->bs.weight <= maximum_weight, n->bs.bits);
    ASSERT_GOOD_GENOTYPE(n->bs.bits);
#endif
    n->n = 1;
    n->left = n->right = NULL;
    return n;
}

void insert(bitstr bs) 
{
   struct node **pcursor = &tree.root;
   struct node *cursor;
#ifdef DIAG
   uint32 depth = 0;
#endif

   while((cursor = *pcursor))
   {
       int d = cmp(cursor->bs,bs);
       if(d < 0)
           pcursor = &cursor->left;
       else if(d > 0)
           pcursor = &cursor->right;
       else
       {
           cursor->n++;
           return;
       }
#ifdef DIAG
       depth++;
#endif
   }

   struct node *new = getnode(bs);

   *pcursor = new;
#if defined(DIAG) || defined(STEPWISE)
   tree.size++;
#ifdef DIAG
   hwm_tree_depth = 
       hwm_tree_depth >= depth?hwm_tree_depth:depth;
#endif
#endif
}

#ifdef STEPWISE
void remove_node(struct node **prev, struct node *child)
{
    struct node *n;

    if((n = child->left))
    {
        *prev = n;
        prev = &child->right;
        while(*prev)
            prev = &(*prev)->left;
        *prev = n->right;
        n->right = child->right;
    } else
        *prev = child->right;

    putnode(child);
    tree.size--;
}

void delete(bitstr bs)
{
    struct node **pcursor = &tree.root;
    struct node *cursor;

    while((cursor = *pcursor))
    {
       int d = cmp(cursor->bs,bs);
       if(d < 0)
           pcursor = &cursor->left;
       else if(d > 0)
           pcursor = &cursor->right;
       else if(cursor->n > 1)
       {
           cursor->n--;
           return;
       }
       else
       {
           remove_node(pcursor,cursor);
           return;
       }
    }
}
#endif

#if defined(STEPWISE) || defined(DIAG)
/* Printing methods */
int bitprint(uint64 *bits,int nl)
{
    uint32 i,j,b;
    uint64 word;

    for(i = 0; i < nwords ; i++)
    {
        word = bits[i];
        for(j = 0 ; j < FULLORPART(i) ; j += allele_size , word >>= allele_size)
        {
            j>0?printf("|"):0;
            for(b=0 ; b < allele_size ; b++)
                printf("%d",!!(word & ( 1UL << (allele_size - b - 1))));
        }
    }
    printf("|");
    printf("%ld",word);
    if(nl) printf("\n");
    return 1;
}

int printbf(struct node *node)
{
    int n = node->n;
    uint32 weight = node->bs.weight;
    bitprint(node->bs.bits,0);
    return printf(": %05u %05u\n",n,weight);
}

void partialdump(struct node *n)
{
    if(n)
    {
        partialdump(n->left);
        printbf(n);
        partialdump(n->right);
    }
}

void dumptree(void)
{
    printf("Tree: \n");
    partialdump(tree.root);
}
#endif

struct _array
{
    bitstr *bs;
    double *w;
    uint32 len;
    uint32 space;
} array;

/* Array methods */

/* TODO Change array to parent_array? */
void initialize_array(void)
{
    array.bs = malloc(sizeof(*array.bs)*0x1000);
    /* Bit strings in one call to malloc */
    uint64 *allbits = malloc(sizeof(uint64)*nwords*0x1000);

    int i;
    for(i=0;i<0x1000;i++,allbits += nwords)
        array.bs[i].bits = allbits;
    array.w    = malloc(sizeof(double)*0x1000);
    array.len  = 0;
    array.space = 0x1000;
}

void increase_array_space(void)
{
    array.bs = realloc(array.bs,sizeof(*array.bs)*array.space*2);
    array.w = realloc(array.w,sizeof(double)*array.space*2);

    uint64 *allbits = malloc(sizeof(uint64)*nwords*array.space);
    uint32 i;
    for(i=array.space;i<2*array.space;i++,allbits += nwords)
        array.bs[i].bits = allbits;

    array.space *= 2;
}

void append_array(struct node* n)
{
    if(array.len > (7*array.space)/8)
        increase_array_space();

    memcpy(array.bs[array.len].bits,n->bs.bits,sizeof(uint64)*nwords);
    /* XXX */
    /* Necessary? */
    array.bs[array.len].weight = n->bs.weight;
    array.w[array.len] = pow(discount,fabs((int)n->bs.weight - env))*n->n;
    if(found_sex < 0) 
        /* Two-fold advantage */
        array.w[array.len] *= 2;
    array.len++;
    assert(array.len <= nindiv);
}

#if defined(DIAG) || defined(STEPWISE)
void dump_array(void)
{
    printf("Array:\n");
    uint32 i;
    for(i = 0 ; i < array.len && (!i || printf("\n")) ; i++)
    {
        bitprint(array.bs[i].bits,0);
        printf(": %6.6f",array.w[i]);
    }
    printf("\n");
}
#endif

/* Globals */
uint32 *no_sex_weights;
uint32 *sex_weights;

void plinearize_and_tally_weights(struct node **pcursor)
{
    struct node* cursor;
    while((cursor = *pcursor))
    {
        plinearize_and_tally_weights(&cursor->left);
        *pcursor = cursor->right;
        int sex_bit = check_sex_bit(cursor->bs);
        if(found_sex < 0 && sex_bit)
            found_sex = array.len;
        append_array(cursor);
        if(found_sex >= 0)
            sex_weights[cursor->bs.weight] += cursor->n;
        else
            no_sex_weights[cursor->bs.weight] += cursor->n ;
        putnode(cursor);
#ifdef DIAG
        assert(( found_sex >= 0 && sex_bit ) ||
                ( found_sex < 0 && !sex_bit ) || !printf("%d %d\n",found_sex,sex_bit));
#endif
    }
}

void linearize_and_tally_weights(void)
{
    array.len = 0;
    found_sex = -1;
    memset(no_sex_weights,0,(maximum_weight + 1)*sizeof(int));
    memset(sex_weights,0,(maximum_weight + 1)*sizeof(int));
#if defined(STEPWISE) || defined(DIAG)
    tree.size = 0;
#endif
    plinearize_and_tally_weights(&tree.root);
    assert(!tree.root);
}

gsl_rng *rng;

void initialize_rng(void)
{
    uint64 seed;
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
#if MACOSX
    int fd;

    INVALID((fd = open("/dev/urandom",O_RDONLY)) < 0,
        "Failed to open \"/dev/urandom\"");

    INVALID(read(fd,&seed,sizeof(seed)) < (ssize_t)sizeof(seed),
#else
    INVALID(getrandom(&seed,sizeof(seed),0) < (ssize_t)sizeof(seed),
#endif
        "Failed to properly initialize random number generator\n")
#if MACOSX
    INVALID(close(fd) < 0,"Failed to close \"/dev/urandom\"");
#endif
    gsl_rng_set(rng,seed);
}

void make_children(uint64 *scratch1,uint32 *choices, uint32 choices_ints)
{
    uint32 i,j;
    bitstr res;
    gsl_ran_discrete_t *tl,*tl_sex = NULL;

    linearize_and_tally_weights();
    tl = gsl_ran_discrete_preproc(array.len,array.w);
    if(found_sex > -1)
        tl_sex = gsl_ran_discrete_preproc(array.len - found_sex,
                array.w + found_sex);

    for(i=0;i<nindiv;i++)
    {
        int k = gsl_ran_discrete(rng,tl);
        if(found_sex < 0 || k < found_sex)
        {
            memcpy(scratch1,array.bs[k].bits,sizeof(uint64)*nwords);
            res.bits = scratch1;
        }
        else
        {
            memset(scratch1,0,sizeof(uint64)*nwords);
            for(j = 0 ; j < choices_ints; j++)
                choices[j] = gsl_rng_get(rng);

            uint64 *dad = array.bs[k].bits;
            uint64 *mom = array.bs[gsl_ran_discrete(rng,tl_sex)+found_sex].bits;
            uint64 *choice;
            uint64 allele;

            for(j =  0; j < nloci ; j++)
            {
                if(choices[j/bitspint] & (1 << ( j % bitspint)))
                    choice = mom;
                else 
                    choice = dad;
                allele = choice[(j*allele_size)/bitspword] >> (j*allele_size) % bitspword;
                allele &= allele_mask;
                scratch1[(j*allele_size)/bitspword]
                    ^= allele << ((j*allele_size) % bitspword);
            }

            res.bits = scratch1;
            /* We may have cleared sex bit, so reset it */
            set_sex_bit(res);
        }
#if defined(DIAG) || defined(STEPWISE)
        ASSERT_GOOD_GENOTYPE(res.bits);
#endif
        mutate(res);
        insert(res);
    }
}

static inline double env_envelope(double x)
{
    if(x > maximum_weight || x < 0)
        return 0;
    x /= maximum_weight/2.;
    x -= 1;
    return 1 - pow(x,ENVELOPE_DEGREE);
}

void pick_new_env(void)
{
    double cand = env + gsl_ran_gaussian_ziggurat(rng,shift_size);
    double here = env_envelope(env);
    double there = env_envelope(cand);
    if(here == 0)
    {
        if(there == 0)
            env = gsl_rng_uniform(rng) < 0.5?env:cand;
        else
            env = cand;
        return;
    }
        
    double chance = there/here;
    if(!isfinite(chance))
    {
        fprintf(stderr,"Has undefined Metropolitan-Hastings probability: exiting\n");
        exit(1);
    }

    if(gsl_rng_uniform(rng) < chance)
        env = cand;
}

/* Mutation data structures */

void initialize_mutation_parameters(char *argv[])
{
    mutation_rate = malloc(sizeof(double)*nalleles*nalleles);
    parse_rates(argv[8]);

    mutation_contrib = malloc(sizeof(uint32)*nalleles);
    parse_contrib(argv[9]);

    mutant_tables = malloc(sizeof(*mutant_tables)*nalleles);
    uint32 i;
    for(i = 0 ; i < nalleles ; i++)
        mutant_tables[i] = gsl_ran_discrete_preproc(nalleles,
                mutation_rate + i*nalleles);

    uint32 max = 0; 
    for(i = 0 ; i < nalleles ; i++)
        max = mutation_contrib[i] > max?mutation_contrib[i]:max;
    maximum_weight = max*nloci;
}

void mutate(bitstr bs)
{
   uint32 i,j;
   uint64 new,scratch,scratch2;

   int sex_bit = check_sex_bit(bs);
   for(i = 0 ; i < nwords ; i++)
   {
       scratch = bs.bits[i];
       new = 0;
       for(j = 0; j < FULLORPART(i); j += allele_size,scratch >>= allele_size)
       {
           scratch2 = gsl_ran_discrete(rng,mutant_tables[scratch & allele_mask]);
#if DIAG
           if((scratch & allele_mask) != scratch2)
               mutation_events++;
#endif
           new ^= scratch2 << j;
       }
       bs.bits[i] = new;
   }

   if(sex_bit)
       set_sex_bit(bs);

   if(sex_change)
   { 
       bs.bits[nwords - 1]
           ^= (uint64)!!(gsl_rng_uniform(rng) < sex_mutation_rate) << residual;
   }
#if defined(DIAG) || defined(STEPWISE)
   ASSERT_GOOD_GENOTYPE(bs.bits);
#endif
}

int main(int argc, char *argv[])
{
    if(argc == 2 && !strcmp(argv[1],"--usage"))
    {
        fprintf(stderr,
                "./exec \\\n"
                "   nloci \\\n"
                "   nalleles \\\n"
                "   shift_rate \\\n"
                "   shift_size \\\n"
                "   discount \\\n"
                "   nosex \\\n"
                "   sex \\\n"
                "   mutation_rate \\\n"
                "   mutation_contrib \\\n"
                "   sex_mutation_rate \\\n"
                "   sex_change \\\n"
                "   [ ngen ]\n");
        exit(0);
    }

    if(argc != 
#if defined(STEPWISE)
            12
#else
            13
#endif
      )
        INVALID(1,"Wrong number of arguments\n");

    bitstr bs;

    nloci = (uint32)strtoul(argv[1],NULL,0);
    INVALID_ARG(nloci == 0,nloci,argv[1]);

    nalleles = (uint32)strtoul(argv[2],NULL,0);
    INVALID_ARG(nalleles < 2,nalleles,argv[1]);

    shift_rate = strtod(argv[3],NULL);
    INVALID_ARG(shift_rate < 0 || shift_rate > 1,shift_rate,argv[2]);

    shift_size = strtod(argv[4],NULL);
    INVALID_ARG(shift_size < 0,shift_size,argv[3]);

    discount = strtod(argv[5],NULL);
    INVALID_ARG(discount < 0 || discount > 1,discount,argv[4]);

    uint32 nosex = (uint32)strtoul(argv[6],NULL,0);
    uint32 sex = (uint32)strtoul(argv[7],NULL,0);

    sex_mutation_rate = strtod(argv[10],NULL);
    sex_change = !!strtol(argv[11],NULL,0);

    initialize_mutation_parameters(argv);

#if !defined(STEPWISE)
    uint32 ngen = (uint32)strtoul(argv[12],NULL,0);
#endif

    /* allele fit in windows of 2,4,8,16,32 or 64 bits */
    allele_size = 0;
    uint32 s = nalleles - 1;
    while(s || bitspword % allele_size)
    {
        allele_size++;
        s >>= 1;
    }

    INVALID(allele_size > bitspword,"Too many alleles\n")

    allele_mask = (1UL << allele_size) - 1;

    /* Extra padding for sex bit */
    nwords = (nloci*allele_size/bitspword) + 1;
    residual = (nloci*allele_size) % bitspword;

    nindiv = nosex + sex;
    INVALID(nindiv == 0,
            "Invalid value for population size: nosex: %s, sex: %s\n",argv[6],argv[7]);

    initialize_rng();
    initialize_array();

    /* Passed arrays for make_children */
    uint64 *scratch = malloc(sizeof(uint64)*nwords);
    uint32 *choices = malloc(bitspint*((nloci + bitspint - 1)/(8*bitspint)));
    uint32 choices_ints = (nloci + 8*bitspint - 1)/(8*bitspint);

    no_sex_weights = malloc((maximum_weight + 1)*sizeof(uint32));
    sex_weights = malloc((maximum_weight + 1)*sizeof(uint32));

    env = maximum_weight/2;

    bs.bits = malloc(sizeof(uint64)*nwords);

#ifndef STEPWISE
    /* Initialize a population randomly */
    uint32 i,j,k;
    size_t rand;
    double *uniform = malloc(sizeof(double)*nalleles);

    for(i = 0 ; i < nalleles ; i++)
        uniform[i] = 1.;
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(nalleles,uniform);

    for(i = 0; i < nindiv ; i++)
    {
        memset(bs.bits,0,sizeof(uint64)*nwords);
        for(j = 0 ; j < nwords ; j++)
        {
            for(k = 0; k < FULLORPART(j); k += allele_size)
            {
                rand = gsl_ran_discrete(rng,table);
                bs.bits[j] ^= rand << k;
            }
        }
        if(i < sex)
            set_sex_bit(bs);
#ifdef DIAG
        ASSERT_GOOD_GENOTYPE(bs.bits);
#endif
        insert(bs);
    }
    gsl_ran_discrete_free(table);
    free(uniform);
    free(bs.bits);
#endif

#if !defined(DIAG) && !defined(STEPWISE)
    for(i = 0 ; i < ngen && (!i || printf("\n")); i++)
    {
        make_children(scratch,choices,choices_ints);

        /* Print weights */
        printf("env: %8.2f\n",env);
        printf("  no_sex:");
        for(j = 0 ; j <= maximum_weight; j++)
            printf(" %u:%u",j,no_sex_weights[j]);
        printf("\n     sex:");
        for(j = 0 ; j <= maximum_weight; j++)
            printf(" %u:%u",j,sex_weights[j]);

        pick_new_env();
    }
    printf("\n");
    return 0;

#elif defined(STEPWISE)
    assert(nloci*allele_size <= 63);
    uint64 x;
    char line[1024];
    for(;;)
    {
        if(fgets(line,1024,stdin) == NULL)
            break;
        if(line[0] == '\n') { dumptree(); continue;}
        if(line[1] == '\n')
        {
            if(line[0] == 'd')
            {
                dumptree();
                linearize_and_tally_weights();
                uint32 i;
                for(i = 0; i<array.len;i++)
                {
                    printf("0x%016lx ",array.bs[i].bits[0]);
                } 
                printf("\n");
                for(i = 0; i<array.len;i++)
                {
                    printf("%f ",array.w[i]);
                }
                printf("\n");
            }
            if(line[0] == 'c')
            {
                make_children(scratch,choices, choices_ints);
                dump_array();
                printf("Found sex?: %d\n",found_sex);
                dumptree();
            }

            continue;
        }

        char *endptr;
        x = strtoul(line,&endptr,2);
        bs.bits = &x;

        if(endptr)
        {
            int i;
            long j;
            if(endptr[0] == '-')
            {
                j = strtol(endptr,NULL,0);
                for(i=0;i< -j ; i++)
                    delete(bs);
            }
            else if(endptr[0] == '+')
            {
                endptr++;
                j = strtol(endptr,NULL,0);
                for(i=0;i<j;i++)
                    insert(bs);
            }
            else
                INVALID(1,"Invalid entry\n");
        }
                
        dumptree();
        printf("    Tree size: %ld\n",tree.size);
        printf("    ncache.size = %ld\n",ncache.size);

    }

#elif defined(DIAG)
    printf("Start\n");
    for(i=0;i<ngen;i++)
    {
        dumptree();
        make_children(scratch,choices,choices_ints);
        dump_array();
        hwm_tree_size = (hwm_tree_size >= tree.size)?hwm_tree_size:tree.size;
    }

    printf("End\n");
    dumptree();
    printf("High water mark tree_size: %u\n",hwm_tree_size);
    printf("High water mark tree depth: %u\n",hwm_tree_depth);
    printf("High water mark cache_size: %u\n",hwm_cache);
    printf("Mutation events: %lu\n",mutation_events);
    return 0;
#endif
}
