#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#if MACOSX
#include<fcntl.h>
#include<unistd.h>
#endif

#include<sys/random.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "this.h"

#define bitspword (8*sizeof(uint64))
#define INVALID(cond,name,arg) \
    if((cond)) \
    { \
        fprintf(stderr,"Invalid value for " #name ": %s\n",arg); \
        exit(1); \
    }
#define ENVELOPE_DEGREE (6)

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
double *mutation_contrib;
uint32 sex_change;

uint32 nindiv;
uint32 nwords;
uint32 residual;
double env;

/* Global for location of sexuals */
int found_sex = -1;

/* High water marks */
#ifdef DIAG
uint32 hwm_cache;
uint32 hwm_mutation_sites_space;
uint32 hwm_power_table_len;
uint64 mutation_events;
uint32 hwm_tree_depth;
#endif

/* Forward declarations */
void mutate(bitstr bs);

/* Bitstring operations */
uint32 weight(uint64 *bits)
{
    uint32 i,res = 0;
    uint64 s = 0;
    for(i = 0; i < nwords - 1; i++)
        for(s=bits[i];s;s >>= 1)
            res += s&1;
    /* Don't count the sex bit */
    for(i=0,s=bits[nwords-1];s && i < residual ;s >>= 1,i++)
        res += s&1;
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

struct tree
{
    struct node *root;
#if defined(STEPWISE) || defined(DIAG)
    size_t size;
#endif
} thetree = {NULL
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
    n->n = 1;
    n->left = n->right = NULL;
    return n;
}

void insert(bitstr bs) 
{
   struct node **pcursor = &thetree.root;
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
   thetree.size++;
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
    thetree.size--;
}

void delete(bitstr bs)
{
    struct node **pcursor = &thetree.root;
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
void printbf(struct node *node)
{
    uint32 i;
    uint64 *bits = node->bs.bits;
    int n = node->n;
    int weight = node->bs.weight;
    for(i = 0 ; i < nwords ; i++)
        printf("%016lx",bits[i]);
    printf(": %05d %05d\n",n,weight);
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
    partialdump(thetree.root);
}
#endif

struct arrays
{
    bitstr *bs;
    double *w;
    uint32 len;
    uint32 space;
} thearray;

/* Array methods */

/* TODO Change array to parent_array? */
void initialize_array(void)
{
    thearray.bs = malloc(sizeof(*thearray.bs)*0x1000);
    /* Bit strings in one call to malloc */
    uint64 *allbits = malloc(sizeof(uint64)*nwords*0x1000);

    int i;
    for(i=0;i<0x1000;i++,allbits += nwords)
        thearray.bs[i].bits = allbits;
    thearray.w    = malloc(sizeof(double)*0x1000);
    thearray.len  = 0;
    thearray.space = 0x1000;
}

void increase_array_space(void)
{
    thearray.bs = realloc(thearray.bs,sizeof(*thearray.bs)*thearray.space*2);
    thearray.w = realloc(thearray.w,sizeof(double)*thearray.space*2);

    uint64 *allbits = malloc(sizeof(uint64)*nwords*thearray.space);
    uint32 i;
    for(i=thearray.space;i<2*thearray.space;i++,allbits += nwords)
        thearray.bs[i].bits = allbits;

    thearray.space *= 2;
}

void append_array(struct node* n)
{

    if(thearray.len > (7*thearray.space)/8)
        increase_array_space();

    memcpy(thearray.bs[thearray.len].bits,n->bs.bits,sizeof(uint64)*nwords);
    /* XXX */
    /* Necessary? */
    thearray.bs[thearray.len].weight = n->bs.weight;
    thearray.w[thearray.len] = pow(discount,fabs((int)n->bs.weight - env))*n->n;
    if(found_sex < 0) 
        /* Two-fold advantage */
        thearray.w[thearray.len] *= 2;
    thearray.len++;
    assert(thearray.len <= nindiv);
}

#if defined(DIAG) || defined(STEPWISE)
void dump_array(void)
{
    printf("Array:\n");
    uint32 i,j;
    for(i = 0 ; i < thearray.len && (!i || printf("\n")) ; i++)
    {
        for(j = 0 ; j < nwords ; j++)
            printf("%016lx",thearray.bs[i].bits[j]);
        printf(": %6.6f",thearray.w[i]);
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
        if(found_sex < 0 && check_sex_bit(cursor->bs))
            found_sex = thearray.len;
        append_array(cursor);
        if(found_sex >= 0)
            sex_weights[cursor->bs.weight] += cursor->n;
        else
            no_sex_weights[cursor->bs.weight] += cursor->n ;
        putnode(cursor);
#ifdef DIAG
        int sex_bit = check_sex_bit(cursor->bs);
        assert(( found_sex >= 0 && sex_bit ) ||
                ( found_sex < 0 && !sex_bit ) );
#endif
    }
}

void linearize_and_tally_weights(void)
{
    thearray.len = 0;
    found_sex = -1;
    memset(no_sex_weights,0,(nloci+1)*sizeof(uint32));
    memset(sex_weights,0,(nloci+1)*sizeof(uint32));
#if defined(STEPWISE) || defined(DIAG)
    thetree.size = 0;
#endif
    plinearize_and_tally_weights(&thetree.root);
    assert(!thetree.root);
}

gsl_rng *rng;

void initialize_rng(void)
{
    uint64 seed;
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
#if MACOSX
    int fd;
    if((fd = open("/dev/urandom",O_RDONLY)) < 0)
    {
        fprintf(stderr,"Failed to open \"/dev/urandom\"");
        exit(1);
    }
    if(read(fd,&seed,sizeof(seed)) < (ssize_t)sizeof(seed))
#else
    if(getrandom(&seed,sizeof(seed),0) < (ssize_t)sizeof(seed))
#endif
    {
        fprintf(stderr,"Failed to properly initialize random number generator");
        exit(1);
    }
#if MACOSX
    if(close(fd) < 0)
    {
        fprintf(stderr,"Failed to close \"/dev/urandom\"");
        exit(1);
    }
#endif
    gsl_rng_set(rng,seed);
}

void make_children(uint64 *scratch1)
{
    uint32 i,j;
    bitstr res;
    gsl_ran_discrete_t *tl,*tl_sex = NULL;

    linearize_and_tally_weights();
    tl = gsl_ran_discrete_preproc(thearray.len,thearray.w);
    if(found_sex > -1)
        tl_sex = gsl_ran_discrete_preproc(thearray.len - found_sex,
                thearray.w + found_sex);

    for(i=0;i<nindiv;i++)
    {
        int k = gsl_ran_discrete(rng,tl);
        if(found_sex < 0 || k < found_sex)
        {
            memcpy(scratch1,thearray.bs[k].bits,sizeof(uint64)*nwords);
            res.bits = scratch1;
        }
        else
        {
            uint64 *dad = thearray.bs[k].bits;
            uint64 *mom = thearray.bs[gsl_ran_discrete(rng,tl_sex)+found_sex].bits;
            /* Where are Dad and Mom different? */
            xor(dad,mom,scratch1);
            /* Where do we get Mom's genes? */
            for(j=0;j<nwords;j++)
            {
                /* gsl_rng_mt19937 delivers only 32 bits per call */
                scratch1[j] &= ~gsl_rng_get(rng);
                scratch1[j] &= ~(gsl_rng_get(rng) << 32);
            }
            xor_me(scratch1,dad);
            res.bits = scratch1;
            /* We may have cleared sex bit, so reset it */
            set_sex_bit(res);
        }
        mutate(res);
        insert(res);
    }
}

static inline double env_envelope(double x)
{
    if(x > nloci || x < 0)
        return 0;
    x /= nloci/2.;
    x -= 1;
    return 1 - pow(x,ENVELOPE_DEGREE);
}

void pick_new_env(void)
{
    double cand = env + gsl_ran_gaussian_ziggurat(rng,shift_size);
    double here = env_envelope(env);
    double there = env_envelope(cand);
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

void alloc_mutation_params(void)
{
    mutation_rates = malloc(sizeof(double)*nalleles*nalleles);
    mutation_contrib = malloc(sizeof(double)*nalleles);
}

struct _mutation_sites
{
    uint32 *where;
    uint32 len;
    uint32 space;
} mutation_sites;

struct _power_table
{
    double *table;
    uint32 len;
} power_table;

void initialize_mutation_sites(void)
{
    mutation_sites.where = malloc(sizeof(uint32)*8);
    mutation_sites.len = 0;
    mutation_sites.space = 8;
#ifdef DIAG
    hwm_mutation_sites_space = mutation_sites.space;
#endif
}

void initialize_power_table(void)
{
    uint32 nloci1 = nloci + sex_change;
    power_table.table = malloc(sizeof(double)*8);
    power_table.len = 8;
    int i;
    for(i = 0 ; i < 8 ; i++)
        power_table.table[i] = pow(1 - mutation_rate,nloci1-i);
#ifdef DIAG
    hwm_power_table_len = 8;
#endif
}

void extend_power_table(void)
{
    uint32 nloci1 = nloci + sex_change;
    power_table.table = realloc(power_table.table,sizeof(double)*(power_table.len+8));
    uint32 i;
    for(i = 0 ; i < 8 && i + power_table.len < nloci1 ; i++)
        power_table.table[i+power_table.len] = pow(1 - mutation_rate,nloci1-i-power_table.len);
    power_table.len += 8;
    if(power_table.len >= nloci1)
        power_table.len = nloci1;
#ifdef DIAG
    hwm_power_table_len = power_table.len;
#endif
}

void widen_mutation_sites(void)
{
    uint32 nloci1 = nloci + sex_change;
    mutation_sites.where = 
        realloc(mutation_sites.where,sizeof(int)*(mutation_sites.space+8));
    mutation_sites.space += 8;
    if(mutation_sites.space > nloci1)
        mutation_sites.space = nloci1;
#ifdef DIAG
    hwm_mutation_sites_space = mutation_sites.space;
#endif
}

void mutate(bitstr bs)
{
   uint32 i,j;
   /* If sex_bit can mutate */
   uint32 nloci1 = nloci + sex_change;

   if(1. - mutation_rate == 1.)
       return;

   /* Mutate every god damn one */
   if(mutation_rate == 1.)
   {
        for(i = 0 ; i < nwords - 1 ; i++)
            bs.bits[i] ^= ~0UL;
        bs.bits[nwords - 1] ^= (1UL << residual) - 1;
        if(sex_change)
            bs.bits[nwords - 1] ^= 1UL << residual;
#ifdef DIAG
        mutation_events += nloci1;
#endif
        return;
   }

   uint32 remaining_sites = nloci1;
   mutation_sites.len = 0;
   for(i = 0 ;i < power_table.len ; i++)
   {
with_bigger_table:
       if(gsl_rng_uniform(rng) < power_table.table[i])
           break;
#ifdef DIAG
       mutation_events++;
#endif
       uint32 where = gsl_rng_uniform_int(rng,remaining_sites--);
       uint32 oldwhere;

       /* Loop until we get list with no repeats, and all 
        * sites are selected with equal probability */

       /* remaining_sites < nloci1 means that we index over *remaining* 
        * sites.  This means, for example, if nloci1 = 8 and we have 6,5,3 in
        * mutations_sites.where so far, then if we choose index 4, 
        * then this should correspond to site 7.  The following nested loop 
        * insures that this happens and we have equal probabilities for remaining
        * sites */
       do
       {
           oldwhere = where;
           for(j = 0 ; j < mutation_sites.len; j++)
           {
               if(where >= mutation_sites.where[j])
               {
                   where++;
                   /* Mark as checked */
                   mutation_sites.where[j] += nloci1;
               }
           }
       } while(oldwhere != where);

       /* Unmark */
       for(j = 0 ; j < mutation_sites.len; j++)
           if(mutation_sites.where[j] >= nloci1)
               mutation_sites.where[j] -= nloci1;

       assert(where < nloci1);
       mutation_sites.where[mutation_sites.len++] = where;
       if(mutation_sites.space < nloci1 && mutation_sites.len >= mutation_sites.space)
           widen_mutation_sites();
   }
   if(i == power_table.len && i < nloci1)
   {
       extend_power_table();
       goto with_bigger_table;
   }
   for(i = 0 ; i < mutation_sites.len ; i++)
   {
       uint32 where = mutation_sites.where[i];
       bs.bits[where/bitspword] ^= (1UL << (where % bitspword));
   }
}

int main(int argc, char *argv[])
{
    if(argc == 2 && !strcmp(argv[1],"--usage"))
    {
        fprintf(stderr,
                "./exec \\\n"
                "   nloci \\\n"
                "   shift_rate \\\n"
                "   shift_size \\\n"
                "   discount \\\n"
                "   nosex \\\n"
                "   sex \\\n"
                "   mutation_rate \\\n"
                "   sex_change \\\n"
                "   [ ngen ]\n");
        exit(0);
    }

    bitstr bs;
    initialize_rng();

    if(argc != 
#if defined(STEPWISE)
            10
#else
            11
#endif
      )
    {
        fprintf(stderr,"Wrong number of arguments\n");
        exit(1);
    }

    nloci = (uint32)strtoul(argv[1],NULL,0);
    INVALID(nloci == 0,nloci,argv[1]);

    nalleles = (uint32)strtoul(argv[2],NULL,0);
    INVALID(nalleles < 2,nalleles,argv[1]);

    shift_rate = strtod(argv[3],NULL);
    INVALID(shift_rate < 0 || shift_rate > 1,shift_rate,argv[2]);

    shift_size = strtod(argv[4],NULL);
    INVALID(shift_size < 0,shift_size,argv[3]);

    discount = strtod(argv[5],NULL);
    INVALID(discount < 0 || discount > 1,discount,argv[4]);

    uint32 nosex = (uint32)strtoul(argv[6],NULL,0);
    uint32 sex = (uint32)strtoul(argv[7],NULL,0);

    alloc_mutation_params();

    parse_rates(argv[8]);
    parse_contrib(argv[9]);

    sex_change = !!strtol(argv[10],NULL,0);
#if !defined(STEPWISE)
    uint32 ngen = (uint32)strtoul(argv[11],NULL,0);
#endif

    uint32 nbits = 0;
    uint32 s = nalleles - 1;
    while(s)
    {
        nbits++;
        s >>= 1;
    }

    /* Extra padding for sex bit */
    nwords = nbits/bitspword + 1;
    residual = nbits % bitspword;

    nindiv = nosex + sex;
    if( (int)nosex < 0 || (int)sex < 0 || nindiv == 0)
    {
        fprintf(stderr,"Invalid value for population size: nosex: %s, sex: %s\n",argv[4],argv[5]);
        exit(1);
    }

    /* initialize structures */
    initialize_array();
    initialize_power_table();
    initialize_mutation_sites();
    uint64 *scratch = malloc(sizeof(uint64)*nwords);
    no_sex_weights = malloc((nloci+1)*sizeof(int));
    sex_weights = malloc((nloci+1)*sizeof(int));

    env = nloci/2;

    bs.bits = malloc(sizeof(uint64)*nwords);
#if !defined(STEPWISE) && !defined(DIAG)
    /* Initialize a population randomly */
    uint32 i,j;
    for(i = 0; i < nindiv ; i++)
    {
        for(j = 0 ; j < nwords ; j++)
            bs.bits[j] = gsl_rng_get(rng) + ( gsl_rng_get(rng) << 32 );
        bs.bits[nwords-1] &= (1UL << residual) - 1;
        if(i < sex)
            set_sex_bit(bs);
        insert(bs);
    }

    for(i = 0 ; i < ngen && (!i || printf("\n")); i++)
    {
        make_children(scratch);

        /* Print weights */
        printf("env: %8.2f\n",env);
        printf("no_sex: 0:%u",no_sex_weights[0]);
        for(j = 1 ; j < nloci + 1 ; j++)
            printf(" %u:%u",j,no_sex_weights[j]);
        printf("\n   sex: 0:%u",sex_weights[0]);
        for(j = 1 ; j < nloci + 1 ; j++)
            printf(" %u:%u",j,sex_weights[j]);

        pick_new_env();
    }
    printf("\n");
    return 0;
#elif defined(STEPWISE)
    assert(nloci <= 63);
    uint64 x;
    char s[1024];
    for(;;)
    {
        if(fgets(s,1024,stdin) == NULL)
            break;
        if(s[0] == '\n') { dumptree(); continue;}
        if(s[1] == '\n')
        {
            if(s[0] == 'd')
            {
                dumptree();
                linearize_and_tally_weights();
                uint32 i;
                for(i = 0; i<thearray.len;i++)
                {
                    printf("0x%016lx ",thearray.bs[i].bits[0]);
                } 
                printf("\n");
                for(i = 0; i<thearray.len;i++)
                {
                    printf("%f ",thearray.w[i]);
                }
                printf("\n");
            }
            if(s[0] == 'c')
            {
                make_children(scratch);
                dump_array();
                printf("Found sex?: %d\n",found_sex);
                dumptree();
            }

            continue;
        }

        char *endptr;
        x = strtoul(s,&endptr,0);
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
        }
                
        dumptree();
        printf("    Tree size: %ld\n",thetree.size);
        printf("    ncache.size = %ld\n",ncache.size);

    }
#elif defined(DIAG)
    uint32 i,j,hwm_treesize = 0;
    for(i=0;i<nindiv;i++)
    {
        for(j=0;j<nwords;j++)
            bs.bits[j] = gsl_rng_get(rng) + (gsl_rng_get(rng) << 32);
        bs.bits[nwords-1] &= (1UL << residual) - 1;
        if(i < sex)
            set_sex_bit(bs);
        insert(bs);
    }

    printf("Start\n");
    dumptree();
    for(i=0;i<ngen;i++)
    {
        make_children(scratch);
        dump_array();
        hwm_treesize = (hwm_treesize >= thetree.size)?hwm_treesize:thetree.size;
    }

    printf("End\n");
    dumptree();
    printf("High water mark tree_size: %u\n",hwm_treesize);
    printf("High water mark tree depth: %u\n",hwm_tree_depth);
    printf("High water mark cache_size: %u\n",hwm_cache);
    printf("High water mark power_table.len: %u\n",hwm_power_table_len);
    printf("High water mark mutation_sites.space: %u\n",hwm_mutation_sites_space);
    printf("Mutation events: %lu\n",mutation_events);
    return 0;
#endif
}
