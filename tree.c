#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<fcntl.h>
#include<unistd.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define bitspword (8*sizeof(uint64))

typedef unsigned long uint64;
typedef unsigned int uint32;

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
double mutation_rate;

uint32 nindiv;
uint32 nwords;
uint32 residual;
int env;

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
    for(i=0,s=bits[nwords-1];s && i < residual ;s >>= 1,i++)
        res += s&1;
    return res;
}

int cmp(bitstr b1, bitstr b2)
{
    int i;
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
    thearray.w[thearray.len] = pow(discount,abs((int)n->bs.weight - env))*n->n;
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
#ifdef DIAG
        if(found_sex > 0)
            assert(check_sex_bit(cursor->bs));
#endif
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
    int fd;
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
    if((fd = open("/dev/urandom",O_RDONLY)) < 0)
    {
        fprintf(stderr,"Failed to open \"/dev/urandom\"");
        exit(1);
    }

    if(read(fd,&seed,sizeof(seed)) < (ssize_t)sizeof(seed))
    {
        fprintf(stderr,"Failed to properly initialize random number generator");
        exit(1);
    }
    if(close(fd) < 0)
    {
        fprintf(stderr,"Failed to close \"/dev/urandom\"");
        exit(1);
    }
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
            xor(dad,mom,scratch1);
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

void shift_env(void)
{
    double x = gsl_rng_uniform(rng);
    if(x < shift_rate/2)
        env--;
    else if(x < shift_rate)
        env++;
    env = (env < 0)?0:(env > (int)nloci)?(int)nloci:env;
}

/* Mutation data structures */
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
    mutation_sites.where = malloc(sizeof(uint32)*10);
    mutation_sites.len = 0;
    mutation_sites.space = 8;
#ifdef DIAG
    hwm_mutation_sites_space = mutation_sites.space;
#endif
}

void initialize_power_table(void)
{
    power_table.table = malloc(sizeof(double)*8);
    power_table.len = 8;
    int i;
    for(i = 0 ; i < 8 ; i++)
        power_table.table[i] = pow(1 - mutation_rate,nloci-i);
#ifdef DIAG
    hwm_power_table_len = 8;
#endif
}

void extend_power_table(void)
{
    power_table.table = realloc(power_table.table,sizeof(double)*(power_table.len+8));
    uint32 i;
    for(i = 0 ; i < 8 && i + power_table.len < nloci ; i++)
        power_table.table[i+power_table.len] = pow(1 - mutation_rate,nloci-i-power_table.len);
    power_table.len += 8;
    if(power_table.len >= nloci)
        power_table.len = nloci;
#ifdef DIAG
    hwm_power_table_len = power_table.len;
#endif
}

void widen_mutation_sites(void)
{
    mutation_sites.where = 
        realloc(mutation_sites.where,sizeof(int)*(mutation_sites.space+8));
    mutation_sites.space += 8;
    if(mutation_sites.space > nloci)
        mutation_sites.space = nloci;
#ifdef DIAG
    hwm_mutation_sites_space = mutation_sites.space;
#endif
}

void mutate(bitstr bs)
{
   uint32 i,j,where,remaining_sites = nloci;

   if(1. - mutation_rate == 1.)
       return;

   /* Mutate every god damn one */
   if(mutation_rate == 1.)
   {
        for(i = 0 ; i < nwords - 1 ; i++)
            bs.bits[i] ^= ~0UL;
        bs.bits[nwords - 1] ^= (1UL << residual) - 1;
#ifdef DIAG
        mutation_events += nloci;
#endif
        return;
   }

   mutation_sites.len = 0;
   for(i = 0 ;i < power_table.len ; i++)
   {
with_bigger_table:
       if(gsl_rng_uniform(rng) < power_table.table[i])
           break;
#ifdef DIAG
       mutation_events++;
#endif
       where = gsl_rng_uniform_int(rng,remaining_sites--);

       for(j = 0 ; j < mutation_sites.len; j++)
       {
           if(where == j)
               where++;
       }
       assert(where < nloci);
       mutation_sites.where[mutation_sites.len++] = where;
       if(mutation_sites.space < nloci && mutation_sites.len >= mutation_sites.space)
           widen_mutation_sites();
   }
   if(i == power_table.len && i < nloci)
   {
       extend_power_table();
       goto with_bigger_table;
   }
   for(i = 0 ; i < mutation_sites.len ; i++)
   {
       int where = mutation_sites.where[i];
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
                "   discount \\\n"
                "   nosex \\\n"
                "   sex \\\n"
                "   mutation_rate \\\n"
                "   [ ngen ]\n");
        exit(0);
    }

    bitstr bs;
    initialize_rng();

    if(argc != 
#if defined(STEPWISE)
            7
#else
            8
#endif
      )
            
    {
        fprintf(stderr,"Wrong number of arguments\n");
        exit(1);
    }
    nloci = (uint32)strtoul(argv[1],NULL,0);
    if((int)nloci <= 0)
    {
        fprintf(stderr,"Invalid value for nloci: %s\n",argv[1]);
        goto out;
    }
    shift_rate = strtod(argv[2],NULL);
    discount = strtod(argv[3],NULL);
    uint32 nosex = (uint32)strtoul(argv[4],NULL,0);
    uint32 sex = (uint32)strtoul(argv[5],NULL,0);
    mutation_rate = strtod(argv[6],NULL);
#if !defined(STEPWISE)
    uint32 ngen = (uint32)strtoul(argv[7],NULL,0);
#endif

    if(mutation_rate < 0 || mutation_rate > 1)
    {
        fprintf(stderr,"Invalid value for mutation rate: %s\n",argv[6]);
        goto out;
    }

    /* Extra padding for sex it */
    nwords = nloci/bitspword + 1;
    residual = nloci % bitspword;

    nindiv = nosex + sex;
    if( (int)nosex < 0 || (int)sex < 0 || nindiv == 0)
    
    {
        fprintf(stderr,"Invalid value for population size: nosex: %s, sex: %s\n",argv[4],argv[5]);
        goto out;
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
        printf("env: %d\n",env);
        printf("no_sex: 0:%u",no_sex_weights[0]);
        for(j = 1 ; j < nloci + 1 ; j++)
            printf(" %u:%u",j,no_sex_weights[j]);
        printf("\n   sex: 0:%u",sex_weights[0]);
        for(j = 1 ; j < nloci + 1 ; j++)
            printf(" %u:%u",j,sex_weights[j]);

        shift_env();
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
                printf("%d\n",found_sex);
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
out:
    exit(1);
}
