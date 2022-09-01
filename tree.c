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
    /* Bitstrings are ordered word-wise little endian.  This
     * is so that the sex bit sorts the array into no_sex
     * genotypes followed by sex genotypes.  */

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
} ncache = {NULL};

struct _tree
{
    struct node *root;
} tree = {NULL};

void putnode(struct node *n)
{
    n->next = ncache.head;
    ncache.head = n;
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
   struct node **pcursor = &tree.root;
   struct node *cursor;

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
   }

   struct node *new = getnode(bs);

   *pcursor = new;
}

#define INITIAL_ARRAY_SIZE (0x1000)
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
    array.bs = malloc(sizeof(*array.bs)*INITIAL_ARRAY_SIZE);
    /* Bit strings in one call to malloc */
    uint64 *allbits = malloc(sizeof(uint64)*nwords*INITIAL_ARRAY_SIZE);

    int i;
    for(i=0;i<INITIAL_ARRAY_SIZE;i++, /* ptr arith */ allbits += nwords)
        array.bs[i].bits = allbits;
    array.w    = malloc(sizeof(double)*INITIAL_ARRAY_SIZE);
    array.len  = 0;
    array.space = INITIAL_ARRAY_SIZE;
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
    }
}

void linearize_and_tally_weights(void)
{
    array.len = 0;
    found_sex = -1;
    memset(no_sex_weights,0,(maximum_weight + 1)*sizeof(int));
    memset(sex_weights,0,(maximum_weight + 1)*sizeof(int));
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

/* Most of the action happens here.  We tear down our tree
 * into an array and build a new generation from it. */
void make_children(uint64 *scratch1,uint32 *choices, uint32 choices_ints)
{
    uint32 i,j;
    bitstr res;
    gsl_ran_discrete_t *tl,*tl_sex = NULL;

    linearize_and_tally_weights();
    tl = gsl_ran_discrete_preproc(array.len,array.w);
    if(found_sex > -1)
        tl_sex = gsl_ran_discrete_preproc(array.len - found_sex,
                /* Ptr arith */ array.w + found_sex);

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
            /* Recombination */

            /* Initialize to zero, then fill in with mom and
             * dad's alleles */
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
            /* Set the sex bit */
            set_sex_bit(res);
        }
        /* Assert that double checks that bits more significant
         * than sex bit are zero */
        assert(( res.bits[nwords-1] & ( (uint64)(-1UL) << (residual + 1) ) )  == 0);

        mutate(res);
        insert(res);
    }
    gsl_ran_discrete_free(tl);
    gsl_ran_discrete_free(tl_sex);

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
        /* Pointer arithmetic in 2nd argument */
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
}

void setBitParameters(void)
{
    /* allele fit in windows of 2,4,8,16,32 or 64 bits */
    assert(nalleles >= 2);
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

}

void initialize_sex_weights(void)
{
    no_sex_weights = malloc((maximum_weight + 1)*sizeof(uint32));
    sex_weights = malloc((maximum_weight + 1)*sizeof(uint32));
}

void unittests(void);

int main(int argc, char *argv[])
{
    if(argc ==2 && !strcmp(argv[1],"unittests"))
    {
        unittests();
        printf("\nAll tests passed\n");
        exit(0);
    }
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
                "   ngen\n");
        exit(0);
    }

    if(argc != 13)
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

    uint32 ngen = (uint32)strtoul(argv[12],NULL,0);

    setBitParameters();

    nindiv = nosex + sex;
    INVALID(nindiv == 0,
            "Invalid value for population size: nosex: %s, sex: %s\n",argv[6],argv[7]);

    initialize_rng();
    initialize_array();

    /* Passed arrays for make_children */
    uint64 *scratch = malloc(sizeof(uint64)*nwords);
    /* Used in recombination: One bit per locus, rounded up to nearest whole
     * number of ints - (nloci + bitspnt -1)/bitspint = # of ints needed to fit
     * nloci bits */
    uint32 *choices = malloc(sizeof(uint32)*((nloci + bitspint - 1)/bitspint));
    uint32 choices_ints = (nloci + bitspint - 1)/bitspint;

    initialize_sex_weights();

    env = maximum_weight/2;

    bs.bits = malloc(sizeof(uint64)*nwords);

    /* Initialize a population randomly */
    uint32 i,j,k;
    size_t rand;
    double *uniform = malloc(sizeof(double)*nalleles);

    for(i = 0 ; i < nalleles ; i++)
        uniform[i] = 1.;
    gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(nalleles,uniform);

    for(i = 0; i < nindiv ; i++)
    {
        /* Initialize to zero so that 
         *  - we can just xor alleles
         *  - sex bit is MSB */
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
        insert(bs);
    }
    gsl_ran_discrete_free(table);
    free(uniform);
    free(bs.bits);

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
}

void unittests(void)
{
    nloci = 52;
    nalleles = 5;

    /* Check basic bit string functions */
    setBitParameters();
    assert( nwords == 4 );
    assert( allele_size == 4 );
    assert( allele_mask == 0x0f );
    assert( residual == 16 );

    uint64 *bits = calloc(sizeof(long)*nwords*2,1);
    uint64 *bits2 = bits + nwords;

    bits[3] = 0x00ff;bits[0] = (uint64)-1;
    bits2[3] = 0x100ff;bits2[0] = (uint64)-1 - ( 0xffUL << 56 );

    bitstr *bs = calloc(sizeof(bitstr),1);
    bitstr *bs2 = calloc(sizeof(bitstr),1);

    bs->bits = bits;
    bs2->bits = bits2;
    assert( cmp(*bs,*bs2) == 1 );
    assert( cmp(*bs, *bs) == 0 );

    assert(check_sex_bit(*bs) == 0);
    assert(check_sex_bit(*bs2) == 1);

    set_sex_bit(*bs);
    assert(cmp(*bs,*bs2) == -1);

    free(bits);free(bs);free(bs2);

    /* Does tree collapse correctly? */

    nindiv = 18;

    typedef uint64 gtype[nwords];
    gtype *gtypes = malloc(sizeof(gtype)*nindiv);

    uint32 perm[] = { 9,17,4,16,1,8,14,15,2,10,12,7,11,0,3,5,6,13 };

    uint32 mc[] = { 0,1,2,3,4 };
    mutation_contrib = malloc(sizeof(uint32)*nalleles);
    memcpy(mutation_contrib,mc,sizeof(uint32)*nalleles);

    bitstr *indvs = malloc(sizeof(bitstr)*(nindiv));

    uint i;
    for(i = 0; i < nindiv; i++)
    {
        gtypes[i][0] = ((uint64)(i/nalleles) << allele_size) + (uint64)(i % nalleles);
        gtypes[i][1] = 0x0404040404040404;
        gtypes[i][2] = 0x4040404040404040;
        gtypes[i][3] = (i % 2)?0x10003:0x3;
    }

    assert(weight(gtypes[1]) == 1 + 32 + 32 + 3);

    for(i = 0; i < nindiv; i++)
        indvs[i].bits = gtypes[i];

    for(i = 0 ; i < nindiv ; i++)
        insert(indvs[perm[i]]);

    initialize_array();
    maximum_weight = ( nalleles - 1 )*nloci;
    initialize_sex_weights();
    discount = 1.01;env = 70.0;
    linearize_and_tally_weights();
#define A(n,sex) (pow(discount, fabs((n) - env))*((sex)?1:2))

    assert( array.w[0] == A(67,0) );
    assert( array.w[1] == A(69,0) );
    assert( array.w[2] == A(71,0) );
    assert( array.w[3] == A(69,0) );
    assert( array.w[4] == A(71,0) );
    assert( array.w[5] == A(69,0) );
    assert( array.w[6] == A(71,0) );
    assert( array.w[7] == A(73,0) );
    assert( array.w[8] == A(71,0) );
    assert( array.w[9] == A(68,1) );
    assert( array.w[10] == A(70,1) );
    assert( array.w[11] == A(68,1) );
    assert( array.w[12] == A(70,1) );
    assert( array.w[13] == A(72,1) );
    assert( array.w[14] == A(70,1) );
    assert( array.w[15] == A(72,1) );
    assert( array.w[16] == A(70,1) );
    assert( array.w[17] == A(72,1) );

    free(gtypes);
    free(indvs);

    /* parsers */
    nloci = 1;
    nalleles = 3;
    setBitParameters();
    mutation_rate = malloc(sizeof(double)*nalleles*nalleles);
    char *s = "(0.125 0.625 0.125)"
              "(0.125 0.5 0.25)"
              "(0.25 0.75 0.25)";
    parse_rates(s);
    assert( mutation_rate[0] ==  0.25 );
    assert( mutation_rate[4] == 0.625 );
    assert( mutation_rate[7] == 0.75  );

    char *t = "(0 1) (2 2) (1 3)";
    free(mutation_contrib);
    mutation_contrib = malloc(sizeof(uint32)*3);
    parse_contrib(t);
    assert( mutation_contrib[0] == 1 );
    assert( mutation_contrib[1] == 3 );
    assert( mutation_contrib[2] == 2 );

    /* Visual display of distribution of env values */
    nloci = 80;
    nalleles = 2;
    setBitParameters();
    free(mutation_contrib);
    mutation_contrib = malloc(sizeof(uint32)*2);
    parse_contrib("(0 0) (1 1)");
    maximum_weight = nloci;
    shift_size = 1.0;
    env = maximum_weight/2.0;
    initialize_rng();

    int *freq = calloc((nloci + 1),sizeof(int));
    for(i = 0 ; i < 10000 ; i++)
    {
        freq[(uint32)env]++;
        pick_new_env();
    }

    int mode = 0;
    for(i = 0; i <= nloci ; i++)
        mode = (freq[i] > mode)?freq[i]:mode;

    mode /= 15;

    int j,k;
    for(j = mode; j >= 0; j--)
    {
        for(k=0 ; k <= (int)nloci ; k++)
        {
            if(freq[k]/15 >= j)
                printf("*");
            else
                printf(" ");
        }
        printf("\n");
    }
    printf("Visual display of distribution of env values\n"
            "(10000 samples, 80 loci, vertical axis scaled)\n");
    free(freq);
    free(mutation_contrib);
    free(mutation_rate);
}
