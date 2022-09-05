#include "recom_fluct_envir_many.h"

/* Marks index of first sexual */
int found_sex = -1;

void parse_rates(char *s);
void parse_contrib(char *s);

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
    array.w[array.len] = pow(discount,fabs(n->bs.weight - env))*n->n;
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
        assert(( res.bits[nwords-1] & (uint64)(-1UL) << (residual + 1) )  == 0);

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

/* Mutation functions */

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

void initialize_population(uint32 sex)
{
    bitstr bs;
    bs.bits = malloc(sizeof(uint64)*nwords);

    /* Initialize a population randomly */
    uint32 i,j,k;
    size_t rand;
    double *uniform = malloc(sizeof(double)*nalleles);

    for(i = 0 ; i < nalleles ; i++) uniform[i] = 1.;
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
}
