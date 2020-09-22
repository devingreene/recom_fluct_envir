#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#include<sys/random.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

typedef unsigned long uint64;
typedef unsigned int uint32;

typedef struct _bitstr
{
    uint64 *bits;
    uint32 weight;
} bitstr ;

/* Globals set by cmdline arguments */
int nloci;
int discount;
int nindiv;
double shift_rate;

int nwords;
int residual;
int env;

uint32 weight(uint64 *bits)
{
    int i,res = 0;
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
    for(i = 0 ; i < nwords ; i++)
    {
        if(b1.bits[i] == b2.bits[i])
            continue;
        return b1.bits[i] < b2.bits[i]?1:-1;
    }
    return 0;
}

void xor(uint64 *p, uint64 *q, uint64 *res)
{
    int i;
    for(i=0;i<nwords;i++)
        res[i] = p[i]^q[i];
}

void xor_me(uint64 *p, uint64 *q)
{
    int i;
    for(i=0;i<nwords;i++)
        p[i] ^= q[i];
}

struct node
{
    bitstr bs;
    int n;
    struct node *left;
    struct node *right;
    /* for cache */
    struct node *next;
};

struct tree
{
    struct node *root;
#if defined(STEPWISE) || defined(DIAG)
    size_t size;
#endif
};

struct tree thetree = {NULL
#if defined(STEPWISE) || defined(DIAG)
    ,0};
#else
};
#endif

struct arrays
{
    bitstr *bs;
    double *w;
    int len;
    int space;
};

struct arrays thearray;

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
    int i;
    for(i=thearray.space;i<2*thearray.space;i++,allbits += nwords)
        thearray.bs[i].bits = allbits;

    thearray.space *= 2;
}

/* Global for location of sexuals */
int found_sex = -1;

void append_array(struct node* n)
{

    if(thearray.len > (7*thearray.space)/8)
        increase_array_space();

    memcpy(thearray.bs[thearray.len].bits,n->bs.bits,sizeof(uint64)*nwords);
    /* TODO */
    /* Necessary? */
    thearray.bs[thearray.len].weight = n->bs.weight;
    thearray.w[thearray.len] = pow(discount,abs(n->bs.weight - env))*n->n;
    if(found_sex < 0) 
        /* Two-fold advantage */
        thearray.w[thearray.len] *= 2;
    thearray.len++;
    assert(thearray.len <= nindiv);
}

struct nodecache
{
    struct node* head;
#if defined(STEPWISE) || defined(DIAG)
    size_t size;
#endif
};

struct nodecache ncache = {NULL
#if defined(STEPWISE) || defined(DIAG)
    ,0} ;
#else
    };
#endif

#if defined(STEPWISE) || defined(DIAG)
size_t hwm_cache = 0;
#endif

void putnode(struct node *n)
{
    n->next = ncache.head;
    ncache.head = n;
#if defined(STEPWISE) || defined(DIAG)
    ncache.size++;
    hwm_cache = (hwm_cache < ncache.size)?ncache.size:hwm_cache;
#endif
}

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
#if defined(STEPWISE) || defined(DIAG)
   thetree.size++;
#endif
}

#if defined(STEPWISE) || defined(DIAG)
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

/* Forward declarations */
int check_sex_bit(bitstr bs);
void set_sex_bit(bitstr bs);

void plinearize(struct node **pcursor)
{
    struct node* cursor;
    while((cursor = *pcursor))
    {
        plinearize(&cursor->left);
        *pcursor = cursor->right;
        if(found_sex < 0 && check_sex_bit(cursor->bs))
            found_sex = thearray.len;
        append_array(cursor);
        putnode(cursor);
    }
}

void linearize(void)
{
    thearray.len = 0;
    found_sex = -1;
#if defined(STEPWISE) || defined(DIAG)
    thetree.size = 0;
#endif
    plinearize(&thetree.root);
    assert(!thetree.root);
}

#if defined(STEPWISE) || defined(DIAG)
void printbf(struct node *node)
{
    int i;
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
    partialdump(thetree.root);
}
#endif

gsl_rng *rng;

void initialize_rng(void)
{
    uint64 seed;
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
    getrandom(&seed,sizeof(seed),0);
    gsl_rng_set(rng,seed);
}

void make_children(uint64 *scratch1)
{
    int i,j;
    bitstr res;
    gsl_ran_discrete_t *tl,*tl_sex = NULL;

    linearize();
    tl = gsl_ran_discrete_preproc(thearray.len,thearray.w);
    if(found_sex > -1)
        tl_sex = gsl_ran_discrete_preproc(thearray.len - found_sex,
                thearray.w + found_sex);

    for(i=0;i<nindiv;i++)
    {
        int k = gsl_ran_discrete(rng,tl);
        if(found_sex < 0 || k < found_sex)
             res = thearray.bs[k];
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
       insert(res);
    }
}

int check_sex_bit(bitstr bs)
{
    return !!(bs.bits[nwords-1] & (1UL << residual));
}

void set_sex_bit(bitstr bs)
{
    bs.bits[nwords-1] |= (1UL << residual);
}

/* Globals */
uint32 *no_sex_weights;
uint32 *sex_weights;

void ptally_weights(struct node **pcursor)
{
    struct node *cursor = *pcursor;
    bitstr bs;
    if(!cursor) return;

    ptally_weights(&(cursor->left));
    bs = cursor->bs;
    if(check_sex_bit(bs))
        sex_weights[cursor->bs.weight] += cursor->n ;
    else
        no_sex_weights[cursor->bs.weight] += cursor->n ;
    ptally_weights(&(cursor->right));
}

void tally_weights()
{
    memset(no_sex_weights,0,(nloci+1)*sizeof(uint32));
    memset(sex_weights,0,(nloci+1)*sizeof(uint32));

    ptally_weights(&thetree.root);
}

void shift_env(void)
{
    double x = gsl_rng_uniform(rng);
    if(x < shift_rate/2)
        env--;
    else if(x < shift_rate)
        env++;
    env = (env < 0)?0:(env > nloci)?nloci:env;
}

int main(int argc, char *argv[])
{
    bitstr bs;
    initialize_rng();

    if(argc != 
#if !defined(STEPWISE) && !defined(DIAG)
            7
#else
            6
#endif
      )
            
    {
        fprintf(stderr,"Wrong number of arguments\n");
        exit(1);
    }
    nloci = (int)strtoul(argv[1],NULL,0);
    if(nloci < 1)
    {
        fprintf(stderr,"Invalid value for nloci\n");
        goto out;
    }
    shift_rate = strtod(argv[2],NULL);
    discount = strtod(argv[3],NULL);
    int nosex = (int)strtoul(argv[4],NULL,0);
    int sex = (int)strtoul(argv[5],NULL,0);

    nwords = nloci/(8*sizeof(uint64)) + 1;
    residual = nloci % (8*sizeof(uint64));
    uint64 *scratch = malloc(sizeof(uint64)*nwords);
    nindiv = nosex + sex;
    if(nindiv < 1)
    {
        fprintf(stderr,"Invalid value for population size\n");
        goto out;
    }

    initialize_array();
    env = nloci/2;

#if !defined(STEPWISE) && !defined(DIAG)
    int ngen = (int)strtoul(argv[6],NULL,0);
    int i,j;
    bs.bits = malloc(sizeof(uint64)*nwords);
    no_sex_weights = malloc((nloci+1)*sizeof(int));
    sex_weights = malloc((nloci+1)*sizeof(int));
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
        tally_weights();
        printf("env: %d\n",env);
        printf("no_sex: 0:%u",no_sex_weights[0]);
        for(j = 1 ; j < nloci + 1 ; j++)
            printf(" %u:%u",j,no_sex_weights[j]);
        printf("\n   sex: 0:%u",sex_weights[0]);
        for(j = 1 ; j < nloci + 1 ; j++)
            printf(" %u:%u",j,sex_weights[j]);

        shift_env();
        make_children(scratch);
    }
    printf("\n");
    return 0;
out:
    exit(1);
#elif defined(STEPWISE)
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
                linearize();
                int i;
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
    int i,j,gens = strtol(argv[1],NULL,0);
    uint64 X[nwords];
    bs.bits = X;
    size_t hwm_treesize = 0;
    for(i=0;i<nindiv;i++)
    {
        for(j=0;j<nwords;j++)
            X[j] = gsl_rng_get(rng) + (gsl_rng_get(rng) << 32);
        X[nwords-1] &= (1 << (nloci + 1) % (8*sizeof(uint64))) - 1;
        insert(bs);
    }

    printf("Start\n");
    dumptree();
    for(i=0;i<gens;i++)
    {
        make_children(scratch);
        hwm_treesize = (hwm_treesize >= thetree.size)?hwm_treesize:thetree.size;
    }

    printf("End\n");
    dumptree();
    printf("High water mark tree_size: %ld\n",hwm_treesize);
    printf("High water mark cache_size: %ld\n",hwm_cache);
#endif
}
