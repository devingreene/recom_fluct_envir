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

#define nwords (nloci/(8*sizeof(uint64)) + 1)
#define residual ((nloci % (8*sizeof(uint64))) + 1)

uint32 nloci;
int env;
uint32 discount;
uint32 nindiv; 
double shift_rate;

uint32 weight(uint64 *bits)
{
    int i,res = 0;
    uint64 s = 0;
    for(i = 0; i < nwords - 1; i++)
        for(s=bits[i];s;s >>= 1)
            res += s&1;
    for(i=residual - 1,s=bits[nwords-1];s && i > 0;s >>= 1,i--)
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
    size_t size;
};

struct tree thetree = {NULL, 0};

struct arrays
{
    bitstr *bs;
    double *w;
    int len;
    int space;
};

struct arrays array;

void initialize_array(void)
{
    array.bs = malloc(sizeof(*array.bs)*0x1000);
    int i;
    uint64 *allbits = malloc(sizeof(uint64)*nwords*0x1000);
    for(i=0;i<0x1000;i++,allbits += nwords)
        array.bs[i].bits = allbits;
    array.w    = malloc(sizeof(double)*0x1000);
    array.len  = 0;
    array.space = 0x1000;
}

void widen_array(void)
{
    array.bs = realloc(array.bs,sizeof(*array.bs)*2*array.space);
    array.w  = realloc(array.w,sizeof(double)*2*array.space);
    int i ;
    uint64 *allbits = malloc(sizeof(uint64)*nwords*array.space);
    for(i=array.space;i<2*array.space;i++,allbits += nwords)
        array.bs[i].bits = allbits;
    array.space *= 2;
}

int found_sex;
void append_array(struct node* n)
{
    if(array.len > 7*array.space/8)
    {
        widen_array();
    }

    memcpy(array.bs[array.len].bits,n->bs.bits,sizeof(uint64)*nwords);
    array.bs[array.len].weight = n->bs.weight;
    array.w[array.len++] = pow(discount,abs(n->bs.weight - env))*n->n;
    if(found_sex == -1) 
        array.w[array.len-1] *= 2;
}

struct nodecache
{
    struct node* head;
    size_t size;
};

struct nodecache ncache = {NULL,0} ;
size_t hwm_cache = 0;

void putnode(struct node *n)
{
    n->next = ncache.head;
    ncache.head = n;
    ncache.size++;
    hwm_cache = (hwm_cache < ncache.size)?ncache.size:hwm_cache;
}

struct node *getnode(bitstr bs)
{
    struct node *n;
    if(ncache.head)
    {
        n = ncache.head;
        ncache.head = n->next;
        assert(ncache.size);
        ncache.size--;
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
   thetree.size++;
}

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

int check_sex_bit(bitstr bs);
void set_sex_bit(bitstr bs);

void plinearize(struct node **pcursor)
{
    struct node* cursor;
    while((cursor = *pcursor))
    {
        plinearize(&cursor->left);
        *pcursor = cursor->right;
        if(found_sex == -1 && check_sex_bit(cursor->bs))
            found_sex = array.len;
        append_array(cursor);

        putnode(cursor);
    }
}

void linearize(void)
{
    array.len = 0;
    found_sex = -1;
    thetree.size = 0;
    plinearize(&thetree.root);
    assert(!thetree.root);
}

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
    tl = gsl_ran_discrete_preproc(array.len,array.w);
    if(found_sex > -1)
        tl_sex = gsl_ran_discrete_preproc(array.len - found_sex,
                array.w + found_sex);

    for(i=0;i<nindiv;i++)
    {
        int k = gsl_ran_discrete(rng,tl);
        if(found_sex == -1 || k < found_sex)
             res = array.bs[k];
        else
        {
            uint64 *dad = array.bs[k].bits;
            uint64 *mom = array.bs[gsl_ran_discrete(rng,tl_sex)].bits;
            xor(dad,mom,scratch1);
            for(j=0;j<nwords;j++)
            {
                scratch1[j] &= ~gsl_rng_get(rng);
                scratch1[j] &= ~(gsl_rng_get(rng) << 32);
            }
            xor_me(scratch1,dad);
            res.bits = scratch1;
            set_sex_bit(res);
        }
       insert(res);
    }
}

int check_sex_bit(bitstr bs)
{
    return bs.bits[nwords-1] & (1UL << (residual - 1));
}

void set_sex_bit(bitstr bs)
{
    bs.bits[nwords-1] |= (1UL << (residual - 1));
}

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
    initialize_array();
    initialize_rng();
    uint64 scratch[nwords];

    if(argc != 7)
    {
        fprintf(stderr,"Wrong number of arguments\n");
        exit(1);
    }
    nloci = (uint32)strtoul(argv[1],NULL,0);
    shift_rate = strtod(argv[2],NULL);
    discount = strtod(argv[3],NULL);
    uint32 nosex = (uint32)strtoul(argv[4],NULL,0);
    uint32 sex = (uint32)strtoul(argv[5],NULL,0);
    uint32 ngen = (uint32)strtoul(argv[6],NULL,0);

    nindiv = nosex + sex;
    env = nloci/2;

#if defined(STEPWISE)
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
#else

    uint32 i,j;
    uint64 X[nwords];
    bs.bits = X;
    no_sex_weights = malloc((nloci+1)*sizeof(uint32));
    sex_weights = malloc((nloci+1)*sizeof(uint32));
    for(i = 0; i < nindiv ; i++)
    {
        for(j = 0 ; j < nwords ; j++)
            X[j] = gsl_rng_get(rng) + ( gsl_rng_get(rng) << 32 );
        X[nwords-1] &= (1UL << (residual - 1)) - 1;
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
#endif
}
