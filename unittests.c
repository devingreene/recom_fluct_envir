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

struct _array
{
    bitstr *bs;
    double *w;
    uint32 len;
    uint32 space;
} array;

uint32 nloci;
double discount;
double shift_rate;
double shift_size;
double *mutation_rate;
uint32 *mutation_contrib;
gsl_ran_discrete_t **mutant_tables;
double sex_mutation_rate;
uint32 sex_change;

struct node;

/* Bitstring operations */
int check_sex_bit(bitstr bs);
void set_sex_bit(bitstr bs);
void insert(bitstr bs) ;
int cmp(bitstr b1, bitstr b2);

void initialize_array(void);
void initialize_rng(void);
double weight(uint64 *bits);
void linearize_and_tally_weights(void);
void make_children(uint64 *scratch1,uint32 *choices, uint32 choices_ints);
void pick_new_env(void);
void initialize_mutation_parameters(char *argv[]);
void setBitParameters(void);
void initialize_sex_weights(void);

uint32 nindiv;
uint32 nwords;
uint32 residual;
uint32 nalleles;
uint32 allele_size;
uint64 allele_mask;

double maximum_weight;
double env;

gsl_rng *rng;

uint32 *no_sex_weights;
uint32 *sex_weights;

void parse_rates(char *s);
void parse_contrib(char *s);

uint32 *traits;
uint32 ntraits;

int main()
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

    /* Does tree build and collapse correctly? */

    nindiv = 19;

    typedef uint64 gtype[nwords];
    gtype *gtypes = malloc(sizeof(gtype)*nindiv);

    uint32 perm[] = { 9,18,17,4,16,1,8,14,15,2,10,12,7,11,0,3,5,6,13 };

    uint32 mc[] = { 0,1,2,3,4 };
    mutation_contrib = malloc(sizeof(uint32)*nalleles);
    memcpy(mutation_contrib,mc,sizeof(uint32)*nalleles);

    uint32 static_traits[] = { 13, 26, 39 };
    traits = &static_traits[0];
    ntraits = 4;

    bitstr *indvs = malloc(sizeof(bitstr)*(nindiv));

    uint i;
    for(i = 0; i < nindiv - 1; i++)
    {
        gtypes[i][0] = ((uint64)(i/nalleles) << allele_size) + (uint64)(i % nalleles);
        gtypes[i][1] = 0x0404040404040404;
        gtypes[i][2] = 0x4040404040404040;
        gtypes[i][3] = (i % 2)?0x10003:0x3;
    }

    // 19th genotype is repeat of 18th
    gtypes[i][0] = gtypes[i-1][0];
    gtypes[i][1] = gtypes[i-1][1];
    gtypes[i][2] = gtypes[i-1][2];
    gtypes[i][3] = gtypes[i-1][3];

    assert(weight(gtypes[1]) == sqrt(1.0*1.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0));

    for(i = 0; i < nindiv; i++)
        indvs[i].bits = gtypes[i];

    for(i = 0 ; i < nindiv ; i++)
        insert(indvs[perm[i]]);

    initialize_array();
    maximum_weight = (double)( nalleles - 1 )*nloci;
    initialize_sex_weights();
    discount = 1.01;env = 70.0;
    linearize_and_tally_weights();
#define A(n,sex) (pow(discount, fabs((n) - env))*((sex)?1:2))
// sqrt(1.0*1.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0)
    assert( array.w[0] == A(sqrt(0.0*0.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[1] == A(sqrt(2.0*2.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[2] == A(sqrt(4.0*4.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[3] == A(sqrt(2.0*2.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[4] == A(sqrt(4.0*4.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[5] == A(sqrt(2.0*2.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[6] == A(sqrt(4.0*4.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[7] == A(sqrt(6.0*6.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[8] == A(sqrt(4.0*4.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),0) );
    assert( array.w[9] == A(sqrt(1.0*1.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[10] == A(sqrt(3.0*3.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[11] == A(sqrt(1.0*1.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[12] == A(sqrt(3.0*3.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[13] == A(sqrt(5.0*5.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[14] == A(sqrt(3.0*3.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[15] == A(sqrt(5.0*5.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[16] == A(sqrt(3.0*3.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );
    assert( array.w[17] == 2*A(sqrt(5.0*5.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0),1) );

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
    maximum_weight = (double)nloci;
    shift_size = 1.0;
    env = maximum_weight/2.0;
    assert(env >= 0);

    /* Fixed rng */
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
    gsl_rng_set(rng,0xdebb1eU);

    uint32 *freq = calloc((nloci + 1),sizeof(int));
    for(i = 10000 ; i; i--)
    {
        freq[(uint32)env]++;
        pick_new_env();
    }

    uint32 mode = 0;
    for(i = 0; i <= nloci ; i++)
        mode = (freq[i] > mode)?freq[i]:mode;

    mode /= 15;

    uint j,k;
    for(j = mode; j; j--)
    {
        for(k=0 ; k <= nloci ; k++)
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

    env = maximum_weight/2.0;
    shift_size = 7.5; // Large SD
    /* Ensure that `env' stays in range
     * ( assert statement in function pick_new_env )
     * */
    for(i = 1000000; i; i--)
        pick_new_env();
    
    /* TODO
     * Recombination, among other things
     * */
}
