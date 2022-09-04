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

#define min(x,y) ((x) > (y)?y:x)
#define INDEX(i,j) (nalleles*(i) + (j))

typedef unsigned long uint64;
typedef unsigned int uint32;

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

/* Calculated paramters */
uint32 nindiv;
uint32 nwords;
uint32 residual;
uint32 nalleles;
uint32 allele_size;
uint64 allele_mask;
uint32 maximum_weight;

/* Environmental paramter */
double env;

/* Array for storing parental generation */
#define INITIAL_ARRAY_SIZE (0x1000)
struct _array
{
    bitstr *bs;
    double *w;
    uint32 len;
    uint32 space;
} array;

/* Initializers */
void initialize_array(void);
void initialize_rng(void);
void initialize_mutation_parameters(char *argv[]);
void initialize_sex_weights(void);
void setBitParameters(void);

/* Bitstring operations */
int check_sex_bit(bitstr bs);
void set_sex_bit(bitstr bs);
void insert(bitstr bs);
int cmp(bitstr b1, bitstr b2);

/* General functions */
uint32 weight(uint64 *bits);
void linearize_and_tally_weights(void);
void make_children(uint64 *scratch1,uint32 *choices, uint32 choices_ints);
void pick_new_env(void);
void mutate(bitstr bs);

/* GSL PRNG */
gsl_rng *rng;

/* Arrays for output */
uint32 *no_sex_weights;
uint32 *sex_weights;
