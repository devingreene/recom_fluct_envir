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
    double weight;
} bitstr ;

/* Globals set by cmdline arguments */
extern uint32 nloci;
extern double discount;
extern double shift_rate;
extern double shift_size;
extern double *mutation_rate;
extern uint32 *mutation_contrib;
extern gsl_ran_discrete_t **mutant_tables;
extern double sex_mutation_rate;
extern uint32 sex_change;

/* Calculated paramters */
extern uint32 nindiv;
extern uint32 nwords;
extern uint32 residual;
extern uint32 nalleles;
extern uint32 allele_size;
extern uint64 allele_mask;
extern double maximum_weight;

/* Environmental paramter */
extern double env;

#define INITIAL_ARRAY_SIZE (0x1000)
extern struct _array
{
    bitstr *bs;
    double *w;
    uint32 len;
    uint32 space;
} array;

/* Initializers */
extern void initialize_array(void);
extern void initialize_rng(void);
extern void initialize_mutation_parameters(char *argv[]);
extern void initialize_sex_weights(void);
extern void setBitParameters(void);
extern void initialize_population(uint32 sex);

/* Bitstring operations */
extern int check_sex_bit(bitstr bs);
extern void set_sex_bit(bitstr bs);
extern void insert(bitstr bs);
extern int cmp(bitstr b1, bitstr b2);

/* General functions */
extern double weight(uint64 *bits);
extern void linearize_and_tally_weights(void);
extern void make_children(uint64 *scratch1,uint32 *choices, uint32 choices_ints);
extern void pick_new_env(void);
extern void mutate(bitstr bs);

/* GSL PRNG */
extern gsl_rng *rng;

/* Arrays for output */
extern uint32 *no_sex_weights;
extern uint32 *sex_weights;

extern uint32 *traits;
extern uint32 ntraits;
