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
