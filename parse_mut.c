#include<stdlib.h>
#include<stdio.h>
#include<ctype.h>
#include<math.h>
#include<limits.h>
#include "this.h"

#define PARSE_ERROR(fmt,...)  \
{ \
    fprintf(stderr,fmt,##__VA_ARGS__); \
    exit(1); \
}

#define GRAPH(fmt,...) \
{ \
    fprintf(stderr,"%s\n",s); \
    int _i; \
    for(_i = 0; _i < (s2 - s) ; _i++) \
        fprintf(stderr," "); \
    fprintf(stderr,"^\n"); \
    PARSE_ERROR(fmt,##__VA_ARGS__) \
}

#define INVALID(fmt,...) \
{ \
    fprintf(stderr,fmt,##__VA_ARGS__); \
    exit(1); \
}


#define ISDELIM(c) ( (c) == ')' || (c) == '(' || (c) == '[' || (c) == ']' )
#define ISLEFT(c) ( (c) == '(' || (c) == '[' )
#define ISRIGHT(c) ( (c) == ')' || (c) == ']' )

// #define extern
extern double *mutation_rate;
extern uint32 *mutation_contrib;
extern uint32 nalleles;

void parse_rates(char *s)
{
    char *s2;
    double rate = strtod(s,&s2);
    while(isspace(*s2)) s2++;
    /* Only a number?  Use it as the rate */
    if(!*s2)
    {
        uint32 i,j;
        for(i = 0 ; i < nalleles ; i++)
        {
            for(j = 0; j < nalleles ; j++)
                if(i != j)
                    mutation_rate[INDEX(i,j)] = rate/(nalleles - 1);
            mutation_rate[INDEX(i,i)] = 1 - rate;
        }
        return;
    }

    /* This must be a matrix, so let's start over */
    int depth = 0;
    int row = -1;
    int col = -1;
    s2 = s;

    while(*s2)
    {
        while(isspace(*s2)) s2++;
        if(ISLEFT(*s2))
        {
            if(depth > 0)
                GRAPH("%s\n","Too deep");
            depth++;
            if(++row >= (int)nalleles)
                GRAPH("%s\n","Too many rows");
            col++;
            s2++;
        }

        else if(ISRIGHT(*s2))
        {
            if(col < (int)nalleles )
                GRAPH("%s\n","Not enough elements in set");
            if(depth == 0)
                GRAPH("%s\n","Too shallow");
            depth--;
            col = -1;
            s2++;
        }

        else if(*s2)
        {
            if(depth == 0)
                GRAPH("%s\n","Token out of set");
            if(row == col)
                /* If diagonal element, skip word */
                while(*s2 && !isspace(*s2) && !ISDELIM(*s2)) s2++;
            else
                mutation_rate[INDEX(row,col)] = strtod(s2,&s2);
            if(!isspace(*s2) && !ISRIGHT(*s2))
                GRAPH("%s\n","Unexpected character or end");
            col++;
            if(col > (int)nalleles)
                GRAPH("%s\n","Too many columns");
        }
    }
    if(depth > 0)
        GRAPH("%s\n","Unexpected end of string");
    if(row < (int)nalleles - 1)
        GRAPH("%s\n","Too few rows");
    /* Fill in the diagonal elements */
    uint32 i,j;
    double sum;
    for(i = 0 ; i < nalleles ; i++)
    {
        sum = 0;
        for(j = 0 ; j < nalleles ; j++)
            if(i != j)
                sum += mutation_rate[INDEX(i,j)];
        mutation_rate[INDEX(i,i)] = 1 - sum;
    }

    for(i = 0 ; i < nalleles*nalleles ; i++)
        if(mutation_rate[i] < 0)
            INVALID("Inferred negative probability values\n");
}

void parse_contrib(char *s)
{
    uint32 i;
    for(i = 0 ; i < nalleles ; i++)
        mutation_contrib[i] = UINT_MAX;
    if(*s == 0)
    {
        if(nalleles > 2)
            PARSE_ERROR("Three or more alleles must have a fitness "
                    "contribution table\n")
        else
        {
            mutation_contrib[0] = 0;
            mutation_contrib[1] = 1;
        }
    }

    char *s2 = s;
    uint32 depth = 0;
    uint32 col = -1;
    while(*s2)
    {
        while(isspace(*s2)) s2++;
        if(ISLEFT(*s2))
        {
            if(depth > 0)
                GRAPH("%s\n","Too deep");
            depth++;
            col++;
            s2++;
        }

        else if(ISRIGHT(*s2))
        {
            if(col >= 2)
                GRAPH("%s\n","Too many columns");
            depth--;
            col = -1;
            s2++;
        }

        else
        {
            if(col == 0)
            {
                uint32 n = strtol(s2,&s2,0);
                if(!isspace(*s2) && !ISRIGHT(*s2))
                    GRAPH("%s\n","Unexpected character or end");
                if(n >= nalleles)
                    GRAPH("No such allele\n");
                if(mutation_contrib[n] != UINT_MAX)
                    GRAPH("Illegal: %u already filled in\n",n);
                mutation_contrib[n] = strtol(s2,&s2,0);
                if(!isspace(*s2) && !ISRIGHT(*s2))
                    GRAPH("%s\n","Unexpected character or end");
            }
            else
                GRAPH("Unexpected character\n");
        }
    }
    for(i = 0 ; i < nalleles ; i++)
        if(mutation_contrib[i] == UINT_MAX)
            PARSE_ERROR("Allele index %u not filled in\n",i);
}

/*
int main(int argc, char *argv[])
{
    if(argc != 4) { fprintf(stderr,"Fuck off\n");exit(1);}
    nalleles = strtol(argv[1],NULL,0);
    mutation_rate = malloc(sizeof(double)*nalleles*nalleles);
    parse_rates(argv[2]);
    uint32 i,j;
    for(i=0 ; i < nalleles && (!i || printf("\n")); i++)
        for(j=0 ; j < nalleles ; j++)
            printf("%f ",mutation_rate[INDEX(i,j)]);
    printf("\n");

    mutation_contrib = malloc(sizeof(double)*nalleles);
    parse_contrib(argv[3]);
    printf("\n");
    for(i = 0 ; i < nalleles ; i++)
        printf("%f ",mutation_contrib[i]);
    printf("\n");
}
*/
