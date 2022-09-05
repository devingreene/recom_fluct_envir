#include<ctype.h>
#include "recom_fluct_envir_many.h"

/*
 * Engine for parsing mutation rates and fitness
 * contributions
 * */

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

#define GO() \
{ \
    while(isspace(*s2)) s2++; \
}

#define ISDELIM(c) ( (c) == ')' || (c) == '(' || (c) == '[' || (c) == ']' )
#define ISLEFT(c) ( (c) == '(' || (c) == '[' )
#define ISRIGHT(c) ( (c) == ')' || (c) == ']' )

void parse_rates(char *s)
{
    char *s2;
    double rate = strtod(s,&s2);
    GO();
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
        GO();
        if(ISLEFT(*s2))
        {
            if(depth > 0)
                GRAPH("Too deep\n");
            depth++;
            if(++row >= (int)nalleles)
                GRAPH("Too many rows\n");
            col++;
            s2++;
        }

        else if(ISRIGHT(*s2))
        {
            if(depth == 0)
                GRAPH("Too shallow\n");
            if(col < (int)nalleles )
                GRAPH("Not enough elements in set\n");
            depth--;
            col = -1;
            s2++;
        }

        else if(*s2)
        {
            if(depth == 0)
                GRAPH("Token out of set\n");
            if(row == col)
                /* If diagonal element, skip word */
                while(*s2 && !isspace(*s2) && !ISDELIM(*s2)) s2++;
            else
                mutation_rate[INDEX(row,col)] = strtod(s2,&s2);
            if(!isspace(*s2) && !ISRIGHT(*s2))
                GRAPH("Unexpected character or end\n");
            col++;
            if(col > (int)nalleles)
                GRAPH("Too many columns\n");
        }
    }
    if(depth > 0)
        GRAPH("Unexpected end of string\n");
    if(row < (int)nalleles - 1)
        GRAPH("Too few rows\n");
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
        INVALID((mutation_rate[i] < 0),"Inferred negative probability values\n");
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
    char *ep;
    uint32 depth = 0;
    uint32 col = -1;
    while(*s2)
    {
        GO();
        if(ISLEFT(*s2))
        {
            if(depth > 0)
                GRAPH("Too deep\n");
            depth++;
            col++;
            s2++;
            GO();
            goto token;
        }

        else if(ISRIGHT(*s2))
        {
            if(depth <= 0)
                GRAPH("Too shallow\n");
            depth--;
            col = -1;
            s2++;
            GO();
        }

        else if(depth == 0 && *s2)
            GRAPH("Unexpected character\n")

        else if(depth == 1)
        {
token:
            if(col >= 1)
                GRAPH("Too many entries\n");
            uint32 n = strtol(s2,&ep,0);
            if(ep == s2)
                GRAPH("Missing or invalid entry\n");
            s2 = ep;
            if(!isspace(*s2) && !ISRIGHT(*s2))
                GRAPH("%s\n","Unexpected character or end");
            if(n >= nalleles)
                GRAPH("No such allele\n");
            if(mutation_contrib[n] != UINT_MAX)
                GRAPH("Illegal: %u already filled in\n",n);
            col++;
            GO();
            mutation_contrib[n] = strtol(s2,&ep,0);
            if(ep == s2)
                GRAPH("Missing or invalid entry\n");
            if(*s2 == '-')
                GRAPH("Negative contributions not accepted\n");
            s2 = ep;
            if(!isspace(*s2) && !ISRIGHT(*s2))
                GRAPH("%s\n","Unexpected character or end");
            col++;
        }
    }
    if(depth != 0)
        GRAPH("Unexpected end of string\n");

    for(i = 0 ; i < nalleles ; i++)
        if(mutation_contrib[i] == UINT_MAX)
            PARSE_ERROR("Allele index %u not filled in\n",i);
}
