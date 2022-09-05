#include "this.h"

int main(int argc, char *argv[])
{
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
                "   traits \\\n"
                "   sex_mutation_rate \\\n"
                "   sex_change \\\n"
                "   ngen\n");
        exit(0);
    }

    if(argc != 14)
        INVALID(1,"Wrong number of arguments\n");

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

    sex_mutation_rate = strtod(argv[11],NULL);
    sex_change = !!strtol(argv[12],NULL,0);

    initialize_mutation_parameters(argv);

    uint32 ngen = (uint32)strtoul(argv[13],NULL,0);

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

    env = maximum_weight/2.0;

    initialize_population(sex);

    uint32 i,j;
    printf("Maximum weight: %0.3f\nInitial env: %0.3f\n",
            maximum_weight, env);
    for(i = 0 ; i < ngen && (!i || printf("\n")); i++)
    {
        make_children(scratch,choices,choices_ints);

        /* Print weights */
        printf("env: %8.2f\n",env);
        printf("  no_sex:");
        for(j = 0; j <= (uint32)maximum_weight; j++)
            printf(" %u:%u",j,no_sex_weights[j]);
        printf("\n     sex:");
        for(j = 0; j <= (uint32)maximum_weight; j++)
            printf(" %u:%u",j,sex_weights[j]);

        pick_new_env();
    }
    printf("\n");
    return 0;
}

