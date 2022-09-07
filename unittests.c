#include "recom_fluct_envir_many.h"

extern void parse_rates(char *s);
extern void parse_contrib(char *s);
extern void parse_traits(char *s);

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

    bitstr *indvs = malloc(sizeof(bitstr)*(nindiv));

    uint32 i;
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

    /* Does parse_traits work correctly? */
    parse_traits(" 13 26 39 ");
    assert( traits[0] == 13
            && traits[1] == 26
            && traits[2] == 39 );
    assert( ntraits == 4 );

    assert( weight(gtypes[1]) == sqrt(1.0*1.0 + 20.0*20.0 + 24.0*24.0 + 23.0*23.0) );

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

    /* Test "null" strings for traits */

    uint32 round;
    char *test_trait_string = (char*)malloc(0x100);;
    for(round = 1; round <= 6; round++)
    {
        uint32 tntraits; // True number of traits
        switch(round)
        {
            case 1:
                {
                    char ex[] = "";
                    strncpy(test_trait_string,ex,sizeof(ex)+1);
                    tntraits = 1;
                }
                break;
            case 2:
                {
                    char ex[] = "60";
                    strncpy(test_trait_string,ex,sizeof(ex)+1);
                    tntraits = 1;
                }
                break;
            case 3:
                {
                    char ex[] = "0";
                    strncpy(test_trait_string,ex,sizeof(ex)+1);
                    tntraits = 2;
                }
                break;
            case 4:
                {
                    char ex[] = "  0  60";
                    strncpy(test_trait_string,ex,sizeof(ex)+1);
                    tntraits = 2;
                }
                break;
            case 5:
                {
                    char ex[] = "  ";
                    strncpy(test_trait_string,ex,sizeof(ex)+1);
                    tntraits = 1;
                }
        }

        free(traits);
        parse_traits(test_trait_string);
        assert( ntraits == tntraits );

        for(i = 0 ; i < nindiv ; i++)
            insert(indvs[perm[i]]);

        linearize_and_tally_weights();
        assert( array.w[0] == A(67,0) );
        assert( array.w[1] == A(69,0) );
        assert( array.w[2] == A(71,0) );
        assert( array.w[3] == A(69,0) );
        assert( array.w[4] == A(71,0) );
        assert( array.w[5] == A(69,0) );
        assert( array.w[6] == A(71,0) );
        assert( array.w[7] == A(73,0) );
        assert( array.w[8] == A(71,0) );
        assert( array.w[9] == A(68,1) );
        assert( array.w[10] == A(70,1) );
        assert( array.w[11] == A(68,1) );
        assert( array.w[12] == A(70,1) );
        assert( array.w[13] == A(72,1) );
        assert( array.w[14] == A(70,1) );
        assert( array.w[15] == A(72,1) );
        assert( array.w[16] == A(70,1) );
        assert( array.w[17] == 2*A(72,1) );

    }
    free(gtypes);
    free(indvs);
    free(test_trait_string);

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

    /* Casual prng */
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

    uint32 j,k;
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
    free(traits);

    env = maximum_weight/2.0;
    shift_size = 7.5; // Large SD
    /* Ensure that `env' stays in range
     * ( assert statement in function pick_new_env )
     * */
    for(i = 1000000; i; i--)
        pick_new_env();

    gsl_rng_free(rng);
    
    /* TODO
     * Recombination, among other things
     * */
}
