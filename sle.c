/*
 *  bu2s.c
 *  
 *
 *  Created by Samuel Melvin Flaxman on 4/10/12.
 *  Copyright 2012-2016 Samuel Melvin Flaxman. All rights reserved.
 *
 */

const char *version = "sle_1.0";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "../MT/dSFMT.h"

int loci_count_0 = 0;
int loci_count_1 = 0;
int loci_count_2 = 0;
int loci_count_many = 0;

// code for using Mersenne Twister RNG
dsfmt_t dsfmt;
#define seedRand(s) dsfmt_init_gen_rand(&dsfmt, s)
#define	randU() dsfmt_genrand_close_open(&dsfmt)
#define	randI() (unsigned)dsfmt_genrand_uint32(&dsfmt)


// constants that have to be set here in source code
#define DIM_DEFAULT 1 // dimensionality of habitat; set to 1 or 2
#define PATCHES_DEFAULT 2 // number of patches across the habitat in a single dimension; see also nPATCHES below, which is the actual number = PATCHES^D
#define TWO_DEME_DEFAULT 1 // if this is 1, then D above should be 1 and PATCHES should be 2, and space will be treated as discrete rather than as continuous
#define DETERMINISTIC_DEFAULT 1
#define INITIAL_POPULATION_SIZE_DEFAULT 5000
#define S_COEFF_DEFAULT 0.02 // mean value of selection coefficients drawn from -exponential distribution
#define FIXED_N_DEFAULT 1 // whether N will be fixed or variable (latter would be density dependent regulation)
#define K_DEFAULT 2500.0 // default carrying capacity per patch; only used if FIXED_N = 0
#define H_DEFAULT 0.5 // dominance coefficient
#define MOSAIC_DEFAULT 0 // landscape as mosaic of patches; if 0, then it's a gradient
#define SD_MOVE_DEFAULT 0.05 // if TWO_DEME is used, this IS the gross migration rate
#define OIRL_DEFAULT 1 // offspring in random locations or not
#define TS_RECORDING_FREQ_DEFAULT 500 // generation interval for recording
#define nSAMPLES_TO_GET_DEFAULT 1000 // how many samples to get

// globals
int DIM = DIM_DEFAULT;
int PATCHES = PATCHES_DEFAULT;
int TWO_DEME = TWO_DEME_DEFAULT;
_Bool DETERMINISTIC = DETERMINISTIC_DEFAULT;
int INITIAL_POPULATION_SIZE = INITIAL_POPULATION_SIZE_DEFAULT;
double S_COEFF = S_COEFF_DEFAULT;
_Bool FIXED_N = FIXED_N_DEFAULT;
double K = K_DEFAULT;
double H = H_DEFAULT;
_Bool MOSAIC = MOSAIC_DEFAULT;
double SD_MOVE = SD_MOVE_DEFAULT;
_Bool OFFSPRING_IN_RANDOM_LOCATIONS = OIRL_DEFAULT;
int TS_RECORDING_FREQ = TS_RECORDING_FREQ_DEFAULT;
int BURN_IN_PERIOD;
int nSAMPLES_TO_GET = nSAMPLES_TO_GET_DEFAULT;
int N;
int nPATCHES;

// functions
void migration(void);
void parseCommandLine(int argc, char *argv[], char *progname);
void recordData(long int t);
void reproduction(void);
void setUpPopulationAndDataFiles(void);
void usage(char *progname);


// beginning of main
int 
main(int argc, char *argv[])
{
    // read in optional command line arguments ...
    char *progname = argv[0];
    parseCommandLine(argc, argv, progname);
    
    setUpPopulationAndDataFiles();
    
    int nSamplesGot = 0;
    long int t = 0;
    
    do {
        
        t++;            // increment generation
        migration();
        reproduction(); // according to fitness; soft selection
        
        if ( ((t % TS_RECORDING_FREQ) == 0) && (t > BURN_IN_PERIOD) ) {
            recordData( t );
            nSamplesGot++;
        }
        
    } while ( nSamplesGot < nSAMPLES_TO_GET );
    
    
    return 0;
}


void
migration(void)
{
    
}


void
parseCommandLine(int argc, char *argv[], char *progname)
{
    _Bool burnInSetByUser = 0;
    int ch;
    while ((ch = getopt(argc, argv, "B:D:d:F:H:K:M:m:N:n:O:P:s:T:?")) != -1) {
        switch (ch) {
            case 'B':
                BURN_IN_PERIOD = atoi(optarg);
                burnInSetByUser = 1;
                break;
            case 'd':
                DIM = atoi(optarg);
                break;
            case 'D':
                DETERMINISTIC = atoi(optarg);
                break;
            case 'F':
                FIXED_N = atoi(optarg);
                break;
            case 'H':
                H = strtod(optarg, (char **)NULL);
                break;
            case 'K':
                K = strtod(optarg, (char **)NULL);
                break;
            case 'M':
                MOSAIC = atoi(optarg);
                break;
            case 'm':
                SD_MOVE = strtod(optarg, (char **)NULL);
                break;
            case 'N':
                INITIAL_POPULATION_SIZE = atoi(optarg);
                break;
            case 'n':
                nSAMPLES_TO_GET = atoi(optarg);
                break;
            case 'O':
                OFFSPRING_IN_RANDOM_LOCATIONS = atoi(optarg);
                break;
            case 'P':
                PATCHES = atoi(optarg);
                break;
            case 's':
                S_COEFF = strtod(optarg, (char **)NULL);
                break;
            case 'T':
                TWO_DEME = atoi(optarg);
                break;
            case '?':
            default:
                usage(progname);
                exit(-1);
        }
    }
    
    if ( !burnInSetByUser )
        BURN_IN_PERIOD = 4 * INITIAL_POPULATION_SIZE;
    
    if ( DIM != 1 && DIM != 2 ) {
        fprintf(stderr, "\nError in parameter choices (-d):\n\tDIM (= %i) must be 1 or 2\n", DIM);
        exit(-1);
    }
    
    if ( H < 0.0 || H > 1.0 ) {
        fprintf(stderr, "\nError in parameter choices (-H):\n\tH (= %f) must be between 0 and 1\n", H);
        exit(-1);
    }
    
    if ( K < 0.0 ) {
        fprintf(stderr, "\nError in parameter choices (-K):\n\tK (= %f) must be > 0.0\n", K);
        exit(-1);
    }
    
    if ( SD_MOVE < 0.0 ) {
        fprintf(stderr, "\nError in parameter choices (-m):\n\tSD_MOVE (= %f) must be > 0.0\n", SD_MOVE);
        exit(-1);
    }
    
    if ( INITIAL_POPULATION_SIZE < (2*pow(PATCHES, DIM)) ) {
        fprintf(stderr, "\nError in parameter choices (-N):\n\tINITIAL_POPULATION_SIZE (= %i) too small!\n", INITIAL_POPULATION_SIZE);
        exit(-1);
    }

    if ( nSAMPLES_TO_GET < 1 ) {
        fprintf(stderr, "\nError in parameter choices (-n):\n\tnSAMPLES_TO_GET (= %i) should be > 0\n", nSAMPLES_TO_GET);
        exit(-1);
    }
    
    if ( PATCHES < 2 ) {
        fprintf(stderr, "\nError in parameter choices (-P):\n\tnPATCHES (= %i) should be >= 2\n", PATCHES);
        exit(-1);
    }
    nPATCHES = pow(PATCHES,DIM);
    
    if ( S_COEFF < 0.0 ) {
        fprintf(stderr, "\nError in parameter choices (-s):\n\tnS_COEFF (= %f) should be >= 0.0\n", S_COEFF);
        exit(-1);
    }
    
}


void
recordData(long int t)
{
    
}


void
reproduction(void)
{
    
}


void
setUpPopulationAndDataFiles(void)
{
    
}


void
usage(char *progname)
{
    fprintf(stdout, "\nOptions and meanings for %s:\n", progname);
    
    fprintf(stdout, "\n\t-B <int>\tlength of burn-in period (generations).\n\t\t\tDefault is 4 * initial population size\n");
    
    fprintf(stdout, "\n\t-D <0 or 1>\tDETERMINISTIC or not (i.e., use random number seed from\n\t\t\tfile or not).  '1' = true.  Default is %i\n", DETERMINISTIC_DEFAULT);
    
    fprintf(stdout, "\n\t-d <1 or 2>\tDimensionality of habitat.\n\t\t\tDefault is %i\n", DIM_DEFAULT);
    
    fprintf(stdout, "\n\t-F <0 or 1>\tFixed population size or not \n\t\t\t'1' = true.  Default is %i.  See also -K below.\n", FIXED_N_DEFAULT);
    
    fprintf(stdout, "\n\t-H <[0,1]>\tDominance coefficient, H. Must be between 0 and 1.\n\t\t\tDefault is %f\n", H_DEFAULT);
    
    fprintf(stdout, "\n\t-K <num>\tCarrying capacity of each deme.  Only used with -F 0.\n\t\t\tMust be > 0.  Default is %f\n", K_DEFAULT);
    
    fprintf(stdout, "\n\t-M <0 or 1>\tHabitat as mosaic of habitat types or not.\n\t\t\t'1' = mosaic.  '0' = gradient.  Default is %i.\n", MOSAIC_DEFAULT);
    
    fprintf(stdout, "\n\t-m <num>\tStandard deviation of migration distances. In two-deme\n\t\t\tmode (see -T below), this is the migration probability.\n\t\t\tOtherwise, the total length of the habitat is defined to\n\t\t\tbe 1.  Default is %f\n", SD_MOVE_DEFAULT);
    
    fprintf(stdout, "\n\t-N <int>\tInitial global population size (all demes combined).\n\t\t\tDefault is %i\n", INITIAL_POPULATION_SIZE_DEFAULT);
    
    fprintf(stdout, "\n\t-n <int>\tNumber of data sampling points at which to record.\n\t\t\tDefault is %i\n", nSAMPLES_TO_GET_DEFAULT);
    
    fprintf(stdout, "\n\t-O <1 or 0>\tOffspring in random locations or not.\n\t\t\t'1' = random within deme; '0' = same location as one\n\t\t\tparent.  Default is %i\n", OIRL_DEFAULT);
    
    fprintf(stdout, "\n\t-P <int>\tNumber of patches along each dimension of habitat.\n\t\t\tDefault is %i\n", PATCHES_DEFAULT);

    fprintf(stdout, "\n\t-s <num>\tSelection coefficient.  Should be > 0.\n\t\t\tDefault is %f\n", S_COEFF_DEFAULT);
    
    fprintf(stdout, "\n\t-T <1 or 0>\tTwo-deme mode or not.  '0' = continuous space.\n\t\t\t'1' = two-deme, discrete-space model.\n\t\t\tDefault is %i\n", TWO_DEME_DEFAULT);
}


