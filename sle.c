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
#define K_DEFAULT 5000.0 // default carrying capacity per patch; only used if FIXED_N = 0
#define H_DEFAULT 0.5 // dominance coefficient
#define MOSAIC_DEFAULT 0 // landscape as mosaic of patches; if 0, then it's a gradient
#define SD_MOVE_DEFAULT 0.05 // if TWO_DEME is used, this IS the gross migration rate
#define OIRL_DEFAULT 1 // offspring in random locations or not
#define TS_RECORDING_FREQ_DEFAULT 500 // generation interval for recording
#define BURN_IN_PERIOD_DEFAULT 20000 // how long to run before recording
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
int BURN_IN_PERIOD = BURN_IN_PERIOD_DEFAULT;
int nSAMPLES_TO_GET = nSAMPLES_TO_GET_DEFAULT;

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
    int ch;
    while ((ch = getopt(argc, argv, "B:D:d:F:H:K:M:m:N:n:O:P:s:T:?")) != -1) {
        switch (ch) {
            case 'B':
                BURN_IN_PERIOD = atoi(optarg);
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
    
}


