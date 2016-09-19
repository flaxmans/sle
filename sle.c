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
#include <gsl_rng.h>            // gnu scientific library //
#include <gsl_randist.h>        // gnu scientific library //
const gsl_rng_type *rngType;    /* generator type */
gsl_rng *rngState; /* rng instance */

int loci_count_0 = 0;
int loci_count_1 = 0;
int loci_count_2 = 0;
int loci_count_many = 0;

// constants or defaults
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

// magic numbers
#define nGENOTYPES 3
#define DIPLOID 2


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
int getAndSetRNGseed(void);
void migration(int *genotypeCounts, int *nInEachPatch);
void parseCommandLine(int argc, char *argv[], char *progname);
void printParametersToFile(long int t, int nSamplesGot, int seed);
void recordData(long int t);
void reproduction(int *genotypeCounts, double *fitnesses, int *nInEachPatch);
int myRound(double x);
int seed_gen(void);
void setUpPopulationAndDataFiles(int *genotypeCounts, int *nInEachPatch, double *fitnesses);
double singleLocusEquilibrium(double s, double m);
void usage(char *progname);


// beginning of main
int 
main(int argc, char *argv[])
{
    // read in optional command line arguments ...
    char *progname = argv[0];
    parseCommandLine(argc, argv, progname);
    
    // set up RNG:
    int seed;
    seed = getAndSetRNGseed();

    // core data structures:
    int genotypeCounts[(nPATCHES * nGENOTYPES)]; // homo, het, other homo for each deme
    int nInEachPatch[nPATCHES]; // count of total individuals in each patch
    double fitnesses[(nPATCHES * nGENOTYPES)];
    
    // initialize population:
    setUpPopulationAndDataFiles( genotypeCounts, nInEachPatch, fitnesses );

    // main loop of work
    int nSamplesGot = 0;
    long int t = 0;
    do {
        
        t++;            // increment generation
        migration( genotypeCounts, nInEachPatch );
        reproduction( genotypeCounts, fitnesses, nInEachPatch ); // according to fitness; soft selection
        
        if ( ((t % TS_RECORDING_FREQ) == 0) && (t > BURN_IN_PERIOD) ) {
            recordData( t );
            nSamplesGot++;
        }
        
    } while ( nSamplesGot < nSAMPLES_TO_GET );
    
    // record metadata:
    printParametersToFile(t, nSamplesGot, seed);
    
    // /*  check test
    int *ipt, i;
    ipt = genotypeCounts;
    for ( i = 0; i < nPATCHES; i++ ) {
        fprintf(stdout, "%i\t%i\t%i\t%i\t%i\n", *ipt, *(ipt+1), *(ipt+2), (*ipt + *(ipt+1) + *(ipt+2)), nInEachPatch[i]);
        ipt += nGENOTYPES;
    }
    // end check test */
    
    return 0;
}


int
getAndSetRNGseed(void)
{
    int seed;
    
    if (DETERMINISTIC) {
        int rcount;
        FILE *fpt;
        fpt = fopen("RnumSeed.txt","r");
        if (fpt == NULL) {
            perror("Can't open RnumSeed.txt");
            exit(-1);
        }
        
        rcount = fscanf(fpt,"%i",&seed);
        if ( rcount ) {
            fclose(fpt);
        }
        else {
            fprintf(stderr, "\n\n\tError! nothing read from file! Exiting!\n\n");
            exit(-1);
        }
        //fprintf(stderr, "\n\nSeed = %i\n\n",seed);
    }
    else {
        seed = seed_gen();
    }
    
    gsl_rng_env_setup(); // set up the environment variables for the RNG
    rngType = gsl_rng_mt19937; // Mersenne Twister// for default you can say T = gsl_rng_default;
    rngState = gsl_rng_alloc(rngType);
    gsl_rng_set(rngState, seed); // initialize
    
    return seed;
}


void
migration(int *genotypeCounts, int *nInEachPatch)
{
    int i, j, leaving[nPATCHES][nGENOTYPES], joining[nPATCHES][nGENOTYPES];
    int *ipt, destination, newTot;
    
    
    if ( !TWO_DEME ) {
        fprintf(stderr, "\nError!  Migration code for scenarios other than two-deme not written yet! Sorry!\n");
        exit(-1);
    }
    
    ipt = genotypeCounts;
    for ( i = 0; i < nPATCHES; i++ ) {
        if (i == 1)
            destination = 0;
        else
            destination = 1;
        for ( j = 0; j < nGENOTYPES; j++ ) {
            leaving[i][j] = gsl_ran_binomial( rngState, SD_MOVE, *(ipt + j) );
            joining[destination][j] = leaving[i][j];
        }
        ipt += nGENOTYPES;
    }
    
    ipt = genotypeCounts;
    newTot = 0;
    for ( i = 0; i < nPATCHES; i++ ) {
        for ( j = 0; j < nGENOTYPES; j++ ) {
            *(ipt + j) = *(ipt + j) + joining[i][j] - leaving[i][j];
            if ( *(ipt + j) < 0 ) {
                fprintf(stderr, "\nError in migration()! Negative individuals!\n");
                exit(-1);
            }
            newTot += *(ipt + j);
        }
        
        // update number in each patch
        nInEachPatch[i] = 0;
        for ( j = 0; j < nGENOTYPES; j++ )
            nInEachPatch[i] = nInEachPatch[i] + *(ipt + j);
        
        // increment pointer
        ipt += nGENOTYPES;
    }
    if ( newTot != N ) {
        fprintf(stderr, "\nError in migration():\n\tNumbers do NOT add up!\n\tnewTot = %i, N = %i\n", newTot, N);
        exit(-1);
    }
    
    
}


void
parseCommandLine(int argc, char *argv[], char *progname)
{
    _Bool burnInSetByUser = 0;
    int ch;
    while ((ch = getopt(argc, argv, "B:D:d:F:H:K:M:m:N:n:O:P:s:T:t:?")) != -1) {
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
            case 't':
                TS_RECORDING_FREQ = atoi(optarg);
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
    else
        N = INITIAL_POPULATION_SIZE;

    if ( nSAMPLES_TO_GET < 1 ) {
        fprintf(stderr, "\nError in parameter choices (-n):\n\tnSAMPLES_TO_GET (= %i) should be > 0\n", nSAMPLES_TO_GET);
        exit(-1);
    }
    
    if ( PATCHES < 2 ) {
        fprintf(stderr, "\nError in parameter choices (-P):\n\tnPATCHES (= %i) should be >= 2\n", PATCHES);
        exit(-1);
    }
    nPATCHES = pow(PATCHES,DIM);
    if ( TWO_DEME ) {
        if ( nPATCHES > 2 ) {
            fprintf(stdout, "\nWarning: nPATCHES set to 2, DIM set to 1, because TWO_DEME mode turned on and overrode other options...\n");
        }
        PATCHES = 2;
        nPATCHES = 2;
        DIM = 1;
    }
    
    if ( S_COEFF < 0.0 ) {
        fprintf(stderr, "\nError in parameter choices (-s):\n\tS_COEFF (= %f) should be >= 0.0\n", S_COEFF);
        exit(-1);
    }
    
    if ( TS_RECORDING_FREQ < 1 ) {
        fprintf(stderr, "\nError in parameter choices (-t):\n\tTS_RECORDING_FREQ (= %i) should be >= 0\n", TS_RECORDING_FREQ);
        exit(-1);
    }
    
}

void
printParametersToFile(long int t, int nSamplesGot, int seed)
{
    FILE *pp;
    
    pp = fopen("parameters.m", "w"); // plain text; ready to read by MATLAB
    
    fprintf(pp, "codeVersion = %s\n", version);
    fprintf(pp, "RnumSeed = %i\n", seed);
    
    fprintf(pp, "BURN_IN_PERIOD = %i;\n", BURN_IN_PERIOD);
    fprintf(pp, "DIM = %i;\n", DIM);
    fprintf(pp, "DETERMINISTIC = %i;\n", DETERMINISTIC);
    fprintf(pp, "FIXED_N = %i;\n", FIXED_N);
    fprintf(pp, "H = %E;\n", H);
    fprintf(pp, "K = %E;\n", K);
    fprintf(pp, "MOSAIC = %i;\n", MOSAIC);
    fprintf(pp, "SD_MOVE = %E;\n", SD_MOVE);
    fprintf(pp, "INITIAL_POPULATION_SIZE = %i;\n", INITIAL_POPULATION_SIZE);
    fprintf(pp, "nSAMPLES_TO_GET = %i;\n", nSAMPLES_TO_GET);
    fprintf(pp, "OFFSPRING_IN_RANDOM_LOCATIONS = %i;\n", OFFSPRING_IN_RANDOM_LOCATIONS);
    fprintf(pp, "PATCHES = %i;\n", PATCHES);
    fprintf(pp, "nPATCHES = %i;\n", nPATCHES);
    fprintf(pp, "S_COEFF = %E;\n", S_COEFF);
    fprintf(pp, "TWO_DEME = %i;\n", TWO_DEME);
    fprintf(pp, "TS_RECORDING_FREQ = %i;\n", TS_RECORDING_FREQ);
    
    fprintf(pp, "\nN = %i;\n", N);
    fprintf(pp, "t = %li;\n", t);
    fprintf(pp, "nSamplesGot = %i;\n", nSamplesGot);
    
    fclose(pp);
}


void
recordData(long int t)
{
    
}


void
reproduction(int *genotypeCounts, double *fitnesses, int *nInEachPatch)
{
    double fitnessWeights[nPATCHES][nGENOTYPES], *dpt, fitTotal, a, b, c;
    int i, j, *ipt;
    int newCounts[(nPATCHES * nGENOTYPES)];
    
    ipt = genotypeCounts;
    dpt = fitnesses;
    
    for ( i = 0; i < nPATCHES; i++ ) {
        fitTotal = 0.0;
        // calculate raw weights:
        for ( j = 0; j < nGENOTYPES; j++ ) {
            fitnessWeights[i][j] = ((double) *(ipt + j)) * (*(dpt + j));
            fitTotal += fitnessWeights[i][j];
        }
        // normalize to sum to one to make them probabilities:
        for ( j = 0; j < nGENOTYPES; j++ )
            fitnessWeights[i][j] = fitnessWeights[i][j] / fitTotal;
        
        if ( nGENOTYPES > 3 ) {
            fprintf(stderr, "\nError in reproduction: scheme assumes three genotypes total!  Sorry!\n");
            exit(-1);
        }
        
        a = fitnessWeights[i][0];
        b = fitnessWeights[i][1];
        c = fitnessWeights[i][2];
        
        *ipt = myRound( ((a * a) + (a * b) + (0.25 * b * b)) * ((double) nInEachPatch[i]) );
        *(ipt + 2) = myRound( ((0.25 * b * b) + (b * c) + (c * c)) * ((double) nInEachPatch[i]) );
        *(ipt + 1) = nInEachPatch[i] - ( *ipt + (*(ipt + 2)) );
        
        // increment pointers:
        ipt += nGENOTYPES;
        dpt += nGENOTYPES;
    }
    
    
}


int
myRound(double x)
{
    return ((int) (x + 0.5));
}


int
seed_gen(void)
{
    /* use calendar time to seed random number generator.  Code adopted from Schildt's textbook */
    
    int stime;
    long ltime;
    FILE *rseed;
    
    /* get the calendar time */
    ltime=time(NULL);
    stime=(unsigned) ltime/2;
    
    
    // generate and store random number seed
    rseed = fopen("RnumSeed.txt","w");
    fprintf(rseed,"%i\n",stime);
    fclose(rseed);
    
    return stime;
}


void
setUpPopulationAndDataFiles(int *genotypeCounts, int *nInEachPatch, double *fitnesses )
{
    int i, j, *ipt, hereNow, totSoFar, midIndex;
    double q1, q0, p1, p0, puse, nPerPatch;
    
    q1 = singleLocusEquilibrium(S_COEFF, SD_MOVE); // accurate only for two deme; frequency of derived allele where it is favored
    if ( p1 <= 0.0 || p1 >= 1.0 ) {
        fprintf(stderr, "\nError in setUpPopulationAndDataFiles():\n\tp1 (= %f) out of bounds!\n", p1);
        exit(-1);
    }
    p1 = 1.0 - q1; // ancestral allele where it is NOT favored
    
    p0 = q1; // ancestral allele where it is favored
    q0 = 1.0 - p0; // derived allele NOT favored
    
    if ( (nPATCHES % 2) == 1 )
        midIndex = (nPATCHES - 1) / 2;
    else
        midIndex = nPATCHES/2;
    
    
    nPerPatch = (double) myRound( ((double) N) / ((double) nPATCHES) );
    fprintf(stdout, "\n%i\t%i\t%f\n", N, nPATCHES, nPerPatch);
    
    ipt = genotypeCounts;
    totSoFar = 0;
    for ( i = 0; i < (nPATCHES-1); i++ ) {
        nInEachPatch[i] = (int) nPerPatch; // number here
        totSoFar += nInEachPatch[i];
        if ( ((nPATCHES % 2) == 1) && (i == midIndex) )
            puse = 0.5;
        else if ( i < (nPATCHES / 2) )
            puse = p0;  // half of habitat good for ancestral
        else
            puse = p1;  // half of habitat good for derived
        
        // Use Hardy-Weinberg equilibrium:
        *ipt = myRound((puse * puse) * nPerPatch); // homozygous ancestral
        *(ipt + 1) = myRound((2.0 * puse * (1.0 - puse)) * nPerPatch); // heterozygous
        hereNow = (*ipt) + (*(ipt + 1));
        *(ipt + 2) = ((int) nPerPatch) - hereNow; // homozygous derived
        ipt += nGENOTYPES; // increment genotype pointer
    }
    nInEachPatch[(nPATCHES - 1)] = N - totSoFar;
    *ipt = myRound((p1 * p1) * ((double) nInEachPatch[(nPATCHES - 1)]));
    *(ipt + 1) = myRound((2.0 * p1 * (1.0 - p1)) * ((double) nInEachPatch[(nPATCHES - 1)]));
    *(ipt + 2) = nInEachPatch[(nPATCHES - 1)] - (*ipt + (*(ipt + 1)));
    
    // /*  check test
    fprintf(stdout, "\n%f\t%f\t%f\t%f\n", p0, q0, p1, q1);
    ipt = genotypeCounts;
    for ( i = 0; i < nPATCHES; i++ ) {
        fprintf(stdout, "%i\t%i\t%i\t%i\t%i\n", *ipt, *(ipt+1), *(ipt+2), (*ipt + *(ipt+1) + *(ipt+2)), nInEachPatch[i]);
        ipt += nGENOTYPES;
    }
    // end check test */
    
    if ( !TWO_DEME ) {
        fprintf(stderr, "\nError in setting up fitnesses: not configured for anything except two deme!  Sorry!\n");
        exit(-1);
    }
    fitnesses[0] = 1.0 + S_COEFF;
    fitnesses[1] = 1.0 + ( (1.0 - H) * S_COEFF );
    fitnesses[2] = 1.0;
    fitnesses[3] = 1.0;
    fitnesses[4] = 1.0 + ( H * S_COEFF );
    fitnesses[5] = 1.0 + S_COEFF;
    
}


double
singleLocusEquilibrium(double s, double m)
{
    // see Mathematica notebook for solution
    return (  (m * (0.5 + 0.75 * s) - 0.125 * s - 0.125 * sqrt(pow(s,2) - 4 * m * pow(s,2) + 4 * pow(m,2) * pow((2 + s),2)))/((-0.25 + m) * s)  );
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
    
    fprintf(stdout, "\n\t-t <int>\tTime sampling recording interval.\n\t\t\tDefault is %i\n", TS_RECORDING_FREQ_DEFAULT);
    
}


