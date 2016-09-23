# sle

### *S*ingle *L*ocus *E*xpectations for selected loci with the fitness and migration scheme of the bu2s model

This is a fast tool for generating samples of the stochastic trajectory of a single locus under migration, selection, and drift when selection is due to ecologically divergent selection acting on two populations.  The useful parts here are (1) estimating the variance in a number of population genetic metrics (the means are usually obtainable from analytic theory) and (2) the mean and variance of fixation times in cases when migration + drift swamps selection and leads to loss.

____________________
###System requirements
* gcc (C compiler)
* GNU Scientific Library
* a Linux/UNIX-style terminal (or terminal emulator)

### Compiling the program
The program can be compiled from source with the following command (in a Linux/UNIX terminal; don't type the dollar sign; that represents the command prompt):

	$ gcc -lm -I/usr/local/include/gsl -L/usr/local/lib -lgsl -lgslcblas -O3 sle.c -o sle

For convenience, this command is also given in the file 'compileCommand.txt' in this repository.


###Explantion of data and metadata files

* parameters.m: lists all values of all parameters and some end state variables.  This is can be called as a Matlab script:
<br>
` > run('./parameters.m')`.  	
Even though the file has equals signs rather than the R assignment operator <br> ('<-'), it can also be read by R with the command 
<br>
` > source('parameters.m')`.

* RnumSeed.txt: the random number seed used to initialize the RNG.

* AlleleFreqTS.txt: this the data output.  This is a space-delimited, plain-text file, with the following header row indicating the columns of data present: time, fixation, s, m, RnumSeed, nSamplesGot, qGlobal, q0, q1, clineWidth, FST, which correspond to the following quantities:

	1. time = generation sample was taken

	2. fixation = 1 if this time sample corresponded to when one allele fixed; = 0 otherwise.
	
	3.  s = selection coefficient
	
	4.  m = migration rate between two demes
	
	5.  RnumSeed = value of random number seed used to seed the RNG (gsl)
	
	6.  nSamplesGot = integer (n) indicating that this is the nth sample from this particular simulation run
	
	7.  qGlobal = global frequency (i.e, both demes combined) of the allele that is best in deme 1. (demes are referred to as "deme 0" and "deme 1")
	
	8. q0 = frequency of allele in deme 0
	
	9. q1 = frequency of allele in deme 1
	
	10. clineWidth = 1 / (q1 - q0)
	
	11. FST = fixation index, computed as (HT - HS)/HT






