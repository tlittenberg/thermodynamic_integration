# thermodynamic_integration
Spline-based integration of likelihood samples from PTMCMC analysis

## Description
The code implements the algorithm described in section 6.1.2 of the paper https://iopscience.iop.org/article/10.1088/0264-9381/32/13/135012. The basic idea is to fit a cubic spline to the evidence integrand of thermodynamic integration. The code uses a bootstrap procedure to estimate the variance and covariance of the likelihood chains from a parallel tempered MCMC. This covariance matricx defines a likelihood function that is used by a trans-dimensionanal Reversible Jump Markov Chain Monte Carlo (RJMCMC) algorithm that fits cubic splines with variable numbers of control points (knots) to the thermodynamic entegrand. Each step of the RJMCMC yields a value for the evidence integral. The end result is a posterior distribution for the evidence, from which point estimates such as the mean and variance of the evidence can be computed.

## Getting Started

### Dependencies
 + cmake
 + gsl

### Installing
`./install.sh /path/to/install/destination`

### Executing the program

Example Usage: `./thermodynamic_integration -f /path/to/file -n Nchains`

Required arguments:

    -f FILE #filename for input chain
    -n INT #number of PT chains

Optional arguments:

    -h #print this message and exit
    -t INT # thin chains by factor INT
    
### Input file format
The program reads ASCII space-delimted tabulated data containing the log-likelihood samples from a parallel tempered MCMC analysis. 
Rows correspond to chain steps, columns correspond to the different tempered chains.
The first row of the file has the temperature for each of the parallel chains with each column increasing from lowest to highest temperature.

     Temp_0    Temp_1    Temp_2    ...
     logL_0[0] logL_1[0] logL_2[0] ... 
     logL_0[1] logL_1[1] logL_2[1] ...
     ...
     ...
     
### Output file format
 + `integrand.dat` : 2 column ASCII, 1st column is `log10 T` and the 2nd column is `<log L>`
 + `evidence_chain.dat` : 2 column ASCII, 1st column is MCMC integration step, 2nd column is inferred evidence.
 + `integrand_quantiles.dat` : 6 column ASCII, 1st column is `log10 T`, 2nd is the median integrand, 3rd and 4th are 25th and 75th quantile, and 5th and 6th are the 5th and 95th quantile.
 + `evidence.dat` : 2 column ASCII, 1st column is the average evidence, 2nd column is the 1-sigma uncertainty on the evidence.
    
## Authors
 + Tyson B. Littenberg 
 + Neil J. Cornish
 
 ## License

This project is licensed under the GNU General Public License v3.0  - see the LICENSE.md file for details
