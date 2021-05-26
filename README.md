# thermodynamic_integration
Spline-based integration of likelihood samples from PTMCMC analysis

## Description

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
    
## Authors
 + Tyson B. Littenberg 
 + Neil J. Cornish
 
 ## License

This project is licensed under the GNU General Public License v3.0  - see the LICENSE.md file for details
