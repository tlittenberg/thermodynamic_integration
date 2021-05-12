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

Required:

    -f FILE #filename for input chain
    -n INT #number of PT chains

Optional:

    -h #print this message and exit
    -t INT # thin chains by factor INT
    
## Authors
 + Tyson B. Littenberg 
 + Neil J. Cornish
 
 ## License

This project is licensed under the GNU General Public License v3.0  - see the LICENSE.md file for details
