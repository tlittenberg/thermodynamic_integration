//
//  main.c
//  thermodynamic_integration
//
//  Created by Littenberg, Tyson B. (MSFC-ZP12) on 5/7/21.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>

#define BUFFER_SIZE 1024 //!< max size of `char` buffers

//chain file structure
struct Chains
{
    size_t NC; //!< number of chains
    size_t NMCMC; //!<number of MCMC samples in chain
    size_t thin; //!<factor for thinning chain
    FILE *chainFile; //!<input chain file
};

void printUsage(const char *program);
void parse_command_line(int argc, char **argv, struct Chains *chain);
void ThermodynamicIntegration(double *temperature, double **loglikelihood, size_t NMCMC, size_t NC, double *logEvidence, double *varLogEvidence);


int main(int argc, char *argv[])
{
    
    //chain structure
    struct Chains *chain = malloc(sizeof(struct Chains));
    
    //array of temperatures
    double *temperature;
    
    //2D array of log likelihood samples
    double **loglikelihood;
    
    //estimated log marginalized likelihood
    double logEvidence;
    
    //variance on estimate of logEvidence
    double varLogEvidence;
    
    
    /*
     Parse command line
     */
    parse_command_line(argc, argv, chain);
    
    /*
     Count lines of chain file
     */
    char* line;
    char lineBuffer[BUFFER_SIZE];
    size_t NMCMC = 0;
    while((line = fgets(lineBuffer, BUFFER_SIZE, chain->chainFile)) != NULL) NMCMC++;
    rewind(chain->chainFile);
    
    //first line is temperatures
    NMCMC--;
    
    //thin chain
    chain->NMCMC = (size_t)floor((double)NMCMC/(double)chain->thin);
    
    /*
     Allocate memory to store input data
     */
    temperature = malloc(chain->NC*sizeof(double));
    loglikelihood = malloc(chain->NC*sizeof(double *));
    for(size_t n=0; n<chain->NC; n++) loglikelihood[n] = malloc(chain->NMCMC*sizeof(double));
    
    /*
     Parse file
     */
    for(size_t n=0; n<chain->NC; n++) fscanf(chain->chainFile,"%lg",&temperature[n]);
    for(size_t m=0; m<NMCMC; m+=chain->thin)
    {
        for(size_t n=0; n<chain->NC; n++)
        {
            fscanf(chain->chainFile,"%lg",&loglikelihood[n][m/chain->thin]);
        }
    }
    
    /*
     Do the thing
     */
    ThermodynamicIntegration(temperature, loglikelihood, chain->NMCMC, chain->NC, &logEvidence, &varLogEvidence);
    
    
    /*
     Tell us how you did
     */
    printf("Spline integration %.2f +/- %.2f\n",logEvidence,sqrt(varLogEvidence));
    FILE *evidence = fopen("evidence.dat","w");
    fprintf(evidence,"%.2f %.2f\n",logEvidence,sqrt(varLogEvidence));
    fclose(evidence);
    
    return 0;
}


void printUsage(const char *program)
{
    fprintf( stdout,"\n");
    fprintf( stdout, "Thermodynamic Integration:\n\n");
    fprintf( stdout, "  Spline-based integration of likelihood samples \n");
    fprintf( stdout, "  from PTMCMC analysis\n");
    fprintf( stdout, "\n");
    fprintf( stdout, "Example Usage: %s -f /path/to/file -n 16 \n", program );
    fprintf( stdout, "\n");
    fprintf( stdout, "  Required:\n");
    fprintf( stdout, "     -f, --file=FILE    filename for input chain\n" );
    fprintf( stdout, "     -n, --nchains=INT  number of PT chains\n" );
    fprintf( stdout, "\n");
    fprintf( stdout, "  Optional:\n");
    fprintf( stdout, "    [-h, --help]        print this message and exit\n" );
    fprintf( stdout, "    [-t, --thin=INT]    thin chains by factor INT\n");
    fprintf( stdout,"\n");
}

void parse_command_line(int argc, char **argv, struct Chains *chain)
{
    chain->thin=1;
    
    int opt = 0;
    while ((opt = getopt(argc, argv, "f:n:t:")) != -1)
    {
        switch(opt)
        {
            case 'f':
                chain->chainFile = fopen(optarg,"r");
                if(chain->chainFile==NULL)
                {
                    printUsage(argv[0]);
                    printf("Error: Cannot open requested file %s\n\n",optarg);
                    exit(1);
                }
                break;
                
            case 'h': // help
                printUsage(argv[0]);
                exit(0);
                break;
                
            case 't':
                chain->thin = (size_t)atoi(optarg);
                break;
                
            case 'n':
                chain->NC = (size_t)atoi(optarg);
                break;
                
            case ':':
                if(optopt == 'f')
                {
                    printUsage(argv[0]);
                    printf("Error: Need file name\n\n");
                    exit(1);
                }
                else if(optopt == 'n')
                {
                    printUsage(argv[0]);
                    printf("Error: Need to specify number of chains\n\n");
                    exit(1);
                }
                break;
                
            default:
                printUsage(argv[0]);
                fprintf(stderr,"Error: Unknown error while parsing options\n\n" );
                exit(1);
        }
    }
    //check that you got what you need
    if(chain->chainFile==NULL)
    {
        printUsage(argv[0]);
        printf("Error: Missing required arguments -f /path/to/filename\n\n");
        exit(1);
    }else if(chain->NC==0)
    {
        printUsage(argv[0]);
        printf("Error: Missing required arguments -n Nchain\n\n");
        exit(1);
    }
}


static void bootcorrposdef(double **data, int NC, int N, double **cmat, double *mu, gsl_rng *seed)
{
    double av;
    double x, y, shrink;
    int i, j, k, cmin, icm, cmax;
    int ii, jj, kk, nn;
    int *cyc;
    int cycles;
    int so, sc, cp;
    int b, bs, bsteps;
    int *cnts;
    double **fmean;
    double *var, *mx;
    
    // This code uses the "threshold bootstrap" method to compute the covariance matrix
    // See http://www.sciencedirect.com/science/article/pii/S0377221700002095 and
    /* @inproceedings{Kim:1993:TBN:256563.256697,
     author = {Kim, Y. B. and Willemain, T. R. and Haddock, J. and Runger, G. C.},
     title = {The Threshold Bootstrap: A New Approach to Simulation Output Analysis},
     booktitle = {Proceedings of the 25th Conference on Winter Simulation},
     series = {WSC '93},
     year = {1993},
     isbn = {0-7803-1381-X},
     location = {Los Angeles, California, USA},
     pages = {498--502},
     numpages = {5},
     url = {http://doi.acm.org/10.1145/256563.256697},
     doi = {10.1145/256563.256697},
     acmid = {256697},
     publisher = {ACM},
     address = {New York, NY, USA},
     } */
    
    // The "shrinkage" method described here is uses to ensure the covariance matrices are
    // postive definite and well balanced https://www.stat.wisc.edu/courses/st992-newton/smmb/files/expression/shrinkcov2005.pdf
    
    bsteps = pow(10,floor(log10(N))-2); //want O(100) samples in each bstep
    if(bsteps < 10)
    {
        printf("Error: Not enough likelihood samples to estimate covariance matrix\n");
        printf("       Reduce chain thinning\n");
        exit(1);
    }
    
    var = malloc(NC*sizeof(double));
    cnts = malloc(NC*sizeof(int));
    
    for(ii=0; ii< NC; ii++)
    {
        mu[ii] = 0.0;
        var[ii] = 0.0;
        for(i=0; i< N; i++)
        {
            mu[ii] += data[ii][i];
            var[ii] += data[ii][i]*data[ii][i];
        }
        mu[ii] /= (double)(N);
        var[ii] /= (double)(N);
        var[ii] -= mu[ii]*mu[ii];
        var[ii] = sqrt(var[ii]);
        
    }
    
    b = 8;
    // if b > 1, then what I'm calling cycles are actually chunks, made up of b cycles
    
    icm = 0;
    do
    {
        kk = 0;
        nn = 0;
        // cycle lengths for individual chains
        cmin = 100000;
        cmax = 0;
        for(ii=0; ii< NC; ii++)
        {
            cycles = 0;
            so = -1;
            if((data[ii][0]-mu[ii]) > 0.0) so = 1;
            cp = 0;
            
            for(i=0; i< N; i++)
            {
                sc = -1;
                if((data[ii][i]-mu[ii]) > 0.0) sc = 1;
                if(sc != so) cp++;   // count the change points
                so = sc;
                if(cp == 2*b)
                {
                    cycles++;
                    cp = 0;
                }
            }
            
            cnts[ii] = cycles;
            if(cycles < cmin)
            {
                cmin = cycles;
                icm = ii;
            }
            if(cycles > cmax)
            {
                cmax = cycles;
            }
            
            kk += cycles;
            nn += cycles*cycles;
            
        }
        
        // looking for between 100 and 1000 chunks
        // This assumes that we have a decent amount of data to begin with
        
        if(cmin <  100) b /= 2;
        if(cmin > 1000) b *= 2;
        
        if(b < 2)
        {
            b = 2;
            cmin = 500;
        }
        
        
    } while(cmin < 100 || cmin > 1000);
    
    // record the change point for each cycle
    // if cyc = i, then one cycle ends at i and the next starts at i+1
    cyc = malloc((N/2)*sizeof(int));
    
    cyc[0] = -1;
    
    
    
    av = 0.0;
    
    for(i=0; i< N; i++)
    {
        av += data[icm][i];    // correlation length determined by most corelated
    }
    av /= (double)(N);
    
    
    cycles = 0;
    so = -1;
    if((data[icm][0]-av) > 0.0) so = 1;
    
    cp = 0;
    
    for(i=0; i< N; i++)
    {
        sc = -1;
        if((data[icm][i]-av) > 0.0) sc = 1;
        if(sc != so) cp++;   // count the change points
        so = sc;
        if(cp == 2*b)
        {
            cycles++;
            cyc[cycles] = i-1;
            cp = 0;
        }
    }
    
    // create bootstrap data sets
    
    fmean = malloc(NC*sizeof(double*));
    for(i=0; i<NC; i++) fmean[i]=malloc(bsteps*sizeof(double));  // averages for each bootstrap at each temperature
    
    for(bs=0; bs< bsteps; bs++)
    {
        for(ii=0; ii< NC; ii++) fmean[ii][bs] = 0.0;
    }
    
    for(bs=0; bs< bsteps; bs++)
    {
        
        k = -1;
        do
        {
            i = (double)(cycles)*gsl_rng_uniform(seed)+1;  // pick a chunk
            for(j=(cyc[i-1]+1); j<= cyc[i]; j++)
            {
                k++;
                if(k < N)
                {
                    for(ii=0; ii< NC; ii++) fmean[ii][bs] += data[ii][j];
                }
            }
            
        }while(k < N);
        
        for(ii=0; ii< NC; ii++) fmean[ii][bs] /= (double)(N);
        
    }
    
    // means and variances of the bootstraps
    mx = malloc(NC*sizeof(double));
    
    for(ii=0; ii< NC; ii++)
    {
        mx[ii] = 0.0;
        var[ii] = 0.0;
        for(bs=0; bs< bsteps; bs++)
        {
            mx[ii] += fmean[ii][bs];
            var[ii] += fmean[ii][bs]*fmean[ii][bs];
        }
        mx[ii] /= (double)(bsteps);
        var[ii] /= (double)(bsteps);
        var[ii] -= mx[ii]*mx[ii];
        var[ii] = sqrt(var[ii]);
        
        //standardize the data
        for(bs=0; bs< bsteps; bs++)
        {
            fmean[ii][bs] = (fmean[ii][bs]-mx[ii])/var[ii];
        }
    }
    
    double **rmat, **rvar;
    
    rmat = malloc(NC*sizeof(double*));
    rvar = malloc(NC*sizeof(double*));
    for(i=0; i<NC; i++)
    {
        rmat[i]=malloc(NC*sizeof(double));
        rvar[i]=malloc(NC*sizeof(double));
    }
    
    rvar[0][0] = 0.0;
    
    // next we compute the correlation r_ij and its variance
    for(ii=0; ii< NC; ii++)
    {
        for(jj=0; jj< NC; jj++)
        {
            rmat[ii][jj] = 0.0;
            for(bs=0; bs< bsteps; bs++)    rmat[ii][jj] += fmean[ii][bs]*fmean[jj][bs];
            rmat[ii][jj] /= (double)(bsteps);
        }
    }
    
    
    for(ii=0; ii< NC; ii++)
    {
        for(jj=0; jj< NC; jj++)
        {
            rvar[ii][jj] = 0.0;
            for(bs=0; bs< bsteps; bs++)    rvar[ii][jj] += pow((fmean[ii][bs]*fmean[jj][bs]-rmat[ii][jj]),2.0);
            rvar[ii][jj] *= (double)(bsteps)/pow((double)(bsteps-1),3.0);
        }
    }
    
    x = 0.0;
    y = 0.0;
    for(ii=0; ii< NC; ii++)
    {
        for(jj=0; jj< NC; jj++)
        {
            if(ii != jj)
            {
                x += rvar[ii][jj];
                y += rmat[ii][jj]*rmat[ii][jj];
            }
        }
    }
    
    shrink = 1.0-x/y;
    if(shrink < 0.0) shrink = 0.0;
    if(shrink > 1.0) shrink = 1.0;
    
    
    for(ii=0; ii< NC; ii++)
    {
        for(jj=0; jj< NC; jj++)
        {
            cmat[ii][jj] = var[ii]*var[jj];
            if(ii!=jj)cmat[ii][jj] *= shrink*rmat[ii][jj];
        }
    }
    
    for(i=0; i<NC; i++)
    {
        free(rmat[i]);
        free(rvar[i]);
        free(fmean[i]);
    }
    free(rmat);
    free(rvar);
    free(fmean);
    free(cyc);
    free(cnts);
    free(var);
    free(mx);
    
    
}

static double log_likelihood(double *x, double *y, double **invC, gsl_spline *spline, gsl_interp_accel *acc, size_t NC)
{
    double *model = malloc(NC*sizeof(double)); //!< interpolated model of data
    double chi2 = 0.0; //!< chi-squared fit of model to data
    
    /* get interpolated value of spline model at each data point */
    for(int i=0; i<NC; i++) model[i] = gsl_spline_eval(spline,x[i],acc);
    
    /* compute chi-squared testing model against covariance matrix from chain samples */
    for(int i=0; i<NC; i++)
    {
        for(int j=0; j<NC; j++)
        {
            chi2 += (y[i]-model[i]) * invC[i][j] * (y[j]-model[j]);
        }
    }
    
    free(model);
    
    /* logL = -chi^2/2 */
    return -0.5*chi2;
}


static double log_prior(double *data, gsl_spline *spline, gsl_interp_accel *acc, int NSP)
{
    double x; //!< independent variable
    double y; //!< dependent variable

    double dx; //!< step size for data samples

    double epsilon = 1e-4; //!< step size for numerical differentiation
    double smooth = 4; //!< degree of smoothing for spline

    double x_plus;  //!< forward difference
    double y_plus;

    double x_minus; //!< backward difference
    double y_minus;
    
    double dydx; //!< first derivative
    double d2ydx2; //!< second derivative
    
    
    double logp = 0.0;
    
    //suppress large derivatives in spline model
    for(int i=1; i< NSP-1; i++)
    {
        x = data[i];
        y = gsl_spline_eval(spline,x,acc);
        
        dx = data[i]-data[i-1];
        
        x_minus = x - epsilon;
        y_minus = gsl_spline_eval(spline,x_minus,acc);
        
        x_plus = x + epsilon;
        y_plus = gsl_spline_eval(spline,x_plus,acc);
        
        dydx   = (y_plus - y_minus)/(2.0*epsilon);
        d2ydx2 = (y_minus + y_plus - 2.0*y)/(epsilon*epsilon);
        
        logp -= (smooth/(double)(NSP)) * (d2ydx2*d2ydx2*dx*dx)/(dydx*dydx);
        
        // mononicity required
        if(dydx < 0.0) logp -= 10.0;
    }
    return logp;
}



void ThermodynamicIntegration(double *temperature, double **loglikelihood, size_t NMCMC, size_t NC, double *logEvidence, double *varLogEvidence)
{
    int Nx, Ny, NSP;
    int mc, i, j, ii;
    double logLx, logLy;
    double logpx, logpy;
    double *ref, *sprd, *datax, *datay;
    int *activex, *activey;
    double alpha, H, av;
    double max, min;
    double x, y, y3, x3, *sigma;
    double trap;
    double *sdatay, *sgm, *srm;
    double *spoints, *sdatax, *tdata, *tpoints;
    double var;
    double **cov, **icov, **chains;
    double base;
    double avr = 0.0;
    
    gsl_spline *cspline = NULL;
    gsl_interp_accel *acc = NULL;
        
    //print logZchain to verify error estimates
    char filename[100];
    sprintf(filename,"evidence_chain.dat");
    FILE *zFile = fopen(filename,"w");
        
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *seed = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set (seed, time(0));
    
    FILE *integrand = fopen("integrand.dat","w");
    for(i=0; i<NC; i++)
    {
        fprintf(integrand,"%lg %lg\n",
                log10(1./temperature[i]),
                gsl_stats_mean(loglikelihood[i], 1, NMCMC));
    }
    fclose(integrand);

    
    // maximum number of spline points - we initially overcover
    NSP = 2*(int)NC-1;
    
    datax  = malloc(NC*sizeof(double));
    datay  = malloc(NC*sizeof(double));
    sigma  = malloc(NC*sizeof(double));
    chains = malloc(NC*sizeof(double *));
    icov = malloc(NC*sizeof(double *));
    cov = malloc(NC*sizeof(double *));
    for(i=0; i<NC; i++)
    {
        chains[i] = malloc(NMCMC*sizeof(double));
        icov[i] = malloc(NC*sizeof(double));
        cov[i] = malloc(NC*sizeof(double));
    }
    
    for(j=0; j< NC; j++) datax[NC-1-j] = -log(temperature[j]);
    
    double avgLogL = 0.0;
    
    for(j=0; j< NMCMC; j++)
    {
        for(i=0; i< NC; i++)
        {
            chains[NC-1-i][j] = loglikelihood[i][j];
            avgLogL += loglikelihood[i][j];
        }
    }
    
    // subtract the mean, helps reduce trapazoid error
    avgLogL /= (double)(NC*NMCMC);
    base = avgLogL*(1.0-exp(datax[0]));
    
    for(j=0; j< NMCMC; j++)
    {
        for(i=0; i< NC; i++)
        {
            chains[i][j] -= avgLogL;
        }
    }
    
    bootcorrposdef(chains, (int)NC, (int)NMCMC, cov, datay, seed);
    
    
    /*
     *
     *  Linear algebra gymnastics with the
     *  correlation matrix
     *
     */
    
    // allocate memory for GSL linear algebra
    
    gsl_matrix *matrix  = gsl_matrix_alloc (NC, NC);
    gsl_matrix *inverse = gsl_matrix_alloc (NC, NC);
    
    gsl_vector *eval = gsl_vector_alloc (NC);
    gsl_matrix *evec = gsl_matrix_alloc (NC, NC);
    
    gsl_permutation *perm = gsl_permutation_alloc(NC);
    
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (NC);
    
    
    // get eigenvalues and eigenvectors of correlation matrix (matrix is destroyed)
    for (i = 0; i < NC; i++)
    {
        sigma[i] = sqrt(cov[i][i]);
        for (j = 0; j < NC; j++)
        {
            gsl_matrix_set(matrix,i,j,cov[i][j]);
        }
    }
        
    int matrix_condition_flag=1;
    
    /*
     If the correlation matrix is singular we condition it by boosting
     the diagonal elements until the matrix is invertable
     */
    while(matrix_condition_flag)
    {
        matrix_condition_flag=0;
        
        // the eigenstuff calculator
        gsl_eigen_symmv(matrix, eval, evec, w);
        
        /*
         if the correlation matrix is up to snuff
         calculate its inverse, otherwise condition
         */
        for(i=0; i<NC; i++)
        {
            if(gsl_vector_get(eval,i) < 0.0)
            {
                fprintf(stderr,"%s:%d WARNING: singluar matrix, conditioning...\n",__FILE__,__LINE__);
                matrix_condition_flag = 1;
            }
            else
            {
                for(j=0; j<NC; j++) gsl_matrix_set(matrix,i,j,cov[i][j]);
            }
        }
        
        if(matrix_condition_flag)
        {
            for(i=0; i<NC; i++)
            {
                cov[i][i] *= 1.05;
                gsl_matrix_set(matrix,i,i,cov[i][i]);
            }
        }
        
    }
    
    // Make LU decomposition of matrix
    gsl_linalg_LU_decomp (matrix, perm, &i);
    
    // Invert the matrix
    gsl_linalg_LU_invert (matrix, perm, inverse);
    
    //get it in a form easily usable by the rest of the code
    for(i=0; i<NC; i++) for(j=0; j<NC; j++) icov[i][j] = gsl_matrix_get(inverse,i,j);
    
    // compute simple trapezoid integrand in beta
    trap = 0.0;
    for(j=1; j< NC; j++)
    {
        trap += 0.5*(exp(datax[j])-exp(datax[j-1]))*(datay[j]+datay[j-1]);
    }
    
    printf("Trapezoid integration in beta %.2f\n", trap+base);
    
    max = datax[NC-1];
    min = datax[0];
    
    ref     = malloc(NSP*sizeof(double));
    sprd    = malloc(NSP*sizeof(double));;
    spoints = malloc(NSP*sizeof(double));;
    tdata   = malloc(NSP*sizeof(double));;
    tpoints = malloc(NSP*sizeof(double));;
    sdatax  = malloc(NSP*sizeof(double));;
    sdatay  = malloc(NSP*sizeof(double));;
    sgm     = malloc(NSP*sizeof(double));;
    srm     = malloc(NSP*sizeof(double));;
    
    activex = malloc(NSP*sizeof(int));
    activey = malloc(NSP*sizeof(int));
    
    for(i=0; i< NSP; i++)
    {
        activex[i] = 1;
        activey[i] = 1;
    }
    
    /*
     Initialize spline model:
     Here we place a spline control point at,
     and bisecting, every data sample
     */
    for(i=0; i<NC; i++) //at the grid points
    {
        j = 2*i;
        sdatax[j]  = datay[i];
        spoints[j] = datax[i];
        sprd[j] = 2.0*sigma[i];
    }
    for(i=0; i<NC-1; i++) //bisecting the grid points
    {
        j=2*i+1;
        sdatax[j]  = 0.5*(datay[i] + datay[i+1]);
        spoints[j] = 0.5*(datax[i] + datax[i+1]);
        sprd[j] = 4.0*(sigma[i] + sigma[i+1]);
    }
    
    //keep copy of original spline points as reference for proposals/priors
    memcpy(ref,sdatax,NSP*sizeof(double));
        
    /*
     *
     *  Initialize likelihood for MCMC
     *
     */

    cspline = gsl_spline_alloc(gsl_interp_cspline, NSP);
    acc = gsl_interp_accel_alloc();
    gsl_spline_init(cspline,spoints,sdatax,NSP);
        
    logLx = log_likelihood(datax, datay, icov, cspline, acc, NC);
    logpx = log_prior(spoints, cspline, acc, NSP);
    
    gsl_spline_free (cspline);
    gsl_interp_accel_free (acc);
    
    
    av    = 0.0;
    var   = 0.0;
    
    /*
     *
     *  The MCMC
     *
     */
    int mcmc_steps = 10000000;  // number of MCMC steps
    int downsample = 1000; //downsample rate
    int trap_steps = 1000; //number of interpolated steps for integration
    
    double **integrand_draws = malloc(trap_steps*sizeof(double *));
    for(i=0; i<trap_steps; i++) integrand_draws[i] = malloc(mcmc_steps/downsample*sizeof(double));
    
    Nx = Ny = NSP;

    for(mc=0; mc<mcmc_steps; mc++)
    {
        logLy=0.0;
        logpy=0.0;
        
        for(i=0; i< NSP; i++)
        {
            sdatay[i] = sdatax[i];
            activey[i] = activex[i];
        }
        
        if(gsl_rng_uniform(seed) > 0.5)   // propose a dimension change
        {
            
            if(gsl_rng_uniform(seed) < 0.5)
                Ny = Nx + 1;
            else
                Ny = Nx - 1;
            
            //skip to next mcmc step if spline points are out of bounds
            if(Ny < 2 || Ny > NSP) continue;
            
            if(Ny < Nx) // propose a kill
            {
                do i = 1 + (int)(gsl_rng_uniform(seed)*(double)(NSP-2)); // pick one to kill (can't kill end points)
                while(activex[i] == 0);  // can't kill it if already dead
                
                activey[i] = 0;
            }
            else
            {
                do i = 1 + (int)(gsl_rng_uniform(seed)*(double)(NSP-2)); // pick one to add
                while(activex[i] == 1);  // can't add if already in use
                
                activey[i] = 1;
                sdatay[i] = ref[i]+10.0*sprd[i]*(1.0-2.0*gsl_rng_uniform(seed));
            }
        }
        else     // within dimension update
        {
            Ny = Nx;
            
            alpha = gsl_rng_uniform(seed);
            
            for(ii=0; ii< NSP; ii++)
            {
                // variety of jump sizes
                if(alpha > 0.8)          sdatay[ii] += sprd[ii]*gsl_ran_gaussian(seed,1);
                else if (alpha > 0.5)    sdatay[ii] += sprd[ii]*1.0e-1*gsl_ran_gaussian(seed,1);
                else if (alpha > 0.2)    sdatay[ii] += sprd[ii]*1.0e-2*gsl_ran_gaussian(seed,1);
                else                     sdatay[ii] += sprd[ii]*1.0e-3*gsl_ran_gaussian(seed,1);
            }
        }
        
        // check that the active points are in range
        for(ii=0; ii< NSP; ii++)
        {
            if(activey[ii] == 1)
            {
                if(sdatay[ii] > ref[ii]+10.0*sprd[ii]) continue; //too large, skip to next iteration
                if(sdatay[ii] < ref[ii]-10.0*sprd[ii]) continue; //too small, skip to next iteration
            }
        }
                
        // store active spline points in tpoints array to pass to interpolator
        i = 0;
        for(ii=0; ii< NSP; ii++)
        {
            if(activey[ii] == 1)  // only use active points
            {
                tpoints[i] = spoints[ii];
                tdata[i] = sdatay[ii];
                i++;
            }
        }
        
        cspline = gsl_spline_alloc(gsl_interp_cspline, Ny);
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(cspline,tpoints,tdata,Ny);
        
        logLy = log_likelihood(datax, datay, icov, cspline, acc, NC);
        logpy = log_prior(spoints, cspline, acc, NSP);
        
        gsl_spline_free (cspline);
        gsl_interp_accel_free (acc);
        
        H = logLy - logLx + logpy - logpx;
        
        alpha = log(gsl_rng_uniform(seed));
        
        if(H > alpha)
        {
            //copy proposal to current state of model
            logLx = logLy;
            logpx = logpy;
            Nx = Ny;
            
            memcpy(sdatax,sdatay,NSP*sizeof(double));
            memcpy(activex,activey,NSP*sizeof(int));
        }
        
        if(mc%downsample == 0)
        {
            i = 0;
            for(ii=0; ii< NSP; ii++)
            {
                if(activex[ii] == 1)  // only use active points
                {
                    tpoints[i] = spoints[ii];
                    tdata[i] = sdatax[ii];
                    i++;
                }
            }
            
            cspline = gsl_spline_alloc(gsl_interp_cspline, Nx);
            acc = gsl_interp_accel_alloc();
            gsl_spline_init(cspline,tpoints,tdata,Nx);
            

            /*
             Compute log evidence using trapezoid integration of
             fine-grid interpolated spline model
             */
            
            //initialize trapezoid integration at first point
            trap = 0.0;
            x3 = min;
            y3 = gsl_spline_eval(cspline,x3,acc);
            
            for(i=1; i<= trap_steps; i++)
            {
                //save current interpolated values
                x = x3;
                y = y3;
                
                //get new interplated values
                x3 = min + (max-min)*(double)(i)/(double)(trap_steps);
                y3 = gsl_spline_eval(cspline,x3,acc);
                
                //store model for diagnostics
                integrand_draws[i-1][mc/downsample] = y3;
                
                //trapezoid integration
                trap += 0.5*(exp(x3)-exp(x))*(y3+y);
            }
            
            gsl_spline_free (cspline);
            gsl_interp_accel_free (acc);

            
            avr += trap;
            var += trap*trap;
            
            fprintf(zFile,"%.12g\n",trap+base);
            
            
            
        }
    }
    fclose(zFile);
    
    x = avr/(double)(mcmc_steps/downsample);
    y = var/(double)(mcmc_steps/downsample);
    
    *logEvidence = x+base;
    *varLogEvidence = y-x*x;
    
    /* print quantiles of reconstructed integrand samples */
    FILE *integrand_quantiles = fopen("integrand_quantiles.dat","w");
    for(i=0; i<trap_steps; i++)
    {
        gsl_sort(integrand_draws[i],1,mcmc_steps/downsample);
        double int_50 = gsl_stats_median_from_sorted_data(integrand_draws[i], 1, mcmc_steps/downsample);
        double int_25 = gsl_stats_quantile_from_sorted_data(integrand_draws[i], 1, mcmc_steps/downsample, 0.25);
        double int_75 = gsl_stats_quantile_from_sorted_data(integrand_draws[i], 1, mcmc_steps/downsample, 0.75);
        double int_05 = gsl_stats_quantile_from_sorted_data(integrand_draws[i], 1, mcmc_steps/downsample, 0.05);
        double int_95 = gsl_stats_quantile_from_sorted_data(integrand_draws[i], 1, mcmc_steps/downsample, 0.95);
        x = min + (max-min)*(double)(i+1)/(double)trap_steps;
        fprintf(integrand_quantiles,"%lg %lg %lg %lg %lg %lg\n",log10(exp(x)),
                int_50+base,
                int_25+base,
                int_75+base,
                int_05+base,
                int_95+base);
    }
    fclose(integrand_quantiles);
    
    gsl_matrix_free(matrix);
    gsl_matrix_free(inverse);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_permutation_free(perm);
    gsl_eigen_symmv_free(w);
    
    
    for(i=0; i<NC; i++)
    {
        free(chains[i]);
        free(icov[i]);
        free(cov[i]);
    }
    
    free(datax);
    free(datay);
    free(sigma);
    
    free(ref);
    free(sprd);
    free(spoints);
    free(tdata);
    free(tpoints);
    free(sdatax);
    free(sdatay);
    free(sgm);
    free(srm);
    
    free(activex);
    free(activey);
    
}
