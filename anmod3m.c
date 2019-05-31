/*	anmod3m.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "anmod3m.h"

int an3()
{
   int error;
   int n;
   
   /* locations of poles */
 	double ta;   /*Pa in Fig. 3 of the paper */
 	double tb;   /*Pb in Fig. 3 of the paper */
	double rgain;     /* location of the pole closest to imaginary axis */
	double nlgain;                /* gain for the control signal */
   double zero_r;                /* Location of zeros */
   int delayn;                   /* forced delay for AN model */
    
    
   error = 0;

   /* overheader, ... will be removed before finish */
   for (n=1; n<=sound_length; n=n+1)
   {
    soundout[n]=2.0*soundin[n];
   }
   /* end of overheader */

    /*=======================================*/
	 /* parameter setup                       */
   /*=======================================*/

	 setparameter(cf, &fp1, &ta, &tb, &rgain, &nlgain, &zero_r, &delayn);
         middleear();
	 controlpath(nlgain);
 	 signalpath(ta, tb, rgain, zero_r);
	 ihczxd2001();

    printf("nlgain = %f\n", nlgain); /* This gain, due to the nonlinear control path, varies with CF */

   return error;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
   double *mxsoundin;
   double *mxsoundout;
   double *para;
   
   int mrows;
   int ncols;
   
   int mnewvec; /*these two are used to create a new vector */
   int nnewvec;
   
   int error;
   int length;
   int n;
   
/* Check for proper number of arguments. */
if(nrhs!=2)
mexErrMsgTxt("Required format: sout = an3m(cf,input);");

if(nlhs!=1)
mexErrMsgTxt("You need one output argument");

	/*  Get the input parameters. */
	if((mxGetM(prhs[0])*mxGetN(prhs[0]))!=1)
	mexErrMsgTxt("The first input must be cf only");
	para = mxGetPr(prhs[0]);
	cf = para[0];

	/*  Get the input */
	mxsoundin = mxGetPr(prhs[1]);

   /*  Get the dimensions of the matrix input in. */
   mrows = mxGetM(prhs[1]);
   ncols = mxGetN(prhs[1]);
   
   /* Check if it is a vector */
   if (mrows!=1 && ncols!=1)
   mexErrMsgTxt("The input sound has to be a vector");
   
   /* Create a new vector with index starting from 1 */
   if (mrows==1)
    {
     mnewvec=mrows; 
     nnewvec=ncols+1;
    } 
   if (ncols==1)
    {
     mnewvec=mrows+1; 
     nnewvec=ncols;
    } 
   soundin = mxGetPr(mxCreateDoubleMatrix(mnewvec, nnewvec, mxREAL));  
   soundout = mxGetPr(mxCreateDoubleMatrix(mnewvec, nnewvec, mxREAL)); 
   meout = mxGetPr(mxCreateDoubleMatrix(mnewvec, nnewvec, mxREAL));
   control_signal = mxGetPr(mxCreateDoubleMatrix(mnewvec, nnewvec, mxREAL));
   ihc_out = mxGetPr(mxCreateDoubleMatrix(mnewvec, nnewvec, mxREAL));
   sout = mxGetPr(mxCreateDoubleMatrix(mnewvec, nnewvec, mxREAL));

   /*  Set the output pointer to the output matrix. */
   plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
   mxsoundout = mxGetPr(plhs[0]);

   length = mrows*ncols;
    
   /* Get the input */
   for (n=1; n<=length; n=n+1)
   {
    soundin[n]=mxsoundin[n-1];
   }
   
   printf("cf=%f, input[1]=%f, length=%d\n",
   			cf, mxsoundin[1], length);

   /*  Call the C subroutine. */
   sound_length = length;
   error = an3();
   
   /* Get the output */
   for (n=1; n<=length; n=n+1)
   {
    mxsoundout[n-1]=sout[n];
   } 
   
   if(error) mexErrMsgTxt("Error in calling the function.");

}
