/* setparameter.cpp */
/* June 07, 2002    */

int setparameter( float  cf, float *pfp1,
                  double *pta, double *ptb, double *prgain, double *pnlgain,
                  double *pzero_r, int *pdelay)
{
 int error_number=1;

 double rgain80;
 double average_control=0.3357;

 *pfp1=1.0854*cf-106.0034;
 *pta=     pow(10, log10(cf)*1.0230 + 0.1607);
 *ptb=     pow(10, log10(cf)*1.4292 - 1.1550) - 1000;
 rgain80=  pow(10, log10(cf)*0.5732 + 1.5220);
 *prgain=  pow(10, log10(cf)*0.4 + 1.9);   /* 75% */

 *pnlgain= (rgain80 - *prgain)/average_control;
 *pzero_r= pow(10, log10(cf)*1.5-0.9 );        
 *pdelay= 25; /* 0.5ms delay for all CFs */

 return(error_number);
}
