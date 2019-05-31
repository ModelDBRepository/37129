#define z_order  1

void controlpath(double nlgain)

{
   extern float cf;
   extern double tdres;
   extern double *control_signal;
   extern double *meout;
   extern int sound_length;

   int i, n;
   double PI=3.14159265358979;
   double fs_bilinear;  /* bilinear frequency =2.0/tdres */

   double x_cf;           /* location of CF on frequency map   */
   double f_shift;        /* frequency shift corresponding to 1.2mm distance along the BM */
   double wbw;            /* bandwidth of the wide-band filter */
   double gain_norm_bp;   /* normalization factor for Bandpass filter   */
   double temp01_bp, temp02_bp;  /* temporary variable for normalization*/
   mycomplex p[7];      /* poles in control space   */
   mycomplex pd[7];     /* poles in discrete domain */
   double z[7];         /* zeros in control space   */
   double zd[7];        /* zeros in discrete domain */
   double kkd;          /* a variable used in gain control */
   double pda[6], pdb[6];
   double goutput[10][4];
   double ginput[10][4];
   double dy;           /* output of digital filter */
   double tempdouble01; /* temporary variables      */
   double tempdouble02;
   double wbw2pi;       /* this is equal to 2*PI*wbw */

   double preal;
   double pimg;

   double control_nl=0.0;

   double p_corner;     /* parameters for the first nonlinearity */
   double p_strength;
   double p_slope;
   double splx;
   double xx;           /* output of the first nonlinearity */
   double acp;          /* parameters for Zhang's(2001) first NL */
   double bcp;
   double ccp;


   double potential;    /* parameters for the second nonlinearity */
   double s0, s1, x0, x1;
   double shift;
   double asym;

   double conlp[5][2];
   double bw_conlp;
   double bw_conlp_2pi;

   /* parameters for the feedback LP filter */
   double fblp[2][2];   /* 1st order LP */
   double bw_fblp;      /* bandwidth of the feedback LP */
   double bw_fblp_2pi;  /* bandwidth multiplied by 2pi  */


   /* bilinear transformation frequency */
   fs_bilinear = 2.0/tdres;

    //%=========================================================
    //% band-pass filter
    //%=========================================================

    //frequency map
    x_cf=11.9*log10(0.8+cf/456);
    //%the peak of the suppression curve is shifted by 1.2mm for all CF
    f_shift=(pow(10,((x_cf+1.2)/11.9))-0.8)*456-cf;

     // wbw controls the bandwidth of the wide-bandpass filter
      wbw=cf/4.0;               

     /* initial locations of poles and zeros */
     p[1].realpart= -2*PI*wbw;
     p[1].imgpart = 2*PI*(cf+f_shift);
     p[2].realpart= -2*PI*wbw;
     p[2].imgpart = -2*PI*(cf+f_shift);

     p[3]=p[1];
     p[4]=p[2];

     p[5]=p[1];
     p[6]=p[2];

     z[1]=0;         /* only one zero is used */
     z[2]=0; 
     z[3]=0; 
     z[4]=0; 
     z[5]=0; 

   for (i=1; i<=8; i=i+1)
     {
      ginput[i][1]=0.0;
      ginput[i][2]=0.0;
      ginput[i][3]=0.0;
      goutput[i][1]=0.0;
      goutput[i][2]=0.0;
      goutput[i][3]=0.0;
     }

   /* setup for the first NL */
       acp=100; 
       bcp=2.5; 
       ccp=0.60; 

   /* setup for the control LP */
   bw_conlp=800;
   bw_conlp_2pi=bw_conlp*2*PI;

   conlp[0][0]=0.0;
   conlp[0][1]=0.0;
   conlp[1][0]=0.0;
   conlp[1][1]=0.0;
   conlp[2][0]=0.0;
   conlp[2][1]=0.0;
   conlp[3][0]=0.0;
   conlp[3][1]=0.0;

   /* setup for the feedback LP */
   bw_fblp = 500;
   bw_fblp_2pi = bw_fblp*2*PI;
   fblp[0][0]=0.0;
   fblp[0][1]=0.0;
   fblp[1][0]=0.0;
   fblp[1][1]=0.0;

   for (n=1; n<=sound_length; n=n+1)  
    {
    /*=========================================================*/
    /* band-pass filter                                        */
    /*=========================================================*/

     ginput[1][3]=ginput[1][2];
     ginput[1][2]=ginput[1][1];
     ginput[1][1]=meout[n];       // meout is X0 in paper Fig. 1.


        /* this part needs to be modified to be a loop,         */
        /* if you want different poles at different locations. */
        /* wbw and wbw2pi are always positive                  */
        wbw2pi=-(p[1].realpart - control_nl);

        /* normalize the gain at cf */
        wbw=wbw2pi/2.0/PI;
        temp01_bp=sqrt(wbw*wbw + f_shift*f_shift);
        temp02_bp=sqrt( (2*cf+f_shift)*(2*cf+f_shift) + (wbw*wbw) );
        gain_norm_bp=2.0*PI*temp01_bp*2.0*PI*(temp02_bp);
        gain_norm_bp=gain_norm_bp*gain_norm_bp*gain_norm_bp;

        /* normalization factor related to zero(s)     */

        gain_norm_bp=gain_norm_bp/
                 pow( sqrt(2*PI*z[1]*2*PI*z[1] + 2*PI*cf*2*PI*cf), z_order );

        kkd=1.0;

     for (i=1; i<=3; i=i+1)          /* 3*2 poles */
      {
       preal=p[i*2].realpart- control_nl;
       pimg=p[i*2].imgpart;

       tempdouble01=(fs_bilinear-preal);
       tempdouble01=tempdouble01*tempdouble01 + pimg*pimg;

       dy=(ginput[i][1] + 2*ginput[i][2] + ginput[i][3]);
       dy=dy+2*goutput[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

       dy=dy-goutput[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);

       dy=dy/tempdouble01;

       ginput[i+1][3]=goutput[i][2];
       ginput[i+1][2]=goutput[i][1];
       ginput[i+1][1]=dy;

       goutput[i][2]=goutput[i][1];
       goutput[i][1]=dy;

      }


      for (i=4; i<=3+z_order; i=i+1)            /* zeros */
      {
       goutput[i][1]= ginput[i][1]*(fs_bilinear-z[i-3])
                      -ginput[i][2]*(fs_bilinear+z[i-3])
                      -goutput[i][1];
       ginput[i+1][2]=ginput[i+1][1];
       ginput[i+1][1]=goutput[i][1];          

      }

     dy=goutput[3+z_order][1];

     dy=dy*kkd*gain_norm_bp;         /* don't forget the gain term */

     /* this dy is X1 in paper fig. 1 */

     /*=========================================================*/
     /*  the first nonlinearity                                 */
     /*=========================================================*/

     /* the first NL of Zhang(2001) */
      if (dy>=0)
       {
        xx=bcp*log(1.0+acp*pow(dy, ccp));
       }
      else
       {
        xx=-bcp*log(1.0+acp*pow(-dy, ccp));
       }
    
     /* this xx is X2 in paper fig. 1 */

     /*==========================================================*/
     /*  the second nonlinearity                                 */
     /*==========================================================*/

        asym=7.0;
        s0 = 8.0;
        x1 = 5.0;
        s1 = 3.0;    
        shift = 1.0/(1.0+asym);

        x0 = s0*log((1.0/shift-1)/(1+exp(x1/s1)));

        potential = 1.0/(1.0+exp(-(xx-x0)/s0)*(1.0+exp(-(xx-x1)/s1)))-shift;

       /* this potential is Y in paper figure 1 */

       potential=potential*nlgain;

       //%==========================================================
       //%   low pass again w/ cut-off freq = 800 Hz
       //%==========================================================

       conlp[0][0]=conlp[0][1];
       conlp[0][1]=potential;

       for (i=1; i<=3; i=i+1)
        {
         conlp[i][1]= ( conlp[i-1][1] + conlp[i-1][0]
                    + conlp[i][0]*(fs_bilinear - bw_conlp_2pi))
                          /(fs_bilinear + bw_conlp_2pi);
        }
       for (i=1; i<=3; i=i+1)
        {
         conlp[i][0] = conlp[i][1];
        }
       control_signal[n]=conlp[3][1]*800*2*PI*800*2*PI*800*2*PI*1.5;

     /*===================================================================*/
     /* feedback low-pass filter                                          */
     /* parameters for the feedback LP filter                             */

       fblp[0][0]=fblp[0][1];
       fblp[0][1]=control_signal[n];

       for (i=1; i<=1; i=i+1)
        {
         fblp[i][1]= ( fblp[i-1][1] + fblp[i-1][0]
                    + fblp[i][0]*(fs_bilinear - bw_fblp_2pi))
                          /(fs_bilinear + bw_fblp_2pi);
        }
       for (i=1; i<=1; i=i+1)
        {
         fblp[i][0] = fblp[i][1];
        }
       control_nl=fblp[1][1]*500*2*PI *10;

  }    // =================================================end of time loop

}



