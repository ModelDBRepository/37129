/* spath_oct2001.cpp */
/* this version intends to remove complex.h because complex calculation */
/* does not seem to be supported by ANSI C     ............... Oct 28, 2001 */

/* July 04, 2002 */
/* name is changed to signalpath.c */

/*Feb 23, 2003 */
/* aa0 is removed */
/* reduced some workload for linear model (don't update pole locations) */

void signalpath(double ta, double tb, double rgain, double zero_r)

  {
   int n, i;

   int order_of_pole;
   int half_order_pole;
   int order_of_zero;

   extern int control_type;
   extern double *meout;
   extern double *soundout;
   extern double *control_signal;
   extern double tdres;
   extern double PI;
   extern float cf;
   extern float fp1;
   extern int sound_length;
   
   double aa;
   double zeroa; /* location of the zeros */
   double fs_bilinear; /* bilinear transformation frequency */
   

   double bfp; /* best frequency location on imaginary axis */

   mycomplex poles[21];   /* can use up to 10 pairs of poles */

   double gain_norm;

   double dy;

   double bpinput[22][4];
   double bpoutput[22][4];
   double preal;
   double pimg;
   
   double tempdouble01;

     
   fs_bilinear = 2.0/tdres;

   /*================ setup the locations of poles and zeros =======*/

   order_of_pole=20;             /* 5 pairs of poles */
   half_order_pole=10;
   order_of_zero=10;

   poles[1].realpart=-rgain;
   poles[1].imgpart=fp1*2*PI;

   poles[5].realpart = poles[1].realpart - ta;
   poles[5].imgpart = poles[1].imgpart - tb;

   poles[3].realpart = (poles[1].realpart + poles[5].realpart) * 0.5;
   poles[3].imgpart =  (poles[1].imgpart +  poles[5].imgpart) * 0.5;

   poles[2] = myconj(poles[1]);
   poles[4] = myconj(poles[3]);
   poles[6] = myconj(poles[5]);

   poles[7] = poles[1];
   poles[8] = poles[2];
   poles[9] = poles[5];
   poles[10]= poles[6];

   for (i=1; i<=10; i=i+1)
    {
     poles[i+10]=poles[i];
    }

   zeroa=-zero_r;

   /*===================== normalize the gain =====================*/

   bfp=2*PI*cf;

   gain_norm=1.0;

   for (n=1; n<=order_of_pole; n=n+1)
    {
     tempdouble01= bfp-poles[n].imgpart;
     gain_norm= gain_norm *
             (tempdouble01*tempdouble01 + poles[n].realpart*poles[n].realpart);
    }


   gain_norm=sqrt(gain_norm);

   gain_norm=gain_norm/pow((sqrt(bfp*bfp+zeroa*zeroa)), order_of_zero);


   for (i=1; i<=half_order_pole + order_of_zero+1; i=i+1)
    {
      bpinput[i][1]=0.0;
      bpinput[i][2]=0.0;
      bpinput[i][3]=0.0;
      bpoutput[i][1]=0.0;
      bpoutput[i][2]=0.0;
      bpoutput[i][3]=0.0;
    }

   /*%==================================================  */
   /*%      time loop begins here                         */
   /*%==================================================  */
   for (n=1; n<=sound_length ; n=n+1)
   {
         aa=-rgain - control_signal[n];
         if (aa>=0)
           {
stayhere:  aa=100;         // you can choose to print an error msg here
            goto stayhere;
           }

       //%========== update pole locations ===================%
       poles[1].realpart=aa;
       poles[1].imgpart=fp1*2*PI;

       poles[5].realpart = poles[1].realpart - ta;
       poles[5].imgpart = poles[1].imgpart - tb;

       poles[3].realpart = (poles[1].realpart + poles[5].realpart) * 0.5;
       poles[3].imgpart =  (poles[1].imgpart +  poles[5].imgpart) * 0.5;

       poles[2] = myconj(poles[1]);
       poles[4] = myconj(poles[3]);
       poles[6] = myconj(poles[5]);

       poles[7] = poles[1];
       poles[8] = poles[2];
       poles[9] = poles[5];
       poles[10]= poles[6];

       for (i=1; i<=10; i=i+1)
        {
         poles[i+10]=poles[i];
        }

     bpinput[1][3]=bpinput[1][2];
     bpinput[1][2]=bpinput[1][1];
     bpinput[1][1]=meout[n];

     for (i=1; i<=half_order_pole; i=i+1)          // 10*2 poles
      {
       preal=poles[i*2].realpart;
       pimg=poles[i*2].imgpart;

       tempdouble01=(fs_bilinear-preal);
       tempdouble01=tempdouble01*tempdouble01 + pimg*pimg;

       dy=(  bpinput[i][1] +
            (1-(fs_bilinear+zeroa)/(fs_bilinear-zeroa))*bpinput[i][2]
             - (fs_bilinear+zeroa)/(fs_bilinear-zeroa)*bpinput[i][3] );

       dy=dy+2*bpoutput[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

       dy=dy-bpoutput[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);

       dy=dy/tempdouble01;

       bpinput[i+1][3]=bpoutput[i][2];
       bpinput[i+1][2]=bpoutput[i][1];
       bpinput[i+1][1]=dy;

       bpoutput[i][2]=bpoutput[i][1];
       bpoutput[i][1]=dy;

      }

     dy=bpoutput[half_order_pole][1]*pow((fs_bilinear-zeroa),10);

     dy=dy*gain_norm;         /* don't forget the gain term */
     soundout[n]=dy;

   /*each loop below is for a pair of poles and one zero */
   /*and for each loop, a zero at -1 (-infinite in control space) is added*/
   /* so that the order of zeros is same as the order of poles */

   soundout[n]=soundout[n]/3;
   /* signal path output is divided by 3 to increase the threshold at CF */

   } /*end of time loop */

 }


