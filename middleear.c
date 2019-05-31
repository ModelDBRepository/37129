/* middleear.c */
/* a simple middle ear model */

/* July 04, 2002 */
/* the output of this program now is meout[] */

/* Sept 26, 2002 */
/* simplify the program */

int middleear()
{
 int error_number=1;

 int pole_order=4;           // two pairs of poles
 int half_pole_order=2;
 int zero_order=2;           // two zeros

 extern double *soundin;
 extern double tdres;
 extern double PI;
 extern double *meout;
 extern int sound_length;
 
 double gain_norm;

 double fs_bilinear; //bilinear transformation frequency
 double kkd;
 double bpinput[4][4];
 double bpoutput[4][4];
 double zerobp;

 double preal;
 double pimg;

 double a1;
 double a2;
 double b1[3];
 double b2[3];
 double dy;
 
 mycomplex p_11, p_12, p_21, p_22; // poles in s plain
 mycomplex p[5];
 
 int n, i, j;

   //%========== locations of poles ===================%

   p_11.realpart=-250*2*PI;
   p_11.imgpart=400*2*PI;

   p_21.realpart = -2000*2*PI;
   p_21.imgpart = 6000*2*PI;

   p_12 = myconj(p_11);
   p_22 = myconj(p_21);

   p[1]=p_11;
   p[2]=p_12;
   p[3]=p_21;
   p[4]=p_22;

   fs_bilinear = 2.0/tdres;

   kkd=1.0;
   for (i=1; i<=half_pole_order; i=i+1)
    {
     kkd=kkd/((fs_bilinear-p[i*2].realpart)*(fs_bilinear-p[i*2].realpart)
              +p[i*2].imgpart*p[i*2].imgpart);
    }

   /* ========= setup zeros =================*/
   zerobp=-200;
   for (i=1; i<=zero_order; i=i+1)
    {
     kkd=kkd*(fs_bilinear - zerobp);
    }

  /*========= initialize some variables======*/
  for (i=1; i<=3; i=i+1)
   {
    for (j=1; j<=3; j=j+1)
     {
      bpinput[i][j]=0.0;
      bpoutput[i][j]=0.0;
     }
   }   

  /*======= normalize gain at 1000Hz========== */
  gain_norm =1.0;
  for (i=1; i<=pole_order; i=i+1)
   {
    gain_norm=gain_norm*sqrt(p[i].realpart*p[i].realpart
                         +(p[i].imgpart-1000*2*PI)*(p[i].imgpart-1000*2*PI));
   }
  for (i=1; i<=zero_order; i=i+1)
   {
    gain_norm=gain_norm/sqrt(2*PI*zerobp*2*PI*zerobp+1000.0*2*PI*1000.0*2*PI);
   }

  /* filter coefficients */

  a1=(1-(fs_bilinear+zerobp)/(fs_bilinear-zerobp));
  a2=- (fs_bilinear+zerobp)/(fs_bilinear-zerobp);
  for (i=1; i<=half_pole_order; i=i+1)
   {
    b1[i]=2*(fs_bilinear*fs_bilinear
              -p[i*2].realpart*p[i*2].realpart-p[i*2].imgpart*p[i*2].imgpart)
           /((fs_bilinear-p[i*2].realpart)*(fs_bilinear-p[i*2].realpart)
              +p[i*2].imgpart*p[i*2].imgpart);
    b2[i]=-((fs_bilinear+p[i*2].realpart)*(fs_bilinear+p[i*2].realpart)
              +p[i*2].imgpart*p[i*2].imgpart)
           /((fs_bilinear-p[i*2].realpart)*(fs_bilinear-p[i*2].realpart)
              +p[i*2].imgpart*p[i*2].imgpart);
   }

   
  for (n=1; n<=sound_length ; n=n+1)
   {
   bpinput[1][3]=bpinput[1][2];
   bpinput[1][2]=bpinput[1][1];
   bpinput[1][1]=soundin[n];

   for (i=1; i<=half_pole_order; i=i+1)
    {
     dy = bpinput[i][1]+a1*bpinput[i][2]+a2*bpinput[i][3];
     dy = dy + b1[i]*bpoutput[i][1] + b2[i]*bpoutput[i][2];

     bpinput[i+1][3]=bpoutput[i][2];
     bpinput[i+1][2]=bpoutput[i][1];
     bpinput[i+1][1]=dy;

     bpoutput[i][2]=bpoutput[i][1];
     bpoutput[i][1]=dy;
    }

   meout[n]=dy*kkd*gain_norm;

  }
   return (error_number);
}   
