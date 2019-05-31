// ihczxd2001.cpp
// from zxd, 2001
//Aug 3, 2001

// sout_average use 10ms to 45ms window now
// NOv 17, 2001

// use universal version of zxd2001, all checked
// ifspike is an extern now
// last version is in /0909/
// Aug28, 2002

void ihczxd2001()
{
  extern double PI;
  extern double *soundout;
  extern double tdres;
  extern double *sout;
  extern double *ihc_out;
  extern float cf;
  extern int ifspike;      // 1--> generate spike
                           // 0--> use sout
  extern int sound_length;

  //extern double sout_average;

  //double sout_sum;   //used for the calculation of sout_average;
  //double tscale;     //used for the calculation of sout_average;

  int i, j;           // loop counters
  double tempdouble;
  double tempA;


  int ihc_lp_order;   // order of the low-pass filter in IHC
  double ihc_lp_gain; // gain of the low-pass filter in IHC
  double ihc_lp_Fc;   // cut-off freq of the low-pass filter in IHC
  double ihc_lp[15];  // tmp for ihc lowpass
  double ihc_c1LP;
  double ihc_c2LP;
  double ihc_lp_c;
  double hcl[15];

  double A0 = 0.1;
  double B = 2000.;
  double C = 1.74;
  double D = 6.87e-9;             
  double dtemp;

  double Ass;
  double PTS = 8.627;
  double Aon;
  double Ar_over_Ast = 6.;
  double Ar;
  double Ast;
  double Pimax = 0.6;
  double spont=50.0;
  double Kcf;
  double Vsat;
  double pst;
  double psl;
  double p1;
  double p2;

  double Prest;
  double CG;
  double gamma1;
  double gamma2;
  double kappa1;
  double kappa2;

  double tauR = 0.002; //human
  double tauST = 0.060;

  double VI0;
  double VI1;
  double VI;

  double alpha;
  double beta;
  double theta1;
  double theta2;
  double theta3;

  double PL;
  double PG;
  double VL;

  double Cirest;
  double CLrest;

  double Vihc;
  double PPI;
  double PPIlast;

  double CInow;
  double CLnow;
  double CIlast;
  double CLlast;

  FILE *psave02;

  if (ifspike==1) Ass = 350;
  else Ass=130; /* read Frank's cmpa.c */

  for (i = 1; i<=sound_length; i++)
  {
    /*/begin Vsp -> Vihc */
    tempdouble = soundout[i];
    if(tempdouble>=0)
    {
      tempA = A0;
    }
    else
    {
      dtemp = pow(-tempdouble,C);
      tempA = -A0*(dtemp+D)/(3*dtemp+D);
    };
    ihc_out[i] = tempA*log(fabs(tempdouble)*B+1.0);

  };

  //-----------------------------------------------------------------------
  //      low pass filter
  //-----------------------------------------------------------------------
   // initialize the parameters of the low-pass
   // from zxd's synapse.c and filters.c
   ihc_lp_order = 7;
   ihc_lp_gain = 1.0;
   ihc_lp_Fc =3800.0;        //this would be 3800Hz if for cat.  4500 for human

   ihc_lp_c = 2.0/tdres;
   ihc_c1LP = ( ihc_lp_c - 2* PI * ihc_lp_Fc ) / ( ihc_lp_c + 2* PI * ihc_lp_Fc );
   ihc_c2LP = 2 * PI * ihc_lp_Fc / (2 * PI * ihc_lp_Fc + ihc_lp_c);

    for (j=1; j<=ihc_lp_order + 1; j=j+1)
       hcl[j]=0.0;

    for(i=1; i<=sound_length; i=i+1)
      {
       /*hc[0] = in[loopSig]*gain;

    for(loopLP=0; loopLP<pOrder; loopLP++)
      hc[loopLP+1] = c1LP * hcl[loopLP+1] + c2LP*(hc[loopLP]+hcl[loopLP]);

    for(loopLP=0; loopLP<=pOrder;loopLP++) hcl[loopLP] = hc[loopLP];
    out[loopSig] = hc[pOrder];   */

       // the above part is copied from zxd
       ihc_lp[1]=ihc_out[i]*ihc_lp_gain;

       for (j=1; j<=ihc_lp_order; j=j+1)
        ihc_lp[j+1] = ihc_c1LP * hcl[j+1] + ihc_c2LP * (ihc_lp[j]+hcl[j]);

       for (j=1; j<=ihc_lp_order+1; j=j+1)
        hcl[j] = ihc_lp[j];

       ihc_out[i]=ihc_lp[ihc_lp_order+1];

       };

  //=======================================================================

  /*/begin the Synapse dynamic */
  Aon = PTS * Ass;
  Ar = (Aon - Ass) * (Ar_over_Ast)/(1. + Ar_over_Ast);  /*/%???? Ar on both sides */
  Ast = Aon - Ass - Ar;
  Prest = Pimax * spont/Aon;
  CG = spont * (Aon - spont)/(Aon * Prest*(1. - spont/Ass));
  gamma1 = CG/spont;
  gamma2 = CG/Ass;

  kappa1 = -( 1. /tauR);
  kappa2 = -( 1. /tauST);

  VI0 = (1.-Pimax/Prest)/
     (gamma1*(Ar*(kappa1-kappa2)/CG/Pimax+kappa2/Prest/gamma1-kappa2/Pimax/gamma2));
  VI1 = (1.-Pimax/Prest)/
     (gamma1*(Ast*(kappa2-kappa1)/CG/Pimax+kappa1/Prest/gamma1-kappa1/Pimax/gamma2));
  VI = (VI0 + VI1)/2.;

  alpha = gamma2/(kappa1 * kappa2);
  beta = -(kappa1 + kappa2) * alpha;
  theta1 = alpha * Pimax/VI;
  theta2 = VI/Pimax;
  theta3 = gamma2 - 1./Pimax;
  PL = (((beta - theta2 * theta3)/theta1) - 1.)*Pimax;
  PG = 1./(theta3 - 1./PL);
  VL = theta1 * PL * PG;
  Cirest = spont/Prest;
  CLrest = Cirest * (Prest + PL)/PL;

  PPIlast = Prest;
  CLlast = CLrest;
  CIlast = Cirest;


  /* from Frank's cmpa.c for species==9 universal */
  Kcf = 2.0+3.0*log10(cf/1000);
  if(Kcf<1.5) Kcf = 1.5;
  Vsat = Pimax*Kcf*20.0*(1+ spont)/(5+ spont);

  tempdouble = log(2)* Vsat/ Prest;
  //if(tempdouble<400)
  if(tempdouble<100)         //zxd used 400, I think 100 is enough and safer.
   {
     pst = log(exp(tempdouble)-1);
   }
  else
   {
    pst = tempdouble;
   }

  psl = Prest*pst/log(2);
  p2 = pst;
  p1 = psl/pst;

  for(i = 1; i<=sound_length; i++)
  {
    Vihc = ihc_out[i];
    tempdouble = p2 * Vihc;
    if (tempdouble < 100)
       {
        PPI = p1 * log(1. + exp(p2 * Vihc)); /*/ soft-rectifier */
       }
    else
       {
        PPI = p1 * tempdouble;
       }

    ihc_out[i] = PPI;

  };


  for(i = 1; i<=sound_length; i++)
  {
    CInow = CIlast + (tdres/VI)*((-PPIlast*CIlast)+PL*(CLlast - CIlast));
    CLnow = CLlast + (tdres/VL)*(-PL*(CLlast-CIlast)+PG*(CG - CLlast));
    PPIlast = ihc_out[i];
    CIlast = CInow;
    CLlast = CLnow;
    sout[i] = CInow*PPIlast;
  };


};



