#define usemiddleear 1         /* 1-> include middle ear model; 0-> don't include middle ear */
#include "setparameter.c"
#include "mycomplex.h"
#include "middleear.c"
#include "controlpath.c"
#include "signalpath.c"
#include "ihczxd2001.c"

#ifdef usespike
 #define uspike 1              /* use spike output */
#else
 #define uspike 0              /* do not use spike output */
#endif
int    ifspike=uspike;

double tdres=2e-5;                  /* sampling time (20ns, by default) */
float  cf;
float  fp1;
int    control_type=1;              /* 0 ->without suppression */
                                    /* 1 ->with suppression    */
double PI=3.14159265358979;
double *soundin;
double *meout;
double *soundout;
double *control_signal;
double *ihc_out;
double *sout;

int sound_length;

