#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MakeLPFIR.h"
/*
main()
{
double coeffs[26];
MakeLPFIR(coeffs, 26, 11000, 50000);
}*/

int MakeLPFIR(double *h, int N, double Fc, double Fs)
{
double *hd;
int n;
double T = 1.0 / Fs;
double wc = 2.0*M_PI*Fc*T;
double tou = (N-1.0)/2.0;
double wn;

hd = (double *) malloc(sizeof(double) * N);
if (hd == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }
      
// compute IDFT
for(n=0; n < N; n++)
   {
   hd[n]=(sin(wc*(n-tou))) / (M_PI*(n-tou));
    
   //check for center 
   if ((n == tou) && ((int)((N / 2)*2) != N))
      {
      hd[n] = wc / M_PI;
      //printf("%d, %d\n",(int)((N / 2)*2), N);
      }
   //printf("hd[%d]=%f\n",n,hd[n]);
   }

// apply window
for(n=0; n < N; n++)
   {
   //wn = 1;
   wn = 0.42 - 0.5 * cos((2 * M_PI * n)/(N-1)) + 0.08 * cos((4 * M_PI * n)/(N-1));
   //wn = 0.54-0.46*cos((2*M_PI*n)/(N-1));
   //wn = (1 - cos((2*M_PI*n) / (N-1)));
   h[n] = hd[n] * wn;
   }

/*for(n=0;n < N;n++)
   {
   printf("h[%d]=%f\n",n,h[n]);
   }*/
free(hd);
return n;
}