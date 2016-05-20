//Main loop taken from http://www.barrgroup.com/Embedded-Systems/How-To/Digital-Filters-FIR-IIR
//coeffs computed in matlab. Only valid for 50ksps data!

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LowPassFilter.h"

/*
 * Sample the input signal (perhaps via A/D).
 */
void LowPassFilter(float *dataStream, unsigned long nSamples, double *filterCoeffs, int N)
   {
   /*
    * Insert the newest sample into an N-sample circular buffer.
    * The oldest sample in the circular buffer is overwritten.
    */
   
   unsigned int idx2;
   static char firsttime = 1;
   //static double h[] = {0,-0.000109655250327,0.000316665647855,0.001676613642539,-0.000476921998763,-0.007259094943430,-0.003189409960007, 0.019323908569853,0.020108432578903,-0.036803451537779,-0.072012095309267,0.053204443013450,0.305238645914910,0.439963839264128,0.305238645914910,0.053204443013450,-0.072012095309267,-0.036803451537779,0.020108432578903,0.019323908569853,-0.003189409960007, -0.007259094943430,-0.000476921998763,0.001676613642539,0.000316665647855,-0.000109655250327, 0};
   static double *x; //circular buffer
   static int oldest = 0;
   double y;
   long idx;
   //long idx2 = -N; //BUG! -N delay should only be for FIRST GO AROUND!
   if(firsttime == 1)
      {   
      firsttime = 0;
      x = (double *) malloc(sizeof(double) * N);
      if (x == NULL)
            {
            printf("Error in malloc\n");
            exit(1);
            }
      memset(x,0,sizeof(double)*N); //zero out the circular buffer
      }
   // printf("LPF\n");
   
   for(idx = 0; idx < nSamples; idx++)
      {
      x[oldest] = dataStream[idx]; 
      
      /*
       * Multiply the last N inputs by the appropriate coefficients.
       * Their sum is the current output.
       */      
      y = 0;
      for (idx2 = 0; idx2 < N; idx2++) 
         { 
         y += filterCoeffs[idx2] * x[(oldest + idx2) % N];         
         }      
      
      /* visualize buffers
      printf("circ: ");
      for(idx2 = 0; idx2 < N; idx2++)
         {         
         printf("%0.2f,", x[(oldest + 1 + idx2) % N]);
         }
      printf("\n");
      
       //spit out what is in the buffer
      printf("raw:  ");
      for(idx2 = 0; idx2 < N; idx2++)
         {         
         printf("%0.2f,", x[idx2]);
         }
      printf("\n");
      */
      
      /*
       * Output the result.
       */
      //if(idx2 >= 0) //delay output so that the buffer can fill 
      //   dataStream[idx2] = y;
      //idx2++;
      dataStream[idx] = y; //you have to live with the delay, otherwise, how to handle multiple calls?
      //if(idx > 30)
      //   exit(1);
      oldest = (oldest + 1) % N;
      }
   }
   
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
