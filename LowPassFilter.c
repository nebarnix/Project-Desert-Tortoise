//Main loop taken from http://www.barrgroup.com/Embedded-Systems/How-To/Digital-Filters-FIR-IIR
//coeffs computed in matlab. Only valid for 50ksps data!

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LowPassFilter.h"

/*
 * Sample the input signal (perhaps via A/D).
 */
void LowPassFilter(float *dataStream, unsigned long nSamples)
   {
   /*
    * Insert the newest sample into an N-sample circular buffer.
    * The oldest sample in the circular buffer is overwritten.
    */
   
   int N = 27,k;
   static char firsttime = 1;
   static double h[] = {0,-0.000109655250327,0.000316665647855,0.001676613642539,-0.000476921998763,-0.007259094943430,-0.003189409960007, 0.019323908569853,0.020108432578903,-0.036803451537779,-0.072012095309267,0.053204443013450,0.305238645914910,0.439963839264128,0.305238645914910,0.053204443013450,-0.072012095309267,-0.036803451537779,0.020108432578903,0.019323908569853,-0.003189409960007, -0.007259094943430,-0.000476921998763,0.001676613642539,0.000316665647855,-0.000109655250327, 0};
   static float x[27];
   static int oldest = 0;
   double y;
   long idx, idx2=-27;
   if(firsttime == 1)
      {   
      firsttime = 0;
      memset(x,0,sizeof(float)*N);
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
      for (k = 0; k < N; k++) 
         { 
         y += h[k] * x[(oldest + k) % N]; 
         } 
      
      oldest = (oldest + 1) % N;
      
      /*
       * Output the result.
       */
      if(idx2 >= 0) //delay output so that the buffer can fill 
          dataStream[idx2] = y;
      idx2++;
      }
   }