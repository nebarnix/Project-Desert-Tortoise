//Main loop taken from http://www.barrgroup.com/Embedded-Systems/How-To/Digital-Filters-FIR-IIR

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include "LowPassFilter.h"

/*
 * Low Pass Filter with zero stuffing interpolator 
 */
 
void LowPassFilterInterp(DECIMAL_TYPE *dataStreamInTime, DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamOut, DECIMAL_TYPE *dataStreamOutTime, unsigned long nSamples, DECIMAL_TYPE *filterCoeffs, int N, int interpFactor)
   {
   /*
   * Insert the newest sample into an N-sample circular buffer.
   * The oldest sample in the circular buffer is overwritten.
   */
   
   unsigned int idxF;
   static char firsttime = 1;
   //static double h[] = {0,-0.000109655250327,0.000316665647855,0.001676613642539,-0.000476921998763,-0.007259094943430,-0.003189409960007, 0.019323908569853,0.020108432578903,-0.036803451537779,-0.072012095309267,0.053204443013450,0.305238645914910,0.439963839264128,0.305238645914910,0.053204443013450,-0.072012095309267,-0.036803451537779,0.020108432578903,0.019323908569853,-0.003189409960007, -0.007259094943430,-0.000476921998763,0.001676613642539,0.000316665647855,-0.000109655250327, 0};
   static DECIMAL_TYPE *x; //circular buffer
   static int oldest = 0;
   DECIMAL_TYPE y; //accumulator
   unsigned long idxO, idxI=0;
   static unsigned char interpCounter=0;
   unsigned long nSamplesInterp = nSamples*interpFactor;
   
   if(firsttime == 1)
      {   
      firsttime = 0;      
      x = (DECIMAL_TYPE *) malloc(sizeof(DECIMAL_TYPE) * N);
      if (x == NULL)
         {
         printf("Error in malloc\n");
         exit(1);
         }
      memset(x,0,sizeof(DECIMAL_TYPE)*N); //zero out the circular buffer
      }
   
   //for(idxO = 0; idxO < (nSamples*interpFactor); idxO++,interpCounter++)
   for(idxO = 0; idxO < nSamplesInterp; idxO++,interpCounter++)
      { 
      if((interpCounter % interpFactor) == 0) //every for loop ALAWYS starts here, with a real data input
         {
         x[oldest] = dataStreamIn[idxI++]; //stuff the current sample into the buffer and push the oldest one out
         interpCounter = 0;
         }
      else
         x[oldest] = 0; //zero stuff the rest
      
      /*
      * Multiply the last N inputs by the appropriate coefficients.
      * Their sum is the current output.
      */      
      
      for (y = 0, idxF = 0; idxF < N; idxF+=interpFactor) //Skip over zero stuffed entries -- should speed things up nicely!
      //for (y = 0, idxF = 0; idxF < N; idxF++) 
         { 
         //y += filterCoeffs[idxF] * x[(oldest + idxF+1) % N]; //follow the X buffer
         y += filterCoeffs[(N-(oldest - idxF+1)) % N] * x[idxF]; //follow the LPF buffer
         //printf("[%d]%f*[%d]%f\t[%d]%f*[%d]%f \n",idxF,filterCoeffs[idxF],(oldest + idxF+1) % N,x[(oldest + idxF+1) % N],N-(oldest - idxF+1) % N,filterCoeffs[N-(oldest - idxF+1) % N],idxF,x[idxF]);
         }  
      //printf("-------\n");
      
      dataStreamOut[idxO] = y; 
      dataStreamOutTime[idxO] = dataStreamInTime[idxI];  
      oldest = (oldest + 1) % N;   
      }
   }

/*
 * Low Pass Filter 
 */
void LowPassFilter(DECIMAL_TYPE *dataStream, unsigned long nSamples, DECIMAL_TYPE *filterCoeffs, int N)
   {
   /*
    * Insert the newest sample into an N-sample circular buffer.
    * The oldest sample in the circular buffer is overwritten.
    */
   
   unsigned int idxF;
   static char firsttime = 1;
   //static double h[] = {0,-0.000109655250327,0.000316665647855,0.001676613642539,-0.000476921998763,-0.007259094943430,-0.003189409960007, 0.019323908569853,0.020108432578903,-0.036803451537779,-0.072012095309267,0.053204443013450,0.305238645914910,0.439963839264128,0.305238645914910,0.053204443013450,-0.072012095309267,-0.036803451537779,0.020108432578903,0.019323908569853,-0.003189409960007, -0.007259094943430,-0.000476921998763,0.001676613642539,0.000316665647855,-0.000109655250327, 0};
   static DECIMAL_TYPE *x; //circular buffer
   static int oldest = 0;
   DECIMAL_TYPE y; //accumulator
   long idxO;
   //long idxF = -N; //BUG! -N delay should only be for FIRST GO AROUND!
   if(firsttime == 1)
      {   
      firsttime = 0;
      x = (DECIMAL_TYPE *) malloc(sizeof(DECIMAL_TYPE) * N);
      if (x == NULL)
            {
            printf("Error in malloc\n");
            exit(1);
            }
      memset(x,0,sizeof(DECIMAL_TYPE)*N); //zero out the circular buffer
      }
   // printf("LPF\n");
   
   for(idxO = 0; idxO < nSamples; idxO++)
      {
      x[oldest] = dataStream[idxO]; 
      
      /*
       * Multiply the last N inputs by the appropriate coefficients.
       * Their sum is the current output.
       */      
      y = 0;
      for (idxF = 0; idxF < N; idxF++) 
         { 
         y += filterCoeffs[idxF] * x[(oldest + idxF+1) % N];         
         }           
      
      /*
       * Output the result.
       */

      dataStream[idxO] = y; //you have to live with the delay, otherwise, how to handle multiple calls?
      oldest = (oldest + 1) % N;
      }
   }
   
int MakeLPFIR(DECIMAL_TYPE *h, int N, DECIMAL_TYPE Fc, DECIMAL_TYPE Fs, int interpFactor)
   {
   DECIMAL_TYPE *hd;
   int n;
   DECIMAL_TYPE T = 1.0 / Fs;
   DECIMAL_TYPE wc = 2.0*M_PI*Fc*T;
   DECIMAL_TYPE tou = (N-1.0)/2.0;
   DECIMAL_TYPE wn;
   
   //Fs = interpFactor*Fs;
   
   hd = (DECIMAL_TYPE *) malloc(sizeof(DECIMAL_TYPE) * N);
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
   
   // apply window function
   for(n=0; n < N; n++)
      {
      //wn = 1;
      wn = 0.42 - 0.5 * cos((2 * M_PI * n)/(N-1)) + 0.08 * cos((4 * M_PI * n)/(N-1));
      //wn = 0.54-0.46*cos((2*M_PI*n)/(N-1));
      //wn = (1 - cos((2*M_PI*n) / (N-1)));
      h[n] = hd[n] * wn * (DECIMAL_TYPE)(interpFactor);
      }
   
   /*for(n=0;n < N;n++)
      {
      printf("h[%d]=%f\n",n,h[n]);
      }*/
   free(hd);
   return n;
   }