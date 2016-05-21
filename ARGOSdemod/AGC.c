#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "AGC.h"

void Squelch(double *dataStream, double *squelchStreamIn, unsigned long nSamples, double squelchThreshold)
{
unsigned long i;
static int lastSquelch=0;
int squelch;
for(i=0; i < nSamples; i++)
   {   
   if(squelchStreamIn[i] < squelchThreshold)
      {
      squelch = 1;
      dataStream[i] = 0;
      }
   else
      squelch = 0;
   if(lastSquelch > squelch)
      printf("Signal!\n");
   lastSquelch = squelch;
   }
}

double StaticGain(double complex *complexData,unsigned int nSamples,double desiredLevel)
   {
   ////n=100000; //number of points to average
   //if nargin < 2 || n==0 
       //n = numel(dataStreamIn);
   unsigned long i;
   double avgLevel=0;
   avgLevel = cabs(complexData[0]);
   for(i=0; i < nSamples; i++)
      {
      avgLevel += cabs(complexData[i]);
      avgLevel /= 2.0;
      }
   //return avgLevel;
      
   return desiredLevel / avgLevel  ;
   }
   
//Automatic Gain Control Block
void NormalizingAGC(double *dataStreamIn, unsigned long nSamples, double AGC_loop_gain)
   {
   //Todo: implement a 'relock' mode for LARGE error values (either adjust gain
   //outright or adjust loop gain
   
   unsigned long idx;
   
   static double gain = 1; //Initial Gain Value
   double desired = 0.6366; //because a sin wave of amplitude 1 has this average absolute value
   double error;
   
   for(idx=0; idx< nSamples; idx++)
      {
      dataStreamIn[idx] *= gain;
      error = desired - (gain * fabs(dataStreamIn[idx]));
      
      //if(error > 1.0) //do something for really large errors
      //    error = -1/error;
      //end
      
      gain = gain + AGC_loop_gain * error ;
      //gaini(idx) = gain;
      }
   }
   
   //Automatic Gain Control Block
void NormalizingAGCC(double complex *dataStreamIn, unsigned long nSamples, double initial, double AGC_loop_gain)
   {
   //Todo: implement a 'relock' mode for LARGE error values (either adjust gain
   //outright or adjust loop gain
   
   unsigned long idx;
   
   static double gain = 0; //Initial Gain Value
   double desired = 5; //because a sin wave of amplitude 1 has this average absolute value
   double error;
   
   if(gain == 0)
      {
      gain = initial;
      }
   
   for(idx=0; idx< nSamples; idx++)
      {
      dataStreamIn[idx] *= gain;
      error = desired - (gain * fabs(dataStreamIn[idx]));
      
      //if(error > 1.0) //do something for really large errors
      //    error = -1/error;
      //end
      
      gain = gain + AGC_loop_gain * error ;
      //gaini(idx) = gain;
      }
   }