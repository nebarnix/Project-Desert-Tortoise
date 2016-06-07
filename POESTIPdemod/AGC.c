#include <complex.h>
#include <math.h>
#include "AGC.h"

double FindSignalAmplitude(double *dataStreamIn, unsigned long nSamples, double alpha)
{
static double average = 0;
unsigned long i;

for(i=0; i < nSamples; i++)
   {    
    average = average * (1.0 - alpha) + alpha*fabs(dataStreamIn[i]);    
   }
return average;
}

double StaticGain(double complex *complexData,unsigned int nSamples,double desiredLevel)
   {
   ////n=100000; //number of points to average
   //if nargin < 2 || n==0 
       //n = numel(dataStreamIn);
   int i;
   double avgLevel=0;
   avgLevel = cabs(complexData[0]);
   for(i=1; i < nSamples; i++)
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