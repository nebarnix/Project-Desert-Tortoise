#include <complex.h>
#include <tgmath.h>
#include <stdio.h>
#include "AGC.h"

DECIMAL_TYPE FindSignalAmplitude(DECIMAL_TYPE *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE alpha)
{
static DECIMAL_TYPE average = 0;
unsigned long i;

for(i=0; i < nSamples; i++)
   {    
   #if USE_FLOATS==1
      average = average * (1.0 - alpha) + alpha*fabsf(dataStreamIn[i]);
   #else
      average = average * (1.0 - alpha) + alpha*fabs(dataStreamIn[i]);
   #endif
   }
return average;
}

//squelch uses a locksignal or perhaps RSSI to determine a cutoff point. If the stream falls below the threshold value it is zero stuffed
//Don't squelch before an AGC for this reason!
void Squelch(DECIMAL_TYPE *dataStream, DECIMAL_TYPE *squelchStreamIn, unsigned long nSamples, DECIMAL_TYPE squelchThreshold)
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
   /*
   if(lastSquelch > squelch)
      printf("Signal!\r");
   else if(lastSquelch < squelch)
      printf("       \r");
   */
   lastSquelch = squelch;
   }
}

DECIMAL_TYPE StaticGain(DECIMAL_TYPE complex *complexData, unsigned int nSamples, DECIMAL_TYPE desiredLevel)
   {
   ////n=100000; //number of points to average
   //if nargin < 2 || n==0 
       //n = numel(dataStreamIn);
   unsigned long i;
   DECIMAL_TYPE avgLevel=0;
   
   #if USE_FLOATS==1
      avgLevel = cabsf(complexData[0]);
   #else
      avgLevel = cabs(complexData[0]);
   #endif
   
   for(i=0; i < nSamples; i++)
      {
      #if USE_FLOATS==1
         avgLevel += cabsf(complexData[i]);
      #else
         avgLevel += cabs(complexData[i]);
      #endif
      
      avgLevel /= 2.0;
      }
   //return avgLevel;
      
   return desiredLevel / avgLevel  ;
   }

//Automatic Gain Control Block (GNUradio based )
void NormalizingAGC(DECIMAL_TYPE *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE initial, DECIMAL_TYPE attack_rate, DECIMAL_TYPE decay_rate)
   {
   
   unsigned long idx;
   
   static DECIMAL_TYPE gain = 1; //Initial Gain Value
   static char firsttime = 1;
   DECIMAL_TYPE rate;
   //double attack_rate = 1e-1; 
   //double decay_rate = 1e-1;
   DECIMAL_TYPE reference = 1.0;  
   DECIMAL_TYPE max_gain = 5000;   
   DECIMAL_TYPE error;
   
   if(firsttime == 1)
      {
      firsttime = 0;
      gain = initial;
      }
   
   for(idx=0; idx< nSamples; idx++)
      {     
   
      //output = input * _gain;
      dataStreamIn[idx] *= gain;
      
      #if USE_FLOATS==1
         error = (fabsf(dataStreamIn[idx])) - reference;
      #else
         error = (fabs(dataStreamIn[idx])) - reference;
      #endif
      
      rate = decay_rate;
      
      #if USE_FLOATS==1
         if(fabsf(error) > gain)
      #else
         if(fabs(error) > gain)
      #endif       
         {
         rate = attack_rate;
         }
         
      gain -= error*rate;
      
      // Not sure about this
      if(gain < 0.0)
         gain = 10e-5;
      
      if(max_gain > 0.0 && gain > max_gain)
         {
         gain = max_gain;
         }     
      }
   }
   
   
   /*
//Automatic Gain Control Block (nebarnix original)
void NormalizingAGC(double *dataStreamIn, unsigned long nSamples, double AGC_loop_gain)
   {
   //Todo: implement a 'relock' mode for LARGE error values (either adjust gain
   //outright or adjust loop gain
   
   unsigned long idx;
   
   static double gain = 1; //Initial Gain Value
   double desired = 10;//6.6366;//0.6366; //because a sin wave of amplitude 1 has this average absolute value
   
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
   }*/
   
   //Automatic Gain Control Block
void NormalizingAGCC(DECIMAL_TYPE complex *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE initial, DECIMAL_TYPE AGC_loop_gain)
   {
   //Todo: implement a 'relock' mode for LARGE error values (either adjust gain
   //outright or adjust loop gain
   
   unsigned long idx;
   
   static DECIMAL_TYPE gain = 0.0; //Initial Gain Value
   DECIMAL_TYPE desired = 5; 
   DECIMAL_TYPE error;
   static char firsttime=1;
   
   if(firsttime == 1)
      {
      firsttime = 0;
      gain = initial;
      }
   
   for(idx=0; idx< nSamples; idx++)
      {
      dataStreamIn[idx] *= gain;
      
      #if USE_FLOATS==1
         error = desired - (gain * fabsf(dataStreamIn[idx]));
      #else
         error = desired - (gain * fabs(dataStreamIn[idx]));
      #endif
      
      
      //if(error > 1.0) //do something for really large errors
      //    error = -1/error;
      //end
      
      gain = gain + AGC_loop_gain * error ;
      //gaini(idx) = gain;
      }
   }