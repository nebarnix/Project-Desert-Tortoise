#include <math.h>
#include "MMClockRecovery.h"

//M&M Clock Recovery Loop (interpolating version!)
unsigned long MMClockRecovery(float *dataStreamIn, unsigned long numSamples, float *dataStreamOut, int Fs, float stepRange, float kp)
   {
   static char firstTime = 1;
   double baud=8320*2-1;
   static double stepSize;
   double stepMax = Fs/(baud-stepRange);
   double stepMin = Fs/(baud+stepRange);
   double currentBit;
   double Error;
   unsigned long count = 1;
   static double nextSample = 0;
   static double sampleLast = 1;
   if(firstTime == 1)
      {
      stepSize= Fs/(baud);
      firstTime = 0;
      }
   
      
   while(nextSample < numSamples)
      {    
      //Stores Bit 
      currentBit  = dataStreamIn[(unsigned int)(round(nextSample))];
      //dataStreamOutTime(count) = dataStreamInTime(round(nextSample));
      dataStreamOut[count] = currentBit;
      //Ind(count)  = nextSample;
      count = count + 1;
      
      //Calculates Error
      Error = sign(sampleLast)*currentBit - sign(currentBit)*sampleLast;
      
      //Updates Step Size
      stepSize = stepSize + kp*Error;
      
      //Limits Step size
      if( stepSize > stepMax )
         stepSize = stepMax;
      
      if( stepSize < stepMin )
         stepSize = stepMin;
      
      
      //Updates nextSample
      nextSample = nextSample + stepSize;
      sampleLast = currentBit;
      }
   nextSample = numSamples - nextSample;
   return count;
   }

int sign(float x) 
   {
   return (x > 0) - (x < 0);
   }