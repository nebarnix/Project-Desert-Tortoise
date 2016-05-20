#include <math.h>
#include "MMClockRecovery.h"

//M&M Clock Recovery Loop (interpolating version!)
unsigned long MMClockRecovery(float *dataStreamIn, unsigned long numSamples, float *dataStreamOut, int Fs, float stepRange, float kp)
   {
   static char firstTime = 1;
   static double baud = 8320*2 - 1;
   static double stepSize=3;
   double stepMax = Fs/(baud-stepRange);
   double stepMin = Fs/(baud+stepRange);
   double currentBit;
   double Error;
   unsigned long count = 0;
   static double nextSample = 0;
   static double sampleLast = 0;
   
   if(firstTime == 1)
      {
      stepSize = Fs/(baud);
      firstTime = 0;
      }
      
   while(round(nextSample) < numSamples)
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
      
      //Limits Step size3
      if( stepSize > stepMax )
         stepSize = stepMax;
      
      if( stepSize < stepMin )
         stepSize = stepMin;
      
      
      //Updates nextSample
      nextSample = nextSample + stepSize;
      sampleLast = currentBit;
      }
   nextSample =  nextSample - numSamples; //roll over to next chunk
   //return count-1; //does this make things better?
   return count;
   }

int sign(float x) 
   {
   return (x > 0) - (x < 0);
   }
