#include <math.h>
#include "MMClockRecovery.h"

//M&M Clock Recovery Loop (interpolating version!)
unsigned long MMClockRecovery(double *dataStreamIn, double *dataStreamInTime,  unsigned long numSamples, double *dataStreamOut, int Fs, double stepRange, double kp)
   {
   static char firstTime = 1;
   static double baud = 8320*2 - 1;
   static double stepSize=3.0;
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
      
   while(rint(nextSample) < numSamples)
      {    
      //Stores Bit 
      currentBit  = dataStreamIn[(unsigned int)(rint(nextSample))];
      //dataStreamOutTime(count) = dataStreamInTime(rint(nextSample));
      dataStreamOut[count] = currentBit;
      dataStreamInTime[count] = dataStreamInTime[(unsigned int)(rint(nextSample))];
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

int sign(double x) 
   {
   return (x > 0) - (x < 0);
   }