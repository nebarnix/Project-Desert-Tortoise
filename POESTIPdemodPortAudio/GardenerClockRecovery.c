#include <math.h>
#include "GardenerClockRecovery.h"

//M&M Clock Recovery Loop (interpolating version!)
unsigned long GardenerClockRecovery(double *dataStreamIn, double *dataStreamInTime,  unsigned long numSamples, double *dataStreamOut, int Fs, double stepRange, double kp)
   {
   static char firstTime = 1;
   static double baud = 8320*2+03;
   static double stepSize;
   //double stepMax = Fs/(baud-stepRange);
   //double stepMin = Fs/(baud+stepRange);
   double currentBit;
   double Error;
   unsigned long count = 0;
   static double nextSample = 0;
   static double prevBit = 0;
   static double halfSample =0;
   
   if(firstTime == 1)
      {
      stepSize = Fs / baud;
      firstTime = 0;
      }
      
   while(rint(nextSample) < numSamples)
      {    
      //Stores Bit 
      currentBit  = dataStreamIn[(unsigned int)(rint(nextSample))];
      halfSample  = dataStreamIn[(unsigned int)(rint(halfSample))];
      //dataStreamOutTime(count) = dataStreamInTime(rint(nextSample));
      dataStreamOut[count] = currentBit;
      dataStreamInTime[count] = dataStreamInTime[(unsigned int)(rint(nextSample))];
      //Ind(count)  = nextSample;
      count = count + 1;
      
      //Calculates Error
      //Error = sign(prevBit)*currentBit - sign(currentBit)*prevBit;
      Error = kp * (currentBit - prevBit) * (halfSample);
      
      //Updates Step Size
      //stepSize = stepSize + kp*Error;
      
      //Limits Step size3
      if( Error > stepRange)
         Error = stepRange;
      else if( Error < -stepRange )
         Error = -stepRange;
           
      
      //Updates nextSample
      //nextSample = nextSample - Error;
      nextSample = (nextSample - Error);
      
      halfSample = nextSample + stepSize/2.0;
      nextSample = nextSample + stepSize;    
      prevBit = currentBit;
      }
   nextSample =  nextSample - numSamples; //roll over to next chunk
   //return count-1; //does this make things better?
   return count;
   }