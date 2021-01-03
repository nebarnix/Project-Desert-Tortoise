#include <tgmath.h>
#include "MMClockRecovery.h"

//M&M Clock Recovery Loop (interpolating version!)
unsigned long MMClockRecovery(DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamInTime,  unsigned long numSamples, DECIMAL_TYPE *dataStreamOut, int Fs, DECIMAL_TYPE baud, DECIMAL_TYPE stepRange, DECIMAL_TYPE kp)
   {
   static char firstTime = 1;
   //static DECIMAL_TYPE baud = 8320*2 - 1; //this is now a parameter
   static DECIMAL_TYPE stepSize=3.0;
   DECIMAL_TYPE stepMax = Fs/(baud-stepRange);
   DECIMAL_TYPE stepMin = Fs/(baud+stepRange);
   DECIMAL_TYPE currentBit;
   DECIMAL_TYPE Error;
   unsigned long count = 0;
   static DECIMAL_TYPE nextSample = 0;
   static DECIMAL_TYPE sampleLast = 0;
   
   if(firstTime == 1)
      {
      stepSize = Fs/(baud);
      firstTime = 0;
      }
   
   #if USE_FLOATS==1   
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
   #else
   while(rintf(nextSample) < numSamples)
         {    
         //Stores Bit 
         currentBit  = dataStreamIn[(unsigned int)(rintf(nextSample))];
         //dataStreamOutTime(count) = dataStreamInTime(rint(nextSample));
         dataStreamOut[count] = currentBit;
         dataStreamInTime[count] = dataStreamInTime[(unsigned int)(rintf(nextSample))];
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
   #endif
   nextSample =  nextSample - numSamples; //roll over to next chunk
   //return count-1; //does this make things better?
   return count;
   }

int sign(DECIMAL_TYPE x) 
   {
   return (x > 0) - (x < 0);
   }