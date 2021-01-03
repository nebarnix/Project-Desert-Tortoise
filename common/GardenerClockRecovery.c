#include <tgmath.h>
#include "GardenerClockRecovery.h"

//M&M Clock Recovery Loop (interpolating version!)
unsigned long GardenerClockRecovery(DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamInTime,  unsigned long numSamples, DECIMAL_TYPE *dataStreamOut, int Fs, DECIMAL_TYPE baud, DECIMAL_TYPE stepRange, DECIMAL_TYPE kp)
   {
   unsigned long count = 0;
   DECIMAL_TYPE currentBit;
   DECIMAL_TYPE Error;
   
   static char firstTime = 1;
   static DECIMAL_TYPE nextSample = 0;
   static DECIMAL_TYPE prevBit = 0;
   static DECIMAL_TYPE halfSample =0;   
   static DECIMAL_TYPE stepSize;
   
   if(firstTime == 1)
      {
      stepSize = Fs / baud;
      firstTime = 0;
      }
   
   #if USE_FLOATS==1
   while(rintf(nextSample) < numSamples)
      {    
      //Stores Bit 
      currentBit  = dataStreamIn[(unsigned int)(rintf(nextSample))];
      halfSample  = dataStreamIn[(unsigned int)(rintf(halfSample))];
      //dataStreamOutTime(count) = dataStreamInTime(rint(nextSample));
      dataStreamOut[count] = currentBit;
      dataStreamInTime[count] = dataStreamInTime[(unsigned int)(rintf(nextSample))];
      
      /*if(count > 0 && count >= (unsigned int)(rint(nextSample)))
         {
         printf("\nWhoops..\n");
         exit(1);
         }*/
      //Ind(count)  = nextSample;
      
      
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
      count ++;
      }
   
   dataStreamInTime[count] = dataStreamInTime[(unsigned int)(rintf(nextSample))];
   #else
   while(rint(nextSample) < numSamples)
      {    
      //Stores Bit 
      currentBit  = dataStreamIn[(unsigned int)(rint(nextSample))];
      halfSample  = dataStreamIn[(unsigned int)(rint(halfSample))];
      //dataStreamOutTime(count) = dataStreamInTime(rint(nextSample));
      dataStreamOut[count] = currentBit;
      dataStreamInTime[count] = dataStreamInTime[(unsigned int)(rint(nextSample))];
      
      /*if(count > 0 && count >= (unsigned int)(rint(nextSample)))
         {
         printf("\nWhoops..\n");
         exit(1);
         }*/
      //Ind(count)  = nextSample;
      
      
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
      count ++;
      }
   dataStreamInTime[count] = dataStreamInTime[(unsigned int)(rint(nextSample))];
   #endif   
   
   
   nextSample =  nextSample - numSamples; //roll over to next chunk
   //return count-1; //returning count-1 makes timing errors go away and makes sense because we end at count++, but doesn't fix time chunking problem
   return count; //returning count makes more frames pass parity, but leads to zeros in the time stream. I don't know why it makes more frames pass parity...
   }