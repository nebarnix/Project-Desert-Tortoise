#include <math.h>
#include <stdio.h>
#include "ManchesterDecode.h"
#include "MMClockRecovery.h"

//Alternate between sample and evaluate phase
//First bit time of the manchester pair is used to evaluate phase
//Second bit time of the manchester pair is used to sample the bit (using both)

unsigned long ManchesterDecode(double *dataStreamIn, double *dataStreamInTime, unsigned long nSymbols, unsigned char *bitStream, double resyncThreshold)
   {
   //convert to bits from raw manchester bits
   //resyncThreshold = 1;
   
   unsigned long idxi, idxo=0;
   static unsigned int clockmod = 0;
   //static unsigned long idxerr=1;
   static double currentSample=0;
   static double prevSample=0;
   static double prevPrevSample=0;
   static unsigned char evenOddCounter=0;
   
   unsigned char currentBit;
   
   
   //for idxi = 2:numel(dataStreamIn)-1    
   for(idxi = 0; idxi < nSymbols; idxi++, evenOddCounter++)
      {      
      prevPrevSample = prevSample;
      prevSample = currentSample;
      currentSample = dataStreamIn[idxi];      
      
      //If not a bit boundary, see if it should be and we're out of sync
      //But only resync on strong bits
      if((evenOddCounter % 2) != clockmod)
         {
         if(sign(prevPrevSample) == sign(prevSample))
            {    
            if(fabs(prevPrevSample) > resyncThreshold && fabs(prevSample) > resyncThreshold)
               {
               //printf("\nResync: %f=%f %ld %ld\n",prevSample,prevPrevSample,idxi,idxo);
               clockmod = (evenOddCounter % 2); //only resync if we have confidence in BOTH bit decisions
               }
            }        
         }
       
      //check for bit boundary, and make decision using the strongest of the
      //two bits. 
      if((evenOddCounter % 2) == clockmod)
         {
         if(fabs((prevSample)) > fabs((currentSample))) //use the strongest symbol to determine bit
            {
            if(prevSample > 0)                
            currentBit = '1';
            else
            currentBit = '0';
            }            
         else
            {
            if(currentSample > 0)                
            currentBit = '0';
            else
            currentBit = '1';                    
            }        
         //bitTime(idxo)=rawTime(idxi);  
         bitStream[idxo] = currentBit;
         dataStreamInTime[idxo] = dataStreamInTime[idxi];
         idxo++;         
         }      
      
      //printf(" I:%ldO:%ld:%ld ",idxi,idxo,evenOddCounter);   
      
      }
   //fprintf([num2str(idxerr) ' errors\n']);
   //printf("\n%ld/2=%ld == %ld\n",idxi,idxi/2,idxo);
   return idxo;
   }