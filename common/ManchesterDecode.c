#include <tgmath.h>
#include <stdio.h>
#include "ManchesterDecode.h"
#include "MMClockRecovery.h"

//Alternate between sample and evaluate phase
//First bit time of the manchester pair is used to evaluate phase
//Second bit time of the manchester pair is used to sample the bit (using both)

unsigned long ManchesterDecode(DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamTime, unsigned long nSymbols, unsigned char *bitStream, DECIMAL_TYPE resyncThreshold)
   {
   //convert to bits from raw manchester bits
   //resyncThreshold = 1;
   
   unsigned long idxi, idxo=0;
   static unsigned int clockmod = 0;
   //static unsigned long idxerr=1;
   static DECIMAL_TYPE currentSample=0;
   static DECIMAL_TYPE prevSample=0;
   static DECIMAL_TYPE prevPrevSample=0;
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
            #if USE_FLOATS==1
               if(fabsf(prevPrevSample) > resyncThreshold && fabsf(prevSample) > resyncThreshold)
                  {
                  //printf("\nResync: %f=%f %ld %ld\n",prevSample,prevPrevSample,idxi,idxo);
                  clockmod = (evenOddCounter % 2); //only resync if we have confidence in BOTH bit decisions
                  }
            #else
               if(fabs(prevPrevSample) > resyncThreshold && fabs(prevSample) > resyncThreshold)
                  {
                  //printf("\nResync: %f=%f %ld %ld\n",prevSample,prevPrevSample,idxi,idxo);
                  clockmod = (evenOddCounter % 2); //only resync if we have confidence in BOTH bit decisions
                  }
            #endif            
            }        
         }
       
      //check for bit boundary, and make decision using the strongest of the
      //two bits. 
      if((evenOddCounter % 2) == clockmod)
         {
         #if USE_FLOATS==1
            if(fabsf((prevSample)) > fabsf((currentSample))) //use the strongest symbol to determine bit
               {
               if(prevSample > 0)                
                  currentBit = '1';
               else
                  currentBit = '0';
               }
         #else
            if(fabs((prevSample)) > fabs((currentSample))) //use the strongest symbol to determine bit
               {
               if(prevSample > 0)                
                  currentBit = '1';
               else
                  currentBit = '0';
               }
         #endif
                     
         else
            {
            if(currentSample > 0)                
               currentBit = '0';
            else
               currentBit = '1';                    
            }        
         //bitTime(idxo)=rawTime(idxi);  
         bitStream[idxo] = currentBit;
         dataStreamTime[idxo] = dataStreamTime[idxi];
         /*if(idxo > 0 && idxo >= idxi)
            {
            printf("\nWhoops! %d >= %d\n", idxo, idxi);               
            exit(1);
            }*/
         idxo++;         
         }      
      
      //printf(" I:%ldO:%ld:%ld ",idxi,idxo,evenOddCounter);   
      
      }
   //fprintf([num2str(idxerr) ' errors\n']);
   //printf("\n%ld/2=%ld == %ld\n",idxi,idxi/2,idxo);
   return idxo;
   }