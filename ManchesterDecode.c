#include <math.h>
#include "ManchesterDecode.h"
#include "MMClockRecovery.h"

unsigned long ManchesterDecode(float *dataStreamIn, unsigned long nSymbols, unsigned char *bitStream, float resyncThreshold)
   {
   //convert to bits from raw manchester bits
   //resyncThreshold = 1;
   
   unsigned long idx, idx2=0;
   static int clockmod = 0;
   static unsigned long idxerr=1;
   static float currentSample=0;
   static float prevSample=0;
   static float prevPrevSample=0;
   static char evenOddCounter=0;
   
   unsigned char currentBit;
   
   
   //for idx = 2:numel(dataStreamIn)-1    
   for(idx = 0; idx < nSymbols; idx++)
      {      
      prevPrevSample = prevSample;
      prevSample = currentSample;
      currentSample = dataStreamIn[idx];      
      
      evenOddCounter++;
      //If not a bit boundary, see if it should be and we're out of sync
      //But only resync on strong bits
      if((evenOddCounter % 2) != clockmod)
         {
         if(sign(prevPrevSample) == sign(prevSample))
            {
            //errx[idxerr]=idx2;
            //erry[idxerr]=(dataStreamIn[idx]);
            idxerr=idxerr+1;   
            if(fabs(prevPrevSample) > resyncThreshold && fabs(prevSample) > resyncThreshold)                
               clockmod = (evenOddCounter % 2); //only resync if we have confidence in BOTH bit decisions                        
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
           //bitTime(idx2)=rawTime(idx);  
         bitStream[idx2] = currentBit;
         idx2++;
         }      
      }
   //fprintf([num2str(idxerr) ' errors\n']);
   return idx2;
   }