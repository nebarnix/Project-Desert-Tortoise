#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ByteSync.h"

#define DOFILEIO
/*
void main()
   {
   //FILE *rawOutFilePtr2;
   //rawOutFilePtr2 = fopen("nul", "wb");
   //char dataStreamBits[] = "";
   char dataStreamBits[] = "1111111111111111111111110110111111110110110111011010001000100010101111001010101111101001000001001111111111111110001011110";
   char dataStreamBits2[] = "000000000011101100000010000000001010101101010101111111111011111111111111111111111111111111111111111111111111";
   //char dataStreamBits[] = "11110110111100010000010000001110100110011001000000000100000100000000101101000011100000000000000000101010111111001110100001010011000100000001111000000000000000000101110001110000001001111100000000000000000000000000000000000000000000000000000000000000000000000000000";
   //char dataStreamBits[] = "11110110111100010000010000001110100110011001111011011110001000011101101111000100000000000001000001000000001011010000111000000000000000001010101111110011101000010100110001000000011110000000000000000001011100011100000010011111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001110000100100000110110010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000001111111111001111011000100101110001010101001100001110110111100010000010000001110100110011001000010000100000100000010000010011000000000000010001110001010111111010111100001100000000100001001111000000000000000000111101111111111001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011000100001101010111000000011000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001100000000000000000001000010000101000000101000010000010000001000000101010100001000111011011110001000001000000111010011001100100010000010000010000011111110000000000011000010101110000001011001110000001000000010000010001000111100000000000000000011111110111111010101000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001101101110101111001011111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000010111111111110111111001010101101110000010101010011000111101101111000100000100000011101001100110010001100001000001000000111111100000000001010110100011100100101100111011100000101000001001000110011110000000000000000001111111011111101010100010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001110011110000000111011110110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
   //char dataStreamBits[] = "1110110111100010000";                      111011011110001000
   int numFrames = FindSyncWords(dataStreamBits, NULL, strlen(dataStreamBits), "0001011110000", 13, NULL);
   numFrames  += FindSyncWords(dataStreamBits2, NULL, strlen(dataStreamBits2), "0001011110000", 13, NULL);
   printf("%d bits and %d Frames\n",strlen(dataStreamBits)+strlen(dataStreamBits2), numFrames);
   //fclose(rawOutFilePtr2);
   }
*/
int FindSyncWords(unsigned char *bitStreamIn, double *bitStreamInTime, unsigned long nSamples,  char *syncWord, unsigned int syncWordLength, FILE *minorFrameFile)
   {
   int idx2;
   static char firstTime = 1;
   static char *historyBufferCirc; //circular buffer
   static int oldest = 0, syncIndicator, frameByteIdx, minorFrameShiftFlag=0, bitIdx;
   int framesFound=0;
   long idx;
   unsigned char zero=0, one=1;
   static unsigned char byte =0;
   
   if(firstTime == 1)
      {
      firstTime = 0;
      //printf("Asking for %d Bytes",sizeof(char) * syncWordLength);
      historyBufferCirc = (char *) malloc(sizeof(char) * syncWordLength);
      if (historyBufferCirc == NULL)
            {
            printf("Error in malloc\n");
            exit(1);
            }
      memset(historyBufferCirc, '.', sizeof(char) * syncWordLength);
      }
   
   //loop through samples
   for(idx = 0; idx < nSamples; idx++)
      {
      //overwrite oldest bit in cir buffer with newest bit   
      historyBufferCirc[oldest] = bitStreamIn[idx];       
      //printf("h[%d]=%c\n",oldest, historyBufferCirc[oldest]);
      syncIndicator = 1;
      
       //since we've already advanced to bit 19 of the frame....
      //for(frameByteIdx=0; frameByteIdx < 103; frameByteIdx++) //minor frames are 103 bytes long
      //11101101 11100010 000
      if(minorFrameShiftFlag == 1)
         {
         if(bitStreamIn[idx]=='0')     
            {
            byte = byte << 1; //This is a zero, just shift
            byte = byte | zero;
            //printf("0 %.2X\n", byte);
            }
         else
            {
            byte = byte << 1; //This is a one, set the bit then shift               
            byte = byte | one;
            //printf("1 %.2X\n", byte);
            }
         
         bitIdx++;   
         if(bitIdx > 7)
            {
            //minorFrame[frameByteIdx]=byte;
            #ifdef DOFILEIO
               fprintf(minorFrameFile,"%.2X ",byte);
            #endif
            
            printf("%.2X ",byte);
            byte = 0;
            bitIdx = 0;
            frameByteIdx++;
            if(frameByteIdx > 8)
               {
               minorFrameShiftFlag = 0;
               #ifdef DOFILEIO
               fprintf(minorFrameFile,"\n");
               #endif
               printf("\n");
               }
            }
         }  
      
      //Look for syncword
      for (idx2 = 0; idx2 < syncWordLength; idx2++) 
         {
         //compare syncword bytes to appropriate circular buffer bytes
         
         if (syncWord[idx2] != historyBufferCirc[(oldest + idx2 + 1) % syncWordLength])
            {
            syncIndicator = 0;
            break;
            }
         }
         
      if(syncIndicator == 1 && minorFrameShiftFlag == 0)
         {
         //gotoxy(1,1);
         #ifdef DOFILEIO
            fprintf(minorFrameFile,"%.5f ",bitStreamInTime[idx]);
            fprintf(minorFrameFile,"%.2X ",0b11100010);
            fprintf(minorFrameFile,"%.2X ",0b11110000);
         #endif
         
         printf("%.5f ",bitStreamInTime[idx]);
         printf("%.2X ",0b11100010);
         printf("%.2X ",0b11110000);
         
         frameByteIdx = 2;
         minorFrameShiftFlag = 1;
         framesFound++;
         bitIdx=0;
         byte=0;
         zero = 0;
         one = 1;
         //bitIdx=;
         }       
      
         
      //Look for Inverse Syncword
      syncIndicator = 0;
      for (idx2 = 0; idx2 < syncWordLength; idx2++) 
         {
         //compare syncword bytes to appropriate circular buffer bytes
         //printf("sw[%d] == hst[%d] , %c == %c, ", idx2,(oldest + idx2 + 1) % syncWordLength,syncWord[idx2],historyBufferCirc[(oldest + idx2) % syncWordLength]);
         
         if (syncWord[idx2] == historyBufferCirc[(oldest + idx2 + 1) % syncWordLength])
            {
            syncIndicator = 0;
            //printf("NO\n\n");
            break;
            }         
         }
         
      if(syncIndicator == 1 && minorFrameShiftFlag == 0)
         {
         //gotoxy(1,1);
         //printf("Backwards?\n");
         #ifdef DOFILEIO
            fprintf(minorFrameFile,"%.5f i ",bitStreamInTime[idx]);
            fprintf(minorFrameFile,"%.2X ",0b11100010);
            fprintf(minorFrameFile,"%.2X ",0b11110000);         
         #endif
         
         printf("%.5f i ",bitStreamInTime[idx]);         
         printf("%.2X ",0b11100010);
         printf("%.2X ",0b11110000);
         
         frameByteIdx = 2;
         minorFrameShiftFlag = 1;
         framesFound++;
         bitIdx=0;
         byte=0;
         zero = 1;
         one = 0;
         //bitIdx=;
         }
     
      //advance oldest bit pointer   
      oldest = (oldest + 1) % syncWordLength;
      }
   return framesFound;   
   }
/*{
                                      1110110111100010000
   static unsigned char syncWord[] = "1110110111100010000";
   //static unsigned char syncWordInverse[] = "0001001000011101111";

   unsigned long idx;
   static int syncIndicator = 0;
   int numSyncWords=0;
   
   for(idx = 0; idx < nSamples; idx++)
      {
      if(bitStreamIn[idx] == syncWord[syncIndicator])
         {
         syncIndicator++;
         if(syncIndicator == 18)
            {
            syncIndicator = 0;
            numSyncWords++;
            }
         }
      else 
         syncIndicator = 0;      
      }
   //printf("\r %d Syncwords",  numSyncWords);    
   return numSyncWords;
   }*/
   //syncWordAllIndex = sort(cat(2,syncWordIndex,syncWordInvIndex));
   
   /*
   for frameIdx=1:numel(syncWordAllIndex)-1
       //See if the frame is normal or inverted bits
       if isempty(find(syncWordInvIndex == syncWordAllIndex(frameIdx),1))
           for frameByteIdx=0:103 //minor frames are 103 bytes long
               byte=0;
               //Start of byte time
               frameTime(frameIdx,frameByteIdx+1)=bitTime(syncWordAllIndex(frameIdx)+frameByteIdx*8);
               //if this is a normal sync word, use normal bits
               for bitIdx=0:7  //bytes are 8 bits long ;)                            
                   if(dataStreamIn(syncWordAllIndex(frameIdx)+frameByteIdx*8+bitIdx)=='0')               
                       byte = bitshift(byte,1); //This is a zero, just shift           
                   else                    
                       byte = bitshift(byte,1); //This is a one, set the bit then shift               
                       byte = bitor(byte,1);              
                   end                
               end        
           minorFrames(frameIdx,frameByteIdx+1)=byte;    
           end
       else //this minor frame is inverted
           for frameByteIdx=0:103 //minor frames are 103 bytes long
               byte=0;
               //Start of byte time
               frameTime(frameIdx,frameByteIdx+1)=bitTime(syncWordAllIndex(frameIdx)+frameByteIdx*8);            
               for bitIdx=0:7  //bytes are 8 bits long ;)                            
                   if(dataStreamIn(syncWordAllIndex(frameIdx)+frameByteIdx*8+bitIdx)=='0')               
                       byte = bitshift(byte,1); //This is a zero, just shift
                       byte = bitor(byte,1);              
                   else                    
                       byte = bitshift(byte,1); //This is a one, set the bit then shift                                  
                   }                
               }        
           minorFrames(frameIdx,frameByteIdx+1)=byte;    
           }
       }
   */
   