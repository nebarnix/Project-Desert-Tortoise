#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ByteSync.h"

/*void main()
   {
   char dataStreamBits[] = "111101101111000100000100000011101001100110010000000001000001000000001011010000111000000000000000001010101111110011101000010100110001000000011110000000000000000001011100011100000010011111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001110000100100000110110010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000001111111111001111011000100101110001010101001100001110110111100010000010000001110100110011001000010000100000100000010000010011000000000000010001110001010111111010111100001100000000100001001111000000000000000000111101111111111001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011000100001101010111000000011000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001100000000000000000001000010000101000000101000010000010000001000000101010100001000111011011110001000001000000111010011001100100010000010000010000011111110000000000011000010101110000001011001110000001000000010000010001000111100000000000000000011111110111111010101000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001101101110101111001011111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000010111111111110111111001010101101110000010101010011000111101101111000100000100000011101001100110010001100001000001000000111111100000000001010110100011100100101100111011100000101000001001000110011110000000000000000001111111011111101010100010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001110011110000000111011110110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011000000000000000011111111111111111111111111111111000000000000000001010101001100001110110111100010000010000001110100110011001001000000100000100000110111110000000000101100100110000000010101110101011010110110111000100100001111000000000000000000111111011111111001010001100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001100000000000000001111111111111001111111111111111111111111111111100101010100100000111011011110001000001000000111010011001100100101000010000010000001101010000000000010111100010011000101011101011011010000101001100010010100111100000000000000000011011100110110110101001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110000000000000000000001110000001000000110110010000000000000000100010101010000000011101101111000100000100000011101001100110010011000001000001000000100001100000000110000000000000011000101110100101111000011000000001001100011110011010110000000011101101011101000010100101000000011001101101101110000000000000000000011101000111100000000000000001001111110110110000000000000000000000000000000000000000000000000101100000110101000000000000000000011000011111001000100011000000100000010000001000001000000100111101111100001110100000000000000001010110001101001000000000000000011111011011010010000000000000000111110111001000100000000000000001000011101010111000000000000000011100100011110010000000000000000111100111011110000000000000000000000000000000000000000000000000000000000000000001010111010110100000000000000000010101000100000000000000000000010111001110000011000011000100001101011000000000011000000000000000001010101001000111110110111100010000010000001110100110011001001110000100000100000001111111010110010001111000100100101010110001100000010000000100000100111001111001000001011000100111111111111010001010011000000001100001100110011000000000000000001010111110000000000000000000000000000111111100100000000000000000000000000000000000000000000000011010110000001110000000000000000110010111011011100000000000100001111101001111100101010000111000000001110100110110000000000000000110100111101101000000000000000000110000010001011000000000000000000111010010011000000000000000000101010101011010100000000000000001101000100010111000000000000000011010010101111100000000000000000000000000000000000000000000000000000000000000000000101010111001000000000000000000101110010000101000000000000001111001000101110010001101101110101111101011111000000000000000000000101010100110110111011011110001000001000000111010011001100101000000010000010000011111100100001111011100100010010100001011000110111000001010000010010100000111100000001000101011111111100111110010101001110000000010001010111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011111111111111111111111111111111100000000000000000100000010011101000000000000000010010101010000000000000000000000110101100000001000000000000000001100101110010111000000000000000000001110101001110000000000000000101101101101110100000000000000000001000010111110000000000000000000000000000000000000000000000000000000000000000000110010111100110000000000000000001000101100100100000000000000100001001001001001111111111110001100000000000011001111111110100011010101010010001111101101111000100000100000011101001100110010100100001000001000000111010000101101000000000000000010000101000000000110101101101110001010010011110000101101001011001110011111110111010101000000000001100111110010000000000000000000110111011100010100000000000000000101000110001001000000000000000000000000000000000000000000000000000010111011010100000000000000000110111101110011111111111111100111111111111111110000000000000000110000010011101000000000000000001111001101111110000000000000000000110110110101010000000000000000000101001111101000000000000000000010110001110110000000000000000010100001010011010000000000000000110101100000001100000000000000000000000000000000000000000000000000000000000000001111001010110111000000000000000000001110101001110000000000000011000001111101011100000000000000110000000001110010111111110101000101010101000101111110110111100010000010000001110100110011001010100000100000100000000010011111111100000000000000000100010111101001110100001010011000101010001111001011000011010100111111100010011101010100100000001001001010001011000000000000000000011011010110110000000000000000100110000011000000000000000000000000000000000000000000000000000010101001001100110000000000000000011110011111101100000110111110110000011011001000000000000000000010000110010010010000000000000000100011100111000100000000000000000110100110110000000000000000000011001010001101100000000000000000110011101100110000000000000000001001010000010111000000000000000010001110010000010000000000000000000000000000000000000000000000000000000000000000110110011001001000000000000000001100001011110001000000000000001011001000101111100011000000001001111111111100001110111001001001000101010100101101111011011110001000001000000111010011001100101011000010000010000010101110100001000011010001000111100001010010100011110000110000000010101100111100101101011010101111111110111111010101010100000000100000101111100000000000000000001101011000000000000000000000000001010101101101110000000000000000000000000000000000000000000000000000111011010110000000000000000010011111110000110001100010000110101100100000001100000000000001000110000000010011000000000000000011111110010101110000000000000000011100110011111000000000000000001010010101100111000000000000000011011001000011010000000000000000100111111100001100000000000000001001111101010101000000000000000000000000000000000000000000000000000000000000000010100111001001010000000000000000000000000000000000000000000000110000000000000000000000000000110000000000000000000000000000000001010101010011001111101101111000100000100000011101001100110010110000001000001000000000010000000000001001101010111000000101111110010000100000001000001011000011110000000000000000001111111011111111010101011000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110110111011000000101111100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011000000000000000000000000000000010101101000101010010111000101000001010101000010011110110111100010000010000001110100110011001011010000100000100000000111010000000011000101010001100011010110011011110000010100000100101101001111000000000000000000111111111111111101010110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001100000000000000001111111111111111111111111111111100000000000000000101010100000000111011011110001000001000000111010011001100101110000010000010000010011000100001110000001010011000010101010010100001101011011011100010111000111100000000000000000011111111111111110101011010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000111111111111100011111111111111111111111111111111010101010001000111101101111000100000100000011101001100110010111100001000001000000111001110100011010100000001001101100101111110111101000010100110001011110011110000000000000000001111111111111011010101110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001101011000000110001100000000100100000000000000000000000000000000101001011111011100000000000000000000111011111000000000000000000001101001110110100000000000000000000100001101010000000000000000000100111100010010000000000000000001011111110000110000000000000000111101100111110000000000000000000000000000000000000000000000000000000000000000001101001110100100000000000000000011010001101001110000000000000011100010110100111100000110111101100000011011001000000000000000010001010101001010011110110111100010000010000001110100110011001100000000100000100000000011101000011000011100000000000110010111111001111100001100000000110000001111000010101011101001111111101111111101010111100000001101010111111010000000000000000010110100000100010000000000000000011011000111110000000000000000000000000000000000000000000000000001001010111000110000000000000000011110101011010100000000000000000000000000000000000000000000000010111010110000110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001100000000000000000001100010000110101101000000001100000000000000000101010100101001111011011110001000001000000111010011001100110001000010000010000001100000001011110000000000010010101101011111101000001000000010000011000100111100000000000000000011111101111111000101100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110000000000000000000110110111011000010101111100000000000000000000010101010010000111101101111000100000100000011101001100110011001000001000001000000100110011111111010101100001001001010101100111001100000101000001001100100011110000000000000000001111111111111111010110001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111111111111111111111111111111110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000010000001111111111111110100001011010110101010101001000001110110111100010000010000001110100110011001100110000100000100000000011100111110000000000000000000000010110011101011010110110111000110011001111000000000000000000111111111111111001011001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110101100000011111111111111110001111111111111111000000000000000001011011101101110000000000000000000011110011110000000000000000000011011011011110000000000000000011000000100010110000000000000000111101101011101000000000000000000100000011101011000000000000000011010010110100000000000000000000000000000000000000000000000000000000000000000000001010111010010000000000000000000011011011011110000000000000001100110110001101010101010010100100000010000010000000010011010011110101010100010101111011011110001000001000000111010011001100110100000010000010000001110011000000000000000000000000000001010111010111010000101001100011010000111100111101101110100010110110110111010101100110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011011110011000001101100100000000000000000000000000011010110000000000000000000000100110011000000000000000000101101110000111100000000000000000010011110110110000000000000000011101010001100000000000000000000101111100011100000000000000000000111001111011111000000000000000000000000000000000000000000000000000000000000000011111111110010100000000000000000100101110110110000000000000000111010111010001010000100110010010000010011000010100001001100111100010101010000101011101101111000100000100000011101001100110011010100001000001000001110111000000000000000000100011111010101110101101111000011000000001101010011110010111001000111011111010011111110010110100000000001011101000100100000000000000000010101101110101100000000000000000100100111001001000000000000000000000000000000000000000000000000110101111011111000000000000000000110001010011101000110001000011010110110000000110000000000000100010010100011010100000000000000000110110101110110000000000000000011000101100000100000000000000000011011110111111000000000000000001001111111010110000000000000000000000000000000010000000000000000110101110000111100000000000000000000000000000000000000000000000000000000000000000110001111001101000000000000000000101100010100000000000000000011011010100000011100110110110100100011010111000001000010101101000001010101000001111110110111100010000010000001110100110011001101100000100000100000111010111000001001001111101011111101010111010011000010000000100000110110001111000000010000000001111111101111110101011010100000001010101000100111000000000000000001110100110101100000000000000000000000101100110000000000000000000000000000000000000000000000000011010111000011110000000000000000001101011011011100011011011101100010010111110000000000000000000000110111000100000000000000000000010011000011001000000000000000001000110011011011000000000000000011101110111100110000000000000000110111110011101100000000000000001010000111010100000000000000000001000111101011100000000000000000000000000000000000000000000000000000000000000000000100110100010000000000000000001000101111010111000000000000001001100100111111110000110011001111000000000000000000000001001001000101010100101010";
   //char dataStreamBits[] = "11110110111100010000010000001110100110011001000000000100000100000000101101000011100000000000000000101010111111001110100001010011000100000001111000000000000000000101110001110000001001111100000000000000000000000000000000000000000000000000000000000000000000000000000";
   //char dataStreamBits[] = "11110110111100010000010000001110100110011001111011011110001000011101101111000100000000000001000001000000001011010000111000000000000000001010101111110011101000010100110001000000011110000000000000000001011100011100000010011111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001110000100100000110110010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000001111111111001111011000100101110001010101001100001110110111100010000010000001110100110011001000010000100000100000010000010011000000000000010001110001010111111010111100001100000000100001001111000000000000000000111101111111111001010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011000100001101010111000000011000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001100000000000000000001000010000101000000101000010000010000001000000101010100001000111011011110001000001000000111010011001100100010000010000010000011111110000000000011000010101110000001011001110000001000000010000010001000111100000000000000000011111110111111010101000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001101101110101111001011111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000010111111111110111111001010101101110000010101010011000111101101111000100000100000011101001100110010001100001000001000000111111100000000001010110100011100100101100111011100000101000001001000110011110000000000000000001111111011111101010100010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001110011110000000111011110110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
   //char dataStreamBits[] = "1110110111100010000";                      111011011110001000
   int numFrames = ByteSyncOnSyncword(dataStreamBits, strlen(dataStreamBits), "1110110111100010000", 19);
   printf("%d bits and %d Frames\n",strlen(dataStreamBits), numFrames);
   }*/

int ByteSyncOnSyncword(unsigned char *bitStreamIn, double *bitStreamInTime, unsigned long nSamples,  char *syncWord, unsigned int syncWordLength, FILE *minorFrameFile)
   {
   int idx2;
   static char firstTime = 1;
   static char *historyBufferCirc; //circular buffer
   static int oldest = 0, syncIndicator, frameByteIdx, minorFrameShiftFlag, bitIdx;
   int framesFound=0;
   long idx;
   unsigned char byte=0, zero=0, one=1;
   
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
      //Enter this loop if we found a syncword last time around            
      if(minorFrameShiftFlag == 1)
         {
         if(bitStreamIn[idx]=='0')     
            {
            byte = byte << 1; //This is a zero, just shift
            byte = byte | zero;
            }
         else
            {
            byte = byte << 1; //This is a one, set the bit then shift               
            byte = byte | one;
            }
         
         bitIdx++;   
         if(bitIdx > 7)
            {
            //minorFrame[frameByteIdx]=byte;
            fprintf(minorFrameFile,"%.2X ",byte);
            byte = 0;
            bitIdx = 0;
            frameByteIdx++;
            if(frameByteIdx > 103)
               {
               minorFrameShiftFlag = 0;
               fprintf(minorFrameFile,"\n");
               }
            }
         }         
         
      //overwrite oldest bit in cir buffer with newest bit   
      historyBufferCirc[oldest] = bitStreamIn[idx];       
      //printf("h[%d]=%c\n",oldest, historyBufferCirc[oldest]);
      syncIndicator = 1;
            
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
         fprintf(minorFrameFile,"%.5f ",bitStreamInTime[idx]);
         fprintf(minorFrameFile,"%.2X ",0b11101101);
         fprintf(minorFrameFile,"%.2X ",0b11100010);
         frameByteIdx = 2;
         minorFrameShiftFlag = 1;
         framesFound++;
         bitIdx=3;
         byte=0;
         zero = 0;
         one = 1;
         //bitIdx=;
         }       
      
         
      //Look for Inverse Syncword
      syncIndicator = 1;
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
         fprintf(minorFrameFile,"%.5f ",bitStreamInTime[idx]);
         fprintf(minorFrameFile,"%.2X ",0b11101101);
         fprintf(minorFrameFile,"%.2X ",0b11100010);
         frameByteIdx = 2;
         minorFrameShiftFlag = 1;
         framesFound++;
         bitIdx=3;
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
   