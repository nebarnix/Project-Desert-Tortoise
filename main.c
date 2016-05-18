#include <unistd.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
//#include <koolplot.h>
#include "wave.h"
#include "AGC.h"
#include "CarrierTrackPLL.h"
#include "LowPassFilter.h"
#include "MMClockRecovery.h"
#include "ManchesterDecode.h"
#include "ByteSync.h"

#define TRUE 1
#define FALSE 0


unsigned int CheckSum(float *dataStreamReal, unsigned long nSamples)
   {
   unsigned int sum=0;
   unsigned long idx;
   for(idx = 0; idx < nSamples; idx++)
      sum += dataStreamReal[idx];
   return sum;
   }

int main(int argc, char **argv) 
   {
   FILE *waveFilePtr;
   FILE *rawOutFilePtr;
   char *filename;
   HEADER header;
   double complex *waveData;
   float *dataStreamReal, *dataStreamSymbols;
   float PLLLock, percentComplete=0;
   float normFactor;
   unsigned long chunksize = 50000, nSamples, i=0, idx, nSymbols, nBits, totalSymbols=0, totalBits=0, totalSamples=0;
   unsigned char *dataStreamBits;
   float Fs;
   int LPF_Order, nFrames=0, totalFrames=0;
   float LPF_Fc;
   double *filterCoeffs;
   //unsigned int CheckSum1=0, CheckSum2=0, CheckSum3=0;
   
   LPF_Order = 26;
   filterCoeffs = malloc(sizeof(double) * LPF_Order);
   if (filterCoeffs == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }
   
   filename = (char*) malloc(sizeof(char) * 1024);
   if (filename == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }
      
   waveData = (double complex*) malloc(sizeof(double complex) * chunksize);
   if (waveData == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }
      
   dataStreamReal = (float*) malloc(sizeof(float) * chunksize);
   if (dataStreamReal == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }
      
   dataStreamSymbols = (float*) malloc(sizeof(float) * chunksize);
   if (dataStreamSymbols == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }
   
   dataStreamBits = (unsigned char*) malloc(sizeof(unsigned char) * chunksize);
   if (dataStreamBits == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }
    
      
   // get file path
   char cwd[1024];
   if (getcwd(cwd, sizeof(cwd)) != NULL) 
      {      
      strcpy(filename, cwd);
      
      // get filename from command line
      if (argc < 2)
         {
         printf("No wave file specified\n");
         return 1;
         }
   
      strcat(filename, "/");
      strcat(filename, argv[1]);
      printf("%s\n", filename);
      }

   // open file
   printf("Opening  file..\n");
   waveFilePtr = fopen(filename, "rb");
   if (waveFilePtr == NULL)
      {
      printf("Error opening file\n");
      exit(1);
      }
      
   rawOutFilePtr = fopen("output.raw", "wb");
   if (rawOutFilePtr == NULL)
      {
      printf("Error opening output file\n");
      exit(1);
      }
   //printf("Filepointer reference: %x\n",(unsigned int)waveFilePtr);
   //printf("Filepointer position: %x\n", ftell(waveFilePtr));
   header = ReadWavHeader(waveFilePtr);
   Fs = (float)header.sample_rate;
   printf("Sample Rate %.2fKHz\n", Fs/1000.0);
   //printHeaderInfo(header);
   //double complex test[] = {2,2,2,2};
   //printf("%d: Avg value: %f\n",i,StaticGain(test, 4, 1.0));
   
   long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);
   
   LPF_Fc = 11000;
   
   
   MakeLPFIR(filterCoeffs, LPF_Order, LPF_Fc, Fs);
   while(!feof(waveFilePtr))
      {
      nSamples = GetComplexWaveChunk(waveFilePtr, header, waveData, chunksize);
      if(i == 0)
         {
         normFactor = StaticGain(waveData, nSamples, 1.0);
         printf("Normalization Factor: %f\n",normFactor);
         }
         
      i+=nSamples;
      
      //normalize
      for(idx=0; idx < nSamples; idx++)
         {
         waveData[idx] *= normFactor;         
         }
      
      PLLLock = CarrierTrackPLL(waveData, dataStreamReal, nSamples, Fs, 4500, 0.1, 0.01, 0.001);      
      //CheckSum1 += CheckSum(dataStreamReal, nSamples);
      LowPassFilter(dataStreamReal, nSamples, filterCoeffs, LPF_Order);
      //CheckSum2 += CheckSum(dataStreamReal, nSamples);
      NormalizingAGC(dataStreamReal, nSamples, 0.00025);
      //CheckSum3 += CheckSum(dataStreamReal, nSamples);
      //nSymbols = MMClockRecovery(dataStreamReal, nSamples, dataStreamSymbols, Fs, 10, 0.15);
      nSymbols = MMClockRecovery(dataStreamReal, nSamples, dataStreamSymbols, Fs, 8, 0.15);
      
      //fwrite(dataStreamReal, sizeof(float), nSamples,rawOutFilePtr);
      //fwrite(dataStreamSymbols, sizeof(float), nSymbols,rawOutFilePtr);
      nBits = ManchesterDecode(dataStreamSymbols, nSymbols, dataStreamBits, 1.0);
      nFrames = FindSyncWords(dataStreamBits, nBits, "1110110111100010000", 19);      
      
      fwrite(dataStreamBits, sizeof(unsigned char), nBits, rawOutFilePtr);
      totalBits += nBits;
      totalFrames += nFrames;
      totalSymbols += nSymbols;
      totalSamples += nSamples;
      if(((float)( i) / num_samples)*100.0 - percentComplete > 0.15)
         {
         percentComplete = ((float)( i) / num_samples)*100.0;
         printf("\r%0.1f%% %ld Ks : %ld Symbols : %ld Bits : %d Frames", ((float)( i) / num_samples)*100.0,totalSamples/1000, totalSymbols, totalBits, totalFrames);      
         }
      
      //fwrite("2", sizeof(unsigned char), 1, rawOutFilePtr);
      /*for(idx=0; idx < nSamples; idx++)
         {
         fVal = (crealf(waveData[i]));
         fwrite(&fVal,sizeof(fVal),1,rawOutFilePtr);
         fVal = (cimagf(waveData[i]));
         fwrite(&fVal,sizeof(fVal),1,rawOutFilePtr);
         }*/
      }
   
   //printf("\nChecksum1=%X Checksum2=%X Checksum3=%X", CheckSum1,CheckSum2,CheckSum3);
   printf("\nAll done! Closing files and exiting.\nENJOY YOUR BITS AND HAVE A NICE DAY\n");
   fclose(rawOutFilePtr);
   fclose(waveFilePtr);
   
   // cleanup before quitting
   free(filename);
   free(dataStreamReal);
   free(dataStreamSymbols);
   free(filterCoeffs);
   free(dataStreamReal);
   free(waveData);
   free(dataStreamBits);
   
   //quit
   fflush(stdout);
   return 0;
   } 
   