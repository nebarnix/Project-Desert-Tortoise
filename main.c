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

#define TRUE 1
#define FALSE 0

int main(int argc, char **argv) 
   {
   FILE *waveFilePtr;
   FILE *rawOutFilePtr;
   char *filename;
   HEADER header;
   double complex *waveData;
   float *dataStreamReal, *dataStreamSymbols;
   float PLLLock;
   float normFactor;
   //float fVal;
   unsigned long chunksize = 50000, nSamples, i=0, idx, nSymbols, nBits;
   unsigned char *dataStreamBits;
   
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
   printHeaderInfo(header);
   //double complex test[] = {2,2,2,2};
   //printf("%d: Avg value: %f\n",i,StaticGain(test, 4, 1.0));
   
   long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);
   
   while(!feof(waveFilePtr))
      {
      nSamples = GetComplexWaveChunk(waveFilePtr, header, waveData, chunksize);
      if(i == 0)
         {
         normFactor = StaticGain(waveData, nSamples, 1.0);
         printf("Normalization Factor: %f\n",normFactor);
         }
         
      i+=nSamples;
      printf("\r%0.1f%% complete", ((float)( i) / num_samples)*100.0);
      //normalize
      for(idx=0; idx < nSamples; idx++)
         {
         waveData[idx] *= normFactor;         
         }
      
      PLLLock = CarrierTrackPLL(waveData, dataStreamReal, nSamples, (float)header.sample_rate, 4500, 0.1, 0.01, 0.001);
      LowPassFilter(dataStreamReal, nSamples);
      NormalizingAGC(dataStreamReal, nSamples, 0.00025);      
      nSymbols = MMClockRecovery(dataStreamReal, nSamples, dataStreamSymbols, header.sample_rate, 10, 0.15);     
      
      //fwrite(dataStreamReal, sizeof(float), nSamples,rawOutFilePtr);
      //fwrite(dataStreamSymbols, sizeof(float), nSymbols,rawOutFilePtr);
      nBits = ManchesterDecode(dataStreamSymbols, nSymbols, dataStreamBits, 1.0);
      fwrite(dataStreamBits, sizeof(unsigned char), nBits, rawOutFilePtr);
      /*for(idx=0; idx < nSamples; idx++)
         {
         fVal = (crealf(waveData[i]));
         fwrite(&fVal,sizeof(fVal),1,rawOutFilePtr);
         fVal = (cimagf(waveData[i]));
         fwrite(&fVal,sizeof(fVal),1,rawOutFilePtr);
         }*/
      }
   
   
   printf("\nAll done! Closing files and exiting.\nENJOY YOUR BITS AND HAVE A NICE DAY\n");
   fclose(rawOutFilePtr);
   fclose(waveFilePtr);
   
   // cleanup before quitting
   free(filename);
   return 0;
   }