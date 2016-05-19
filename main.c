#include <ctype.h>
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

//int spaceCraftID, dayNum, 

unsigned int CheckSum(unsigned char *dataStreamReal, unsigned long nSamples)
   {
   unsigned int sum=0;
   unsigned long idx;
   for(idx = 0; idx < nSamples; idx++)
      {
      sum += (dataStreamReal[idx]);
      //printf("%.2X %ld,",dataStreamReal[idx]);
      }      
   //printf("\n");
   return sum;
   }

int main(int argc, char **argv) 
   {
   FILE *waveFilePtr=NULL;
   FILE *rawOutFilePtr=NULL;
   FILE *minorFrameFile=NULL;
   
   HEADER header;
   
   unsigned long chunksize = 100000, nSamples, i=0, idx, nSymbols, nBits, totalSymbols=0, totalBits=0, totalSamples=0;
   
   float *dataStreamReal=NULL, *dataStreamSymbols=NULL;
   float Fs;
   float LPF_Fc;   
   float PLLLock, percentComplete=0;
   float normFactor=0;
   
   double *filterCoeffs=NULL;
   double complex *waveData=NULL;
   
   unsigned int CheckSum1=0, CheckSum2=0, CheckSum3=0;
   int LPF_Order, nFrames=0, totalFrames=0,c;
   
   unsigned char *dataStreamBits=NULL;
   char *filename=NULL;
   char *cvalue = NULL;   
   
   while ((c = getopt (argc, argv, "nc:")) != -1)
      {
      switch (c)
         {
         case 'n':
            cvalue = optarg;
            normFactor = atoi(cvalue);
            printf("Static Gain Override %f\n",normFactor);
            break;
         case 'c':
            cvalue = optarg;
            chunksize = atoi(cvalue);
            printf("Using %ld chunksize\n",chunksize);
            break;
         case '?':
            if (optopt == 'c' || optopt == 'n')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
               fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
               fprintf (stderr,
               "Unknown option character `\\x%x'.\n",
               optopt);
               return 1;
         default:
            abort ();
         }
      }
   
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
      strcat(filename, argv[optind]);
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
      
   
   minorFrameFile = fopen("minorFrame.txt", "w");
   if (minorFrameFile == NULL)
      {
      printf("Error opening output file\n");
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
   long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);
   
   printf("Sample Rate %.2fKHz. Total samples %ld\n", Fs/1000.0,num_samples);
   //printHeaderInfo(header);
   //double complex test[] = {2,2,2,2};
   //printf("%d: Avg value: %f\n",i,StaticGain(test, 4, 1.0));
   
   LPF_Fc = 11000;   
   MakeLPFIR(filterCoeffs, LPF_Order, LPF_Fc, Fs);
   while(!feof(waveFilePtr))
      { //num_samples
      nSamples = GetComplexWaveChunk(waveFilePtr, header, waveData, chunksize);
      
      if(i == 0 && normFactor == 0)
         {
         normFactor = StaticGain(waveData, nSamples, 1.0);
         //normFactor = 0.000041;
         printf("Normalization Factor: %f\n",normFactor);
         }
         
      i+=nSamples;
      
      //normalize
      for(idx=0; idx < nSamples; idx++)
         {
         //if(feof(waveFilePtr))
            
         waveData[idx] *= normFactor;         
         }
      CheckSum1 += CheckSum((unsigned char*)waveData, sizeof(*waveData) * nSamples);
      PLLLock = CarrierTrackPLL(waveData, dataStreamReal, nSamples, Fs, 4500, 0.1, 0.01, 0.001);      
      //CheckSum1 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      LowPassFilter(dataStreamReal, nSamples, filterCoeffs, LPF_Order);
      //CheckSum2 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      NormalizingAGC(dataStreamReal, nSamples, 0.00025);
      //CheckSum3 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      //nSymbols = MMClockRecovery(dataStreamReal, nSamples, dataStreamSymbols, Fs, 10, 0.15);
      nSymbols = MMClockRecovery(dataStreamReal, nSamples, dataStreamSymbols, Fs, 8, 0.15);
      
      //fwrite(dataStreamReal, sizeof(float), nSamples,rawOutFilePtr);
      //fwrite(dataStreamSymbols, sizeof(float), nSymbols,rawOutFilePtr);
      nBits = ManchesterDecode(dataStreamSymbols, nSymbols, dataStreamBits, 1.0);
      nFrames = FindSyncWords(dataStreamBits, nBits, "1110110111100010000", 19, minorFrameFile);      
      
      fwrite(dataStreamBits, sizeof(unsigned char), nBits, rawOutFilePtr);
      totalBits += nBits;
      totalFrames += nFrames;
      totalSymbols += nSymbols;
      totalSamples += nSamples;
      if((((float)( i) / num_samples)*100.0 - percentComplete > 0.15) || feof(waveFilePtr))
         {
         percentComplete = ((float)( i) / num_samples)*100.0;
         printf("\r%0.1f%% %0.3f Ks : %ld Symbols : %ld Bits : %d Frames", ((float)( i) / num_samples)*100.0,(totalSamples)/1000.0, totalSymbols, totalBits, totalFrames);      
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
   