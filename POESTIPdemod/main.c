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
   //FILE *rawOutFilePtr=NULL;
   FILE *minorFrameFile=NULL;
   
   HEADER header;
   
   unsigned long chunkSize = 100000, nSamples, i=0, idx, nSymbols, nBits, totalSymbols=0, totalBits=0, totalSamples=0;
   
   double *dataStreamReal=NULL, *dataStreamSymbols=NULL;
   double Fs;
   double LPF_Fc;   
   double PLLLock, percentComplete=0;
   double normFactor=0;
   
   double *filterCoeffs=NULL, *waveDataTime=NULL;
   double complex *waveData=NULL;
   
   
   unsigned int CheckSum1=0, CheckSum2=0, CheckSum3=0;
   int LPF_Order, nFrames=0, totalFrames=0,c;
   
   unsigned char *dataStreamBits=NULL;
   char *filename=NULL;
   //char *cvalue = NULL;   
   
   while ((c = getopt (argc, argv, "n:c:")) != -1)
      {
      switch (c)
         {
         case 'n':
            if(optarg == NULL)
               {
               printf("Static gain unspecified");
               return 1;
               }
            normFactor = atof(optarg);
            printf("Static Gain Override %f\n",normFactor);
            break;
         case 'c':
            if(optarg == NULL)
               {
               printf("Chucksize unspecified");
               return 1;
               }
            chunkSize = atoi(optarg);
            printf("Using %ld chunkSize\n",chunkSize);
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
   filename = (char*) malloc(sizeof(char) * 1024);      
   waveData = (double complex*) malloc(sizeof(double complex) * chunkSize);      
   waveDataTime   = (double *) malloc(sizeof(double ) * chunkSize);
   dataStreamReal = (double*) malloc(sizeof(double) * chunkSize);      
   dataStreamSymbols = (double*) malloc(sizeof(double) * chunkSize);   
   dataStreamBits = (unsigned char*) malloc(sizeof(unsigned char) * chunkSize);
   
   if (dataStreamBits == NULL || 
      filterCoeffs == NULL || 
      filename == NULL || 
      waveDataTime == NULL ||
      waveData  == NULL || 
      dataStreamReal == NULL || 
      dataStreamSymbols  == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }    
      
   // get file path
   char cwd[1024];
   if (getcwd(cwd, sizeof(cwd)) != NULL) 
      {      
      strcpy(filename, argv[optind]);
      
      // get filename from command line
      if (argc < 2)
         {
         printf("No wave file specified\n");
         return 1;
         }
      printf("%s\n", filename);
      }

   // open files
   printf("Opening IO files..\n");
   waveFilePtr = fopen(filename, "rb");
   minorFrameFile = fopen("minorFrame.txt", "w");
  
   if (minorFrameFile == NULL ||
      waveFilePtr == NULL)
      {
      printf("Error opening output files\n");
      exit(1);
      }
      
   /*
   rawOutFilePtr = fopen("output.raw", "wb");
   if (rawOutFilePtr == NULL)
      {
      printf("Error opening output file\n");
      exit(1);
      }*/

   header = ReadWavHeader(waveFilePtr);
   Fs = (double)header.sample_rate;
   long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);
   
   printf("Sample Rate %.2fKHz. Total samples %ld\n", Fs/1000.0,num_samples);

   LPF_Fc = 11000;   
   MakeLPFIR(filterCoeffs, LPF_Order, LPF_Fc, Fs);
   while(!feof(waveFilePtr))
      { 
      nSamples = GetComplexWaveChunk(waveFilePtr, header, waveData, waveDataTime, chunkSize);
      
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
          waveData[idx] *= normFactor;         
         }
      //CheckSum0 += CheckSum((unsigned char*)waveData, sizeof(*waveData) * nSamples);
      PLLLock = CarrierTrackPLL(waveData, dataStreamReal, nSamples, Fs, 4500, 0.1, 0.01, 0.001);      
      //CheckSum1 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      LowPassFilter(dataStreamReal, nSamples, filterCoeffs, LPF_Order);
      //CheckSum2 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      NormalizingAGC(dataStreamReal, nSamples, 0.00025);
      //CheckSum3 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      nSymbols = MMClockRecovery(dataStreamReal, waveDataTime, nSamples, dataStreamSymbols, Fs, 9, 0.15);      
      
      //fwrite(dataStreamReal, sizeof(double), nSamples,rawOutFilePtr);
      
      nBits = ManchesterDecode(dataStreamSymbols, waveDataTime, nSymbols, dataStreamBits, 1.0);
      nFrames = FindSyncWords(dataStreamBits, waveDataTime, nBits, "1110110111100010000", 19, minorFrameFile);      
      
      totalBits += nBits;
      totalFrames += nFrames;
      totalSymbols += nSymbols;
      totalSamples += nSamples;
      if((((double)( i) / num_samples)*100.0 - percentComplete > 0.15) || feof(waveFilePtr))
         {
         percentComplete = ((double)( i) / num_samples)*100.0;
         printf("\r");
         //printf("\n");
         printf("%0.1f%% %0.3f Ks : %0.1f Sec: %ld Sym : %ld Bits : %d Frames", ((double)( i) / num_samples)*100.0,(totalSamples)/1000.0, waveDataTime[0], totalSymbols, totalBits, totalFrames);
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
   //fclose(rawOutFilePtr);
   fclose(waveFilePtr);
   fclose(minorFrameFile);
   
   
   // cleanup before quitting
   free(filename);
   free(dataStreamReal);
   free(dataStreamSymbols);
   free(filterCoeffs);
   free(dataStreamReal);
   free(waveData);
   free(dataStreamBits);
   
   //quit
   //fflush(stdout);
   return 0;
   } 
   