#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
//#include <koolplot.h>
#include "../common/wave.h"
#include "../common/AGC.h"
#include "../common/CarrierTrackPLL.h"
#include "../common/LowPassFilter.h"
//#include "../common/MMClockRecovery.h"
#include "../common/GardenerClockRecovery.h"
#include "../common/ManchesterDecode.h"
#include "../common/ByteSync.h"

#define DSP_GDNR_ERR_LIM            (0.1) //was 0.1
#define DSP_GDNR_GAIN               (3.0) //was 2.5
#define DSP_BAUD                    (400*2)
#define DSP_AGC_ATCK_RATE           (1e-1)//(0.5e-1) //was (1e-1) //attack is when gain is INCREASING (weak signal)
#define DSP_AGC_DCY_RATE            (2e-1) //was (1e-1) //decay is when the gain is REDUCING (strong signal)
#define DSP_PLL_LOCK_ALPA           (0.004)

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

const char *get_filename_ext(const char *filename) 
   {
   const char *dot = strrchr(filename, '.');
   if(!dot || dot == filename)
      return "";
   return dot + 1;
   }
   
int main(int argc, char **argv) 
   {
   FILE *inFilePtr=NULL;
   FILE *rawOutFilePtr=NULL;
   FILE *minorFrameFile=NULL;
   
   HEADER header;
   
   unsigned long chunkSize = 10000, nSamples, i=0, idx, nSymbols, nBits, totalSymbols=0, totalBits=0, totalSamples=0;
   
   double *dataStreamReal, *dataStreamSymbols, *lockSignalStream;
   double Fs;
   double LPF_Fc;   
   double avgPhase, percentComplete=0;
   double normFactor=0;
   
   double *filterCoeffs=NULL, *waveDataTime=NULL;
   double complex *waveData=NULL;
   
   
   //unsigned int CheckSum1=0, CheckSum2=0, CheckSum3=0;
   int LPF_Order, nFrames=0, totalFrames=0,c;
   
   unsigned char *dataStreamBits=NULL;
   char *inFileName=NULL;
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
   printf("Using %ld chunkSize\n",chunkSize);
   LPF_Order = 50;
   filterCoeffs = malloc(sizeof(double) * LPF_Order);   
   inFileName = (char*) malloc(sizeof(char) * 1024);      
   waveData = (double complex*) malloc(sizeof(double complex) * chunkSize);      
   waveDataTime   = (double *) malloc(sizeof(double ) * chunkSize);
   dataStreamReal = (double*) malloc(sizeof(double) * chunkSize);
   lockSignalStream = (double*) malloc(sizeof(double) * chunkSize);
   dataStreamSymbols = (double*) malloc(sizeof(double) * chunkSize);   
   dataStreamBits = (unsigned char*) malloc(sizeof(unsigned char) * chunkSize);
   
   if (dataStreamBits == NULL || 
      filterCoeffs == NULL || 
      inFileName == NULL || 
      waveDataTime == NULL ||
      waveData  == NULL || 
      dataStreamReal == NULL || 
      lockSignalStream == NULL ||
      dataStreamSymbols  == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }    
      
   // get file path
   char cwd[1024];
   if (getcwd(cwd, sizeof(cwd)) != NULL) 
      {      
      strcpy(inFileName, argv[optind]);
      
      // get inFileName from command line
      if (argc < 2)
         {
         printf("No wave file specified\n");
         return 1;
         }
      printf("%s\n", inFileName);
      }

   // open files
   printf("Opening IO files..\n");
   inFilePtr = fopen(inFileName, "rb");
   minorFrameFile = fopen("packets.txt", "w");
  
   if (minorFrameFile == NULL ||
      inFilePtr == NULL)
      {
      printf("Error opening output files\n");
      exit(1);
      }
      
   
   rawOutFilePtr = fopen("output.raw", "wb");
   if (rawOutFilePtr == NULL)
      {
      printf("Error opening output file\n");
      exit(1);
      }
      
   if(strcasecmp(get_filename_ext(inFileName),"wav") == 0)
      header = ReadWavHeader(inFilePtr);
   else
      {
      printf("RAW files not yet supported :(\n");
      exit(1);
      }
   //else if(strcasecmp(get_filename_ext(inFileName),"raw") == 0)
     // header = ReadRawHeader(inFilePtr);
   
   Fs = (double)header.sample_rate;
   long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);
   
   printf("Sample Rate %.2fKHz and %d bits per sample. Total samples %ld\n", Fs/1000.0, header.bits_per_sample ,num_samples);

   LPF_Fc = 700;   
   MakeLPFIR(filterCoeffs, LPF_Order, LPF_Fc, Fs, 1);
   while(!feof(inFilePtr))
      { 
      nSamples = GetComplexWaveChunk(inFilePtr, header, waveData, waveDataTime, chunkSize);
      
      if(i == 0 && normFactor == 0)
         {
         normFactor = StaticGain(waveData, nSamples, 1);
         //normFactor = 1;
         printf("Normalization Factor: %f\n",normFactor);
         }
         
      i+=nSamples;
      /*
      //normalize
      for(idx=0; idx < nSamples; idx++)
         {
          waveData[idx] *= normFactor;         
         }*/
      NormalizingAGCC(waveData, nSamples, normFactor, 0.001);
            
      avgPhase = CarrierTrackPLL(waveData, dataStreamReal, lockSignalStream, nSamples, Fs, 550, 1, DSP_PLL_LOCK_ALPHA, 0.015, 0.015);      
      (void)avgPhase;      
      LowPassFilter(dataStreamReal, nSamples, filterCoeffs, LPF_Order);      
      NormalizingAGC(dataStreamReal, nSamples, DSP_AGC_ATCK_RATE, DSP_AGC_DCY_RATE);      
      //Squelch(dataStreamReal, lockSignalStream, nSamples, 0.25);
      //nSymbols = MMClockRecovery(dataStreamReal, waveDataTime, nSamples, dataStreamSymbols, Fs, 3, 0.15);
      nSymbols = GardenerClockRecovery(dataStreamReal, waveDataTime, nSamples, dataStreamSymbols, Fs, DSP_BAUD, DSP_GDNR_ERR_LIM, DSP_GDNR_GAIN);
      
      //fwrite(dataStreamReal, sizeof(double), nSamples,rawOutFilePtr);
      
      nBits = ManchesterDecode(dataStreamSymbols, waveDataTime, nSymbols, dataStreamBits, 0.5);
      // fwrite(dataStreamBits, sizeof(char), nBits,rawOutFilePtr);
      nFrames = ByteSyncOnSyncword(dataStreamBits, waveDataTime, nBits, "0001011110000", 13, 7, 0, minorFrameFile);      
      
      totalBits += nBits;
      totalFrames += nFrames;
      totalSymbols += nSymbols;
      totalSamples += nSamples;
      if((((double)( i) / num_samples)*100.0 - percentComplete > 0.15) || feof(inFilePtr))
         {
         percentComplete = ((double)( i) / num_samples)*100.0;
         printf("\r");
         //printf("\n");
         printf("%0.1f%% %0.3f Ks : %0.1f Sec: %ld Sym : %ld Bits : %d Frames", ((double)( i) / num_samples)*100.0,(totalSamples)/1000.0, waveDataTime[0], totalSymbols, totalBits, totalFrames);
         }
      
      
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
   fclose(inFilePtr);
   fclose(minorFrameFile);
   
   
   // cleanup before quitting
   free(inFileName);
   free(dataStreamReal);
   free(dataStreamSymbols);
   free(filterCoeffs);
   free(dataStreamReal);
   free(waveData);
   free(dataStreamBits);
   free(waveDataTime);
   
   //quit
   //fflush(stdout);
   return 0;
   } 
   