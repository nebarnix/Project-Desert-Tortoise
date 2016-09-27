#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../common/wave.h"
#include "../common/AGC.h"
#include "../common/CarrierTrackPLL.h"
#include "../common/LowPassFilter.h"
//#include "../common/MMClockRecovery.h"
#include "../common/GardenerClockRecovery.h"
#include "../common/ManchesterDecode.h"
#include "ByteSync.h"

#define TRUE 1
#define FALSE 0
#define DEFAULT_CHUNKSIZE  (10000)

#define DSP_MAX_CARRIER_DEVIATION   (4500.0) //was (4500)

//Wave file 
//#define DSP_PLL_LOCK_THRESH         (0.1) //(0.10)
//#define DSP_PLL_ACQ_GAIN            (0.01) //(0.005)  
//#define DSP_PLL_TRCK_GAIN           (0.001) //(0.0015)

//Port audio
#define DSP_PLL_LOCK_THRESH         (0.025) //(0.10)
#define DSP_PLL_ACQ_GAIN            (0.025) //(0.005)  
#define DSP_PLL_TRCK_GAIN           (0.0013) //(0.0015)
#define DSP_PLL_LOCK_ALPHA          (0.00005) //
#define DSP_SQLCH_THRESH            (0.01) //was (0.25)
#define DSP_MM_MAX_DEVIATION        (10.0) //was (3.0)
#define DSP_MM_GAIN                 (0.15)
#define DSP_GDNR_ERR_LIM            (0.1) //was 0.1
#define DSP_GDNR_GAIN               (3.0) //was 2.5
#define DSP_BAUD                    (8320*2+0.3)
#define DSP_MCHSTR_RESYNC_LVL       (0.75) //was (0.5)
#define DSP_AGC_ATCK_RATE           (10.0e-3)//      //attack is when gain is INCREASING (weak signal)
#define DSP_AGC_DCY_RATE            (20.0e-3) //     //decay is when the gain is REDUCING (strong signal)
//#define DSP_AGC_GAIN               (.00025)
#define DSP_AGCC_GAIN               (0.00025)
#define DSP_LPF_FC                  (11000.0)//(was (11000)
#define DSP_LPF_INTERP              (3) //interpolation order
#define DSP_LPF_ORDER               (26*DSP_LPF_INTERP) //was (26)
//#define DSP_LPF_ORDER               (26) //was (26)

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

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
   FILE *waveFilePtr;
   FILE *rawOutFilePtr;
   FILE *minorFrameFile;
   
   HEADER header;
   
   unsigned long idx,chunkSize = DEFAULT_CHUNKSIZE, nSamples, i=0, nSymbols, nBits, totalSymbols=0, totalBits=0, totalSamples=0;
   
   double tempVal;
   double Fs;
   double averagePhase, percentComplete=0;
   double normFactor=0;
      
   double *waveDataTime=NULL,*dataStreamLPF=NULL,*dataStreamReal=NULL, *dataStreamSymbols=NULL, *dataStreamLPFTime=NULL;
   
   double *filterCoeffs=NULL;
   double complex *waveData=NULL;
   char qualityString[20];
   
   unsigned int CheckSum0=0, CheckSum1=0, CheckSum2=0, CheckSum3=0;
   int nFrames=0, totalFrames=0,c;
   
   unsigned char *dataStreamBits=NULL;
   char *filename=NULL;
   char outFileName[100];
   char outputRawFiles=0;   
   
   const char *build_date = __DATE__;
   printf("Project Desert Tortoise: Wave file NOAA TIP Demodulator by Nebarnix.\nBuild date: %s\n",build_date);
   
   while ((c = getopt (argc, argv, "rn:c:")) != -1)
      {
      switch (c)
         {
         case 'r':
            outputRawFiles = 1;
            printf("Outputting Debugging Raw Files\n");
            break;            
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
   
   filterCoeffs = malloc(sizeof(double) * DSP_LPF_ORDER);
   filename = (char*) malloc(sizeof(char) * 1024);      
   waveData = (double complex*) malloc(sizeof(double complex) * chunkSize);      
   waveDataTime   = (double *) malloc(sizeof(double ) * chunkSize);
   dataStreamReal = (double*) malloc(sizeof(double) * chunkSize);
   dataStreamLPF = (double*) malloc(sizeof(double) * chunkSize*DSP_LPF_INTERP);
   dataStreamLPFTime = (double*) malloc(sizeof(double) * chunkSize*DSP_LPF_INTERP);   
   dataStreamSymbols = (double*) malloc(sizeof(double) * chunkSize);   
   dataStreamBits = (unsigned char*) malloc(sizeof(unsigned char) * chunkSize);
   
   if (dataStreamBits == NULL || 
      filterCoeffs == NULL || 
      filename == NULL || 
      waveDataTime == NULL ||
      waveData  == NULL || 
      dataStreamReal == NULL || 
      dataStreamLPF == NULL ||
      dataStreamLPFTime == NULL ||
      dataStreamSymbols  == NULL)
      {
      printf("Error in malloc\n");
      exit(1);
      }    
      
   // get file path
   char cwd[1024];
   if (getcwd(cwd, sizeof(cwd)) != NULL) 
      {      
            
      // get filename from command line
      if (argc < 2)
         {
         printf("No wave file specified\n");
         return 1;
         }
      strcpy(filename, argv[optind]);
      printf("%d %s\n",argc, filename);
      }

   // open files
   printf("Opening IO files..\n");
   waveFilePtr = fopen(filename, "rb");
   time_t t = time(NULL);
   struct tm tm = *localtime(&t);
   //printf("now: %d-%d-%d %d:%d:%d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
   snprintf(outFileName, 100,"minorFrames_%4d%02d%02d_%02d%02d%02d.txt",tm.tm_year + 1900,tm.tm_mon + 1,tm.tm_mday,tm.tm_hour, tm.tm_min, tm.tm_sec);
   minorFrameFile = fopen(outFileName, "w");
  
   if (minorFrameFile == NULL ||
      waveFilePtr == NULL)
      {
      printf("Error opening output files\n");
      exit(1);
      }
      
   if(outputRawFiles == 1)
      {
      rawOutFilePtr = fopen("output.raw", "wb");
      if (rawOutFilePtr == NULL)
         {
         printf("Error opening output file\n");
         exit(1);
         }
      }

   header = ReadWavHeader(waveFilePtr);
   Fs = (double)header.sample_rate;
   long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);
   
   printf("Sample Rate %.2fKHz and %d bits per sample. Total samples %ld\n", Fs/1000.0, header.bits_per_sample ,num_samples);
    
   MakeLPFIR(filterCoeffs, DSP_LPF_ORDER, DSP_LPF_FC, Fs*DSP_LPF_INTERP, DSP_LPF_INTERP);
   
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
      //for(idx=0; idx < nSamples; idx++)
      //   {
      //    waveData[idx] *= normFactor;         
      //   }
      CheckSum0 += CheckSum((unsigned char*)waveData, sizeof(*waveData) * nSamples);
      NormalizingAGCC(waveData, nSamples, normFactor, DSP_AGCC_GAIN);
      
      /*
      if(outputRawFiles == 1)
         {
         //fwrite(waveData, sizeof(complex double), nSamples,rawOutFilePtr);
         for(idx=0; idx < nSamples; idx++)
            {
            tempVal = creal(waveData[idx]);
            fwrite(&tempVal, sizeof(double), 1,rawOutFilePtr);
            tempVal = cimag(waveData[idx]);
            fwrite(&tempVal, sizeof(double), 1,rawOutFilePtr);
            }
         }*/
      
      CheckSum1 += CheckSum((unsigned char*)waveData, sizeof(*waveData) * nSamples);
      //averagePhase = CarrierTrackPLL(waveData, dataStreamReal, NULL, nSamples, Fs, 4500, 0.1, 0.01, 0.001);      
      averagePhase = CarrierTrackPLL(waveData, dataStreamReal, NULL, nSamples, Fs, DSP_MAX_CARRIER_DEVIATION, DSP_PLL_LOCK_THRESH, DSP_PLL_LOCK_ALPHA, DSP_PLL_ACQ_GAIN, DSP_PLL_TRCK_GAIN);
      CheckSum2 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      
      
      //LowPassFilter(dataStreamReal, nSamples, filterCoeffs, DSP_LPF_ORDER);
      LowPassFilterInterp(waveDataTime, dataStreamReal, dataStreamLPF, dataStreamLPFTime, nSamples, filterCoeffs, DSP_LPF_ORDER, DSP_LPF_INTERP);
      
      /*if(outputRawFiles == 1)         
         fwrite(dataStreamLPF, sizeof(double), nSamples*DSP_LPF_INTERP,rawOutFilePtr);*/
      
      CheckSum3 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      
      //NormalizingAGC(dataStreamReal, nSamples, 0.00025);
      NormalizingAGC(dataStreamLPF, nSamples*DSP_LPF_INTERP, DSP_AGC_ATCK_RATE, DSP_AGC_DCY_RATE);
      /*if(outputRawFiles == 1)         
         fwrite(dataStreamLPF, sizeof(double), nSamples*DSP_LPF_INTERP,rawOutFilePtr);*/
      //CheckSum3 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      //nSymbols = MMClockRecovery(dataStreamReal, waveDataTime, nSamples, dataStreamSymbols, Fs, 9, 0.15);      
      
      nSymbols = GardenerClockRecovery(dataStreamLPF, dataStreamLPFTime, nSamples*DSP_LPF_INTERP, dataStreamSymbols, Fs*DSP_LPF_INTERP, DSP_BAUD, DSP_GDNR_ERR_LIM, DSP_GDNR_GAIN);            
      if(outputRawFiles == 1)         
         fwrite(dataStreamLPF, sizeof(double), nSymbols,rawOutFilePtr);
      
      nBits = ManchesterDecode(dataStreamSymbols, waveDataTime, nSymbols, dataStreamBits, 1.0);
      //int ByteSyncOnSyncword(unsigned char *bitStreamIn, double *bitStreamInTime, unsigned long nSamples,  char *syncWord, unsigned int syncWordLength, FILE *minorFrameFile);
      nFrames = ByteSyncOnSyncword(dataStreamBits, waveDataTime, nBits, "1110110111100010000", 19, minorFrameFile);      
      
      totalBits += nBits;
      totalFrames += nFrames;
      totalSymbols += nSymbols;
      totalSamples += nSamples;
      if((((double)( i) / num_samples)*100.0 - percentComplete > 0.15) || feof(waveFilePtr))
         {
         percentComplete = ((double)( i) / num_samples)*100.0;
         printf("\r");
         //printf("\n");
         averagePhase = 10.0 * log10( pow(1.5708 - averagePhase,2));
         
         if(averagePhase > -4.3)
            snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_GREEN,averagePhase,ANSI_COLOR_RESET);
         else if(averagePhase > -5)
            snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_YELLOW,averagePhase,ANSI_COLOR_RESET);
         else if(averagePhase > -6)
            snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_YELLOW,averagePhase,ANSI_COLOR_RESET);
         else
            snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_RED,averagePhase,ANSI_COLOR_RESET);
         
         printf("%0.1f%% %0.3f Ks : %0.1f Sec: %ld Sym : %ld Bits : %d Frames : %s   ", ((double)( i) / num_samples)*100.0,(totalSamples)/1000.0, waveDataTime[0], totalSymbols, totalBits, totalFrames, qualityString);
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
   
   printf("\nChecksum0=%X Checksum1=%X Checksum2=%X Checksum3=%X", CheckSum0, CheckSum1,CheckSum2,CheckSum3);
   printf("\nAll done! Closing files and exiting.\nENJOY YOUR BITS AND HAVE A NICE DAY\n");
   
   if(outputRawFiles == 1)   
      fclose(rawOutFilePtr);
   
   if (fclose(waveFilePtr)) { printf("error closing file."); exit(-1); }
   if (fclose(minorFrameFile)) { printf("error closing file."); exit(-1); }
   
   // cleanup before quitting
   free(filename);
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
   