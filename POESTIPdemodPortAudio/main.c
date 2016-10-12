#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>

#if defined (__WIN32__)
  #include <windows.h>
  #ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
     #define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
  #endif
#endif

#include "../common/AGC.h"
#include "../common/CarrierTrackPLL.h"
#include "../common/LowPassFilter.h"
//#include "MMClockRecovery.h"
#include "../common/GardenerClockRecovery.h"
#include "../common/ManchesterDecode.h"
#include "../common/portaudio.h"
#include "ByteSync.h"

//#define SAMPLE_RATE         (48000)
#define SAMPLE_RATE         (50000)
#define PA_SAMPLE_TYPE      paFloat32
#define NUM_CHANNELS        (2)
#define DEFAULT_CHUNKSIZE  (1000)

#define DSP_MAX_CARRIER_DEVIATION   (4500.0) //was (4500)
#define DSP_PLL_LOCK_THRESH         (0.025) //(0.10)
#define DSP_PLL_ACQ_GAIN            (0.025) //(0.005)  
#define DSP_PLL_TRCK_GAIN           (0.0013) //(0.0015)
#define DSP_PLL_LOCK_ALPHA          (0.00005)
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

#define TRUE 1
#define FALSE 0

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define QUALITY_SHIT -20 
#define QUALITY_LOW -6
#define QUALITY_MEDIUM -5
#define QUALITY_GOOD -4.3

//#define RAW_OUTPUT_FILES

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
   //PortAduio Variables
   PaStreamParameters inputParameters;
   PaStream *stream;
   PaError err;
   
   //Files we will use
   #ifdef RAW_OUTPUT_FILES
      FILE *rawOutFilePtr=NULL;
      FILE *rawOutFilePtr2=NULL;
   #endif
   
   FILE *minorFrameFile=NULL;
   
   unsigned long chunkSize = DEFAULT_CHUNKSIZE, i=0, idx, idx2, nSymbols, nBits, totalSymbols=0, totalBits=0, totalSamples=0;
   
   float *waveFrame; //data from input stream
   float realPart, imagPart;
   
   double *waveDataTime=NULL,*dataStreamLPF=NULL,*dataStreamReal=NULL, *dataStreamSymbols=NULL, *dataStreamLPFTime=NULL;
   double Fs = SAMPLE_RATE;
   double LPF_Fc;   
   double averagePhase;
   double normFactor=0;
   double Time=0;
   
   double *filterCoeffs=NULL;
   double complex *waveData=NULL;
   char qualityString[20];
   
   unsigned int CheckSum1=0, CheckSum2=0, CheckSum3=0;
   int nFrames=0, totalFrames=0,c;
   
   unsigned char *dataStreamBits=NULL;
   char outFileName[100];   
   
   //const char *build_date = __DATE__;
   //enable ANSI color formatting for (newer) windows builds!
   #if defined (__WIN32__)
      SetConsoleMode(GetStdHandle(STD_OUTPUT_HANDLE), ENABLE_VIRTUAL_TERMINAL_PROCESSING | ENABLE_PROCESSED_OUTPUT | ENABLE_WRAP_AT_EOL_OUTPUT);
   #endif
   printf("Project Desert Tortoise: Realtime NOAA TIP Demodulator by Nebarnix.\nBuild date: %s\n",__DATE__);
   printf("Please tune SDR input to NOAA-15/18 137.350Mhz or NOAA-19 137.770Mhz\n");
   printf("Please set bandwidth to >=32khz, filter order lowest possible, RAW output mode\n");
   
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
      
   printf("Using %ld chunkSize\n",chunkSize);
   
   
   //Allocate the memory we will need
   filterCoeffs = malloc(sizeof(double) * DSP_LPF_ORDER);         
   waveData = (double complex*) malloc(sizeof(double complex) * chunkSize);      
   waveFrame = (float *) malloc(sizeof(float ) * chunkSize*2);
   waveDataTime   = (double *) malloc(sizeof(double ) * chunkSize);
   dataStreamReal = (double*) malloc(sizeof(double) * chunkSize);
   dataStreamLPF = (double*) malloc(sizeof(double) * chunkSize*DSP_LPF_INTERP);
   dataStreamLPFTime = (double*) malloc(sizeof(double) * chunkSize*DSP_LPF_INTERP);   
   dataStreamSymbols = (double*) malloc(sizeof(double) * chunkSize);   
   dataStreamBits = (unsigned char*) malloc(sizeof(unsigned char) * chunkSize);
   
   if (dataStreamBits == NULL || 
      filterCoeffs == NULL ||  
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
      
   //Initiate audio stream
   err = Pa_Initialize();
   
   if( err != paNoError ) goto error;
   
   inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
   
   if (inputParameters.device == paNoDevice) 
      {
      fprintf(stderr,"Error: No default input device.\n");
      goto error;
      }
      
   inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
   inputParameters.channelCount = NUM_CHANNELS;
   inputParameters.sampleFormat = PA_SAMPLE_TYPE;
   inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency ;
   inputParameters.hostApiSpecificStreamInfo = NULL;
   
   err = Pa_OpenStream(
            &stream,
            &inputParameters,
            NULL,
            SAMPLE_RATE,
            chunkSize,
            0, /* paClipOff, */  /* we won't output out of range samples so don't bother clipping them */
            NULL,
            NULL );
   
   if( err != paNoError ) goto error;
   
   err = Pa_StartStream( stream );
   if( err != paNoError ) goto error;  

   // open files
   printf("Opening Output files..\n");
   
   
   time_t t = time(NULL);
   struct tm tm = *localtime(&t);
   //printf("now: %d-%d-%d %d:%d:%d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
   snprintf(outFileName, 100,"minorFrames_%4d%02d%02d_%02d%02d%02d.txt",tm.tm_year + 1900,tm.tm_mon + 1,tm.tm_mday,tm.tm_hour, tm.tm_min, tm.tm_sec);
   minorFrameFile = fopen(outFileName, "w");
  
   if (minorFrameFile == NULL)
      {
      printf("Error opening output files\n");
      exit(1);
      }
      
   #ifdef RAW_OUTPUT_FILES
   rawOutFilePtr = fopen("output.raw", "wb");
   if (rawOutFilePtr == NULL)
      {
      printf("Error opening output file\n");
      exit(1);
      }
      
   rawOutFilePtr2 = fopen("outputbits.raw", "wb");
   if (rawOutFilePtr2 == NULL)
      {
      printf("Error opening output file\n");
      exit(1);
      }
   #endif
   
   LPF_Fc = DSP_LPF_FC;   
   MakeLPFIR(filterCoeffs, DSP_LPF_ORDER, DSP_LPF_FC, Fs*DSP_LPF_INTERP, DSP_LPF_INTERP);

   
while(!kbhit())
      { 
      err = Pa_ReadStream(stream, waveFrame, chunkSize);
      if(err != paNoError)
         {
         printf("overflow :(\n");
         //goto error;
         }
         
      idx2 = 0;
      //convert from stereo IQ to complex stream
      for(idx = 0; idx < chunkSize*2; idx+=2)
         {
         realPart = waveFrame[idx];
         imagPart = waveFrame[idx+1];
         waveData[idx2] = realPart + imagPart * I;
         Time += (1/Fs);
         waveDataTime[idx2] = Time;
         idx2++;
         //if(idx % 100 == 0) printf("\r%6.3f\t\t\t%6.3f",realPart, imagPart);
         }
      
      if(i == 0 && normFactor == 0)
         {
         normFactor = StaticGain(waveData, chunkSize, 1.0);
         //normFactor = 0.000041;
         printf("Normalization Factor: %f\n",normFactor);
         }
         
      i+=chunkSize;
      
      NormalizingAGCC(waveData, chunkSize, normFactor, DSP_AGCC_GAIN);
      
      //CheckSum0 += CheckSum((unsigned char*)waveData, sizeof(*waveData) * nSamples);
      averagePhase = CarrierTrackPLL(waveData, dataStreamReal, NULL, chunkSize, Fs, DSP_MAX_CARRIER_DEVIATION, DSP_PLL_LOCK_THRESH, DSP_PLL_LOCK_ALPHA, DSP_PLL_ACQ_GAIN, DSP_PLL_TRCK_GAIN);      
      //CheckSum1 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      
      //LowPassFilter(dataStreamReal, chunkSize, filterCoeffs, DSP_LPF_ORDER);
      LowPassFilterInterp(waveDataTime, dataStreamReal, dataStreamLPF, dataStreamLPFTime, chunkSize, filterCoeffs, DSP_LPF_ORDER, DSP_LPF_INTERP);
      
      //CheckSum2 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      //SNR = 10.0 * log10(pow(signalBefore,2) / pow((signalBefore - signalAfter) ,2 )) - SNROffset;
      //NormalizingAGC(dataStreamReal, nSamples, 0.00025);
      NormalizingAGC(dataStreamLPF, chunkSize*DSP_LPF_INTERP, DSP_AGC_ATCK_RATE, DSP_AGC_DCY_RATE);
      //CheckSum3 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      
      //#ifdef RAW_OUTPUT_FILES
      //   fwrite(dataStreamLPF, sizeof(double), chunkSize*DSP_LPF_INTERP,rawOutFilePtr);
      //#endif
      
      //nSymbols = MMClockRecovery(dataStreamLPF, dataStreamLPFTime, chunkSize*DSP_LPF_INTERP, dataStreamSymbols, Fs*DSP_LPF_INTERP, DSP_MM_MAX_DEVIATION, DSP_MM_GAIN);
      nSymbols = GardenerClockRecovery(dataStreamLPF, dataStreamLPFTime, chunkSize*DSP_LPF_INTERP, dataStreamSymbols, Fs*DSP_LPF_INTERP, DSP_BAUD, DSP_GDNR_ERR_LIM, DSP_GDNR_GAIN);
      
      #ifdef RAW_OUTPUT_FILES
         fwrite(dataStreamSymbols, sizeof(double), nSymbols,rawOutFilePtr);
      #endif
      
            //ManchesterDecode(double *dataStreamIn, double *dataStreamInTime, unsigned long nSymbols, unsigned char *bitStream, double resyncThreshold)
      nBits = ManchesterDecode(dataStreamSymbols, dataStreamLPFTime, nSymbols, dataStreamBits, DSP_MCHSTR_RESYNC_LVL);
      
      #ifdef RAW_OUTPUT_FILES
         fwrite(dataStreamBits, sizeof(char), nBits,rawOutFilePtr2);
      #endif
      
      //it would be nice if this could output some quality information using averagePhase <3
      nFrames = ByteSyncOnSyncword(dataStreamBits, dataStreamLPFTime, nBits, "1110110111100010000", 19, minorFrameFile);      
      
      totalBits += nBits;
      totalFrames += nFrames;
      totalSymbols += nSymbols;
      totalSamples += chunkSize;
      
      //Print updates     
      printf("\r");
      
      averagePhase = 10.0 * log10( pow(1.5708 - averagePhase,2));
      
      if(averagePhase > QUALITY_GOOD)
         snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_GREEN,averagePhase,ANSI_COLOR_RESET);
      else if(averagePhase > QUALITY_MEDIUM)
         snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_YELLOW,averagePhase,ANSI_COLOR_RESET);
      else if(averagePhase > QUALITY_LOW)
         snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_YELLOW,averagePhase,ANSI_COLOR_RESET);
      else //QUALITY_SHIT
         snprintf(qualityString, 20,"%s%02.1fQ%s", ANSI_COLOR_RED,averagePhase,ANSI_COLOR_RESET);
      
      //SUPRESS OUTPUT HERE
      printf("%0.3f Ks : %0.1f Sec: %ld Sym : %ld Bits : %d Frames : %s   ",(totalSamples)/1000.0, waveDataTime[0], totalSymbols, totalBits, totalFrames, qualityString);
   
      
      //fwrite("2", sizeof(unsigned char), 1, rawOutFilePtr);
      /*for(idx=0; idx < nSamples; idx++)
         {
         fVal = (crealf(waveData[i]));
         fwrite(&fVal,sizeof(fVal),1,rawOutFilePtr);
         fVal = (cimagf(waveData[i]));
         fwrite(&fVal,sizeof(fVal),1,rawOutFilePtr);
         }*/
      }
   
   getch(); //clear the buffer
   //printf("\nChecksum1=%X Checksum2=%X Checksum3=%X", CheckSum1,CheckSum2,CheckSum3);
   printf("\nAll done! Closing files and exiting.\nENJOY YOUR BITS AND HAVE A NICE DAY\n");
   
   err = Pa_CloseStream( stream );
   if( err != paNoError ) goto error;
   
   Pa_Terminate();
   free(waveFrame);
   free(waveData);
   //return 0;
 
    
   #ifdef RAW_OUTPUT_FILES
      fclose(rawOutFilePtr);
      fclose(rawOutFilePtr2);
   #endif
   
   fclose(minorFrameFile);
   
   
   // cleanup before quitting
   free(dataStreamSymbols);
   free(filterCoeffs);
   free(dataStreamReal);
   free(dataStreamLPF);
   free(waveData);
   free(dataStreamBits);
   free(waveDataTime);
   
   //quit
   //fflush(stdout);
   return 0;
   
   error:
      #ifdef RAW_OUTPUT_FILES
         fclose(rawOutFilePtr);
         fclose(rawOutFilePtr2);
      #endif
      
      fclose(minorFrameFile);
      
      // cleanup before quitting
      
      free(dataStreamReal);
      free(dataStreamSymbols);
      free(filterCoeffs);
      free(dataStreamLPF);
      free(waveData);
      free(dataStreamBits);
      free(waveDataTime);
      Pa_Terminate();
      fprintf( stderr, "An error occured while using the portaudio stream\n" );
      fprintf( stderr, "Error number: %d\n", err );
      fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
      free(waveFrame);
      free(waveData);
      return -1;
   } 
   