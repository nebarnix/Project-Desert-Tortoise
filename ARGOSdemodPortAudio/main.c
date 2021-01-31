#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>
#include <conio.h>

#include "config.h" //float vs DECIMAL_TYPE and terminal type set here

//#include "../common/metadata.h"
#include "../common/AGC.h"
#include "../common/CarrierTrackPLL.h"
#include "../common/LowPassFilter.h"
//#include "../common/MMClockRecovery.h"
#include "../common/GardenerClockRecovery.h"
#include "../common/ManchesterDecode.h"
#include "../common/portaudio.h"
#include "ByteSync.h"

#define TRUE 1
#define FALSE 0

//port audio stream params
#define SAMPLE_RATE         (48000)
#define PA_SAMPLE_TYPE      paFloat32
#define NUM_CHANNELS        (2)

#define DEFAULT_CHUNKSIZE  (2400)
//small chunk sizes (<2000) lead to weird time axis problems AND (separately) buffer overflows!
//TODO: Figure out what the heck small chunk sizes mess up
//BUT we want chunk sizes under the minorframe time for meta tracking. 
//There are 10 minor frames per second means chunks should be under Sample rate/10 or 4800

#define DSP_MAX_CARRIER_DEVIATION   (550.0) //was (550.0)
//new values have units of radians per second
//#define DSP_PLL_LOCK_THRESH         (0.1)      //was (1.00) //AR0.1) (doesn't really do anything if the track and acquire gains are the same)
#define DSP_PLL_LOCK_ALPHA          (3.1831)   //(0.004) at 5khz
#define DSP_PLL_ACQ_GAIN            (16)     //0.0015 //was (0.015) at 5khz (16 decoded the best at 5khz) 
#define DSP_PLL_TRCK_GAIN           DSP_PLL_ACQ_GAIN       // 0.0015 //was (0.015) at 5khz
//#define DSP_SQLCH_THRESH            (0.15)     //was (0.25) //was (0.15) //AR0.3 at 5khz

//old values, dependant on sample rate
 #define DSP_PLL_LOCK_THRESH         (0.1) //was (1.00) //AR0.1) (doesn't really do anything if the track and acquire gains are the same)
// #define DSP_PLL_LOCK_ALPHA          (0.004)
// #define DSP_PLL_ACQ_GAIN            (0.015) //0.0015 //was (0.015) 
// #define DSP_PLL_TRCK_GAIN           (0.015) // 0.0015 //was (0.015)
 #define DSP_SQLCH_THRESH            (0.15) //was (0.25) //was (0.15) //AR0.3

//new values have units of radians per second
#define DSP_AGC_ATCK_RATE           (79.5775) //(1e-1)//(0.5e-1) //was (1e-1) //attack is when gain is INCREASING (weak signal)
#define DSP_AGC_DCY_RATE            (159.1549) //(2e-1) //was (1e-1) //decay is when the gain is REDUCING (strong signal)
//old values, dependant on sample rate
//#define DSP_AGC_ATCK_RATE           (1e-1)//(0.5e-1) //was (1e-1) //attack is when gain is INCREASING (weak signal)
//#define DSP_AGC_DCY_RATE            (2e-1) //was (1e-1) //decay is when the gain is REDUCING (strong signal)

#define DSP_AGCC_GAIN               (0.0015) //(0.0005) //was (0.001)
#define DSP_LPF_FC                  (700) //(was (700)
#define DSP_LPF_ORDER               (50) //was (50)

//#define DSP_MM_MAX_DEVIATION        (2.0) //was (3.0)
//#define DSP_MM_GAIN                 (0.2) //was (0.15)
#define DSP_GDNR_ERR_LIM            (0.1) //was 0.1
#define DSP_GDNR_GAIN               (3.0) //was 2.5
#define DSP_BAUD                    (400*2)
#define DSP_MCHSTR_RESYNC_LVL       (0.5) //was (0.5)

#define RAW_OUTPUT_FILES

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
   
   DECIMAL_TYPE *dataStreamReal=NULL, *dataStreamSymbols=NULL, *lockSignalStream=NULL;
   DECIMAL_TYPE Fs = SAMPLE_RATE;
   //DECIMAL_TYPE LPF_Fc;   
   DECIMAL_TYPE averagePhase=0;
   DECIMAL_TYPE normFactor=0;
   DECIMAL_TYPE Time=0;
   
   DECIMAL_TYPE *filterCoeffs=NULL, *waveDataTime=NULL;
   DECIMAL_TYPE complex *waveData=NULL;
   
   
   //unsigned int CheckSum1=0, CheckSum2=0, CheckSum3=0;
   int nFrames=0, totalFrames=0,c;
   
   unsigned char *dataStreamBits=NULL;
   char outFileName[100];
   
   //const char *build_date = __DATE__;
   printf("Project Desert Tortoise: Realtime ARGOS Demodulator by Nebarnix.\nBuild date: %s\n",__DATE__);
   
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
            if(chunkSize != DEFAULT_CHUNKSIZE)
               printf("Override: Using %ld chunkSize\n",chunkSize);
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
      
   if(chunkSize == DEFAULT_CHUNKSIZE)   
      printf("Using default %ld chunkSize\n",chunkSize);
   
   //Allocate the memory we will need
   filterCoeffs = malloc(sizeof(DECIMAL_TYPE) * DSP_LPF_ORDER);        
   waveData = (DECIMAL_TYPE complex*) malloc(sizeof(DECIMAL_TYPE complex) * chunkSize);    
   waveFrame = (float *) malloc(sizeof(float ) * chunkSize*2);
   waveDataTime   = (DECIMAL_TYPE *) malloc(sizeof(DECIMAL_TYPE ) * chunkSize);
   dataStreamReal = (DECIMAL_TYPE*) malloc(sizeof(DECIMAL_TYPE) * chunkSize);
   lockSignalStream = (DECIMAL_TYPE*) malloc(sizeof(DECIMAL_TYPE) * chunkSize);
   dataStreamSymbols = (DECIMAL_TYPE*) malloc(sizeof(DECIMAL_TYPE) * chunkSize);   
   dataStreamBits = (unsigned char*) malloc(sizeof(unsigned char) * chunkSize);
   
   if (dataStreamBits == NULL || 
      filterCoeffs == NULL ||  
      waveDataTime == NULL ||
      waveData  == NULL || 
      dataStreamReal == NULL || 
      lockSignalStream == NULL ||
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
   snprintf(outFileName, 100,"packets_%4d%02d%02d_%02d%02d%02d.txt",tm.tm_year + 1900,tm.tm_mon + 1,tm.tm_mday,tm.tm_hour, tm.tm_min, tm.tm_sec);
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
   
   //LPF_Fc = DSP_LPF_FC;   
   MakeLPFIR(filterCoeffs, DSP_LPF_ORDER, DSP_LPF_FC, Fs, 1);
   
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
         normFactor = StaticGain(waveData, chunkSize, 1);
         //normFactor = 1;
         printf("Normalization Factor: %f\n",normFactor);
         }
         
      i+=chunkSize;
     
      //NormalizingAGCC(waveData, chunkSize, normFactor, DSP_AGCC_GAIN);//this is no longer neccesary because the PLL is now normalizing
      
      //averagePhase = CarrierTrackPLL(waveData, dataStreamReal, lockSignalStream, chunkSize, Fs, DSP_MAX_CARRIER_DEVIATION, DSP_PLL_LOCK_THRESH, DSP_PLL_LOCK_ALPHA, DSP_PLL_ACQ_GAIN, DSP_PLL_TRCK_GAIN);
      averagePhase = CarrierTrackPLL(waveData, dataStreamReal, lockSignalStream, chunkSize, Fs, DSP_MAX_CARRIER_DEVIATION, DSP_PLL_LOCK_THRESH, DSP_PLL_LOCK_ALPHA*(2.0*M_PI/Fs), DSP_PLL_ACQ_GAIN*(2.0*M_PI/Fs), DSP_PLL_TRCK_GAIN*(2.0*M_PI/Fs));
      (void)averagePhase;      
      //fwrite(dataStreamReal, sizeof(DECIMAL_TYPE), chunkSize,rawOutFilePtr);
      
      LowPassFilter(dataStreamReal, chunkSize, filterCoeffs, DSP_LPF_ORDER);      
         
      NormalizingAGC(dataStreamReal, chunkSize,  normFactor, DSP_AGC_ATCK_RATE*(2.0*M_PI/Fs), DSP_AGC_DCY_RATE*(2.0*M_PI/Fs));
      
      #ifdef RAW_OUTPUT_FILES
         fwrite(dataStreamReal, sizeof(DECIMAL_TYPE), chunkSize,rawOutFilePtr);
      #endif
      
      Squelch(dataStreamReal, lockSignalStream, chunkSize, DSP_SQLCH_THRESH);
      //fwrite(dataStreamReal, sizeof(DECIMAL_TYPE), chunkSize,rawOutFilePtr);
      
      //nSymbols = MMClockRecovery(dataStreamReal, waveDataTime, chunkSize, dataStreamSymbols, Fs, DSP_MM_MAX_DEVIATION, DSP_MM_GAIN);
      nSymbols = GardenerClockRecovery(dataStreamReal, waveDataTime, chunkSize, dataStreamSymbols, Fs, DSP_BAUD, DSP_GDNR_ERR_LIM, DSP_GDNR_GAIN);
      
      //#ifdef RAW_OUTPUT_FILES
      //fwrite(dataStreamSymbols, sizeof(DECIMAL_TYPE), nSymbols,rawOutFilePtr);
      //#endif
      
      nBits = ManchesterDecode(dataStreamSymbols, waveDataTime, nSymbols, dataStreamBits, DSP_MCHSTR_RESYNC_LVL);
      
      //#ifdef RAW_OUTPUT_FILES
      //   fwrite(dataStreamBits, sizeof(char), nBits,rawOutFilePtr2);
      //#endif
      
      //nFrames = ByteSyncOnSyncword(dataStreamBits, waveDataTime, nBits, "0001011110000", minorFrameFile);
      nFrames = FindSyncWords(dataStreamBits, waveDataTime, nBits, "0001011110000", 13, minorFrameFile);
      
      totalBits += nBits;
      totalFrames += nFrames;
      totalSymbols += nSymbols;
      totalSamples += chunkSize;

      }
   getch(); //clear the buffer
   //printf("\nChecksum1=%X Checksum2=%X Checksum3=%X", CheckSum1,CheckSum2,CheckSum3);
   //printf("\nAll done! Closing files and exiting.\nENJOY YOUR BITS AND HAVE A NICE DAY\n");
   
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
   if (fclose(minorFrameFile)) { printf("error closing file."); exit(-1); }
   
   if(totalFrames == 0) 
   {
      printf("\n\nNone bits found :(\nRemoving output file and exiting.\nMAY YOU HAVE MORE BETTER BITS ANOTHER DAY\n");
      remove(outFileName); //NO MORE ZERO SIZED FILE LITTER
   }
   else
      printf("\nAll done! Closing files and exiting.\nENJOY YOUR BITS AND HAVE A NICE DAY\n");
   
   // cleanup before quitting
   free(dataStreamSymbols);
   free(filterCoeffs);
   free(dataStreamReal);
   free(waveData);
   free(dataStreamBits);
   free(waveDataTime);
   free(lockSignalStream);
   
   //quit
   //fflush(stdout);
   return 0;
   
   error:
      #ifdef RAW_OUTPUT_FILES
      fclose(rawOutFilePtr);
      fclose(rawOutFilePtr2);
      #endif
      fclose(minorFrameFile);
      remove(outFileName); //NO MORE ZERO SIZED FILE LITTER
      // cleanup before quitting
      
      free(dataStreamReal);
      free(dataStreamSymbols);
      free(filterCoeffs);
      free(waveData);
      free(dataStreamBits);
      free(waveDataTime);
      free(lockSignalStream);
      Pa_Terminate();
      fprintf( stderr, "An error occured while using the portaudio stream\n" );
      fprintf( stderr, "Error number: %d\n", err );
      fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
      free(waveFrame);
      free(waveData);
      return -1;
   } 
   