//POES TIP telemetry demodulator - 135.350Mhz and 137.770Mhz
//This version uses IQ data (RAW mode in SDRSharp) over a soundcard. 
//Nebarnix 2020

#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>
#include <conio.h>

#include "config.h" //float vs double and terminal type set here

//#include "../common/metadata.h"
#include "../common/AGC.h"
#include "../common/CarrierTrackPLL.h"
#include "../common/LowPassFilter.h"
//#include "MMClockRecovery.h"
#include "../common/GardenerClockRecovery.h"
#include "../common/ManchesterDecode.h"
#include "../common/portaudio.h"
#include "ByteSync.h"

#define SAMPLE_RATE         (48000)
//#define SAMPLE_RATE         (50000)
#define PA_SAMPLE_TYPE      paFloat32
#define NUM_CHANNELS        (2)

#define TRUE 1
#define FALSE 0
#define DEFAULT_CHUNKSIZE  (2400)

//small chunk sizes (<2000) lead to weird time axis problems AND (separately) buffer overflows!
//TODO: Figure out what the heck small chunk sizes mess up
//BUT we want chunk sizes under the minorframe time for meta tracking. 
//There are 10 minor frames per second means chunks should be under Sample rate/10 or 4800

#define DSP_MAX_CARRIER_DEVIATION   (4500.0) //was (4500)

//new values have units of radians per second
#define DSP_PLL_ACQ_GAIN            198.9437  
#define DSP_PLL_TRCK_GAIN           10.3451
#define DSP_PLL_LOCK_ALPHA          0.3979


//old values, dependant on sample rate
 
//#define DSP_PLL_ACQ_GAIN            (0.025) //(0.005)  
//#define DSP_PLL_TRCK_GAIN           (0.0013) //(0.0015)
//#define DSP_PLL_LOCK_ALPHA          (0.00005)

#define DSP_SQLCH_THRESH            (0.05) //was (0.25)
////#define DSP_PLL_LOCK_THRESH         (0.025) //(0.10) seems to lock on way too easily to nothing at all
#define DSP_PLL_LOCK_THRESH         (0.10)
#define DSP_PLL_SWEEP_THRESH         (0.05) //(0.10) //0.025 //0.05

#define DSP_MM_MAX_DEVIATION        (10.0) //was (3.0)
#define DSP_MM_GAIN                 (0.15)
#define DSP_GDNR_ERR_LIM            (0.1) //was 0.1
#define DSP_GDNR_GAIN               (3.0) //was 2.5
#define DSP_BAUD                    (8320*2+0.3)
#define DSP_MCHSTR_RESYNC_LVL       (0.75) //was (0.5)

//new values have units of radians per second
#define DSP_AGC_ATCK_RATE           (79.5775) //(1e-1)//(0.5e-1) //was (1e-1) //attack is when gain is INCREASING (weak signal)
#define DSP_AGC_DCY_RATE            (159.1549) //(2e-1) //was (1e-1) //decay is when the gain is REDUCING (strong signal)
//#define DSP_AGC_ATCK_RATE           (238.7324) //(1.0e-2)//      //attack is when gain is INCREASING (weak signal)
//#define DSP_AGC_DCY_RATE            (477.4648) //(2.0e-2) //     //decay is when the gain is REDUCING (strong signal)
//old values, dependant on sample rate (TODO: why are these different from the wav version)
//#define DSP_AGC_ATCK_RATE           (1.0e-2)//      //attack is when gain is INCREASING (weak signal)
//#define DSP_AGC_DCY_RATE            (2.0e-2) //     //decay is when the gain is REDUCING (strong signal)

//#define DSP_AGC_GAIN               (.00025)
#define DSP_AGCC_GAIN               (0.00025)
#define DSP_LPF_FC                  (11000.0)//(was (11000)
#define DSP_LPF_INTERP              (3) //interpolation order
#define DSP_LPF_ORDER               (26*DSP_LPF_INTERP) //was (26)

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
   
   DECIMAL_TYPE *waveDataTime=NULL,*dataStreamLPF=NULL,*dataStreamReal=NULL, *dataStreamSymbols=NULL, *dataStreamLPFTime=NULL, *lockSignalStream=NULL;
   DECIMAL_TYPE Fs = SAMPLE_RATE;
   //DECIMAL_TYPE LPF_Fc;   
   DECIMAL_TYPE averagePhase;
   DECIMAL_TYPE normFactor=0;
   DECIMAL_TYPE Time=0;
   
   DECIMAL_TYPE *filterCoeffs=NULL;
   DECIMAL_TYPE complex *waveData=NULL;
   
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
            if(chunkSize != DEFAULT_CHUNKSIZE)
               printf("Override: Using %ld chunkSize\n",chunkSize);
            break;
         case '?':
            /*if (optopt == 'c' || optopt == 'n')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
               fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else*/
            fprintf (stderr, "-n <int> : set initial audio gain\n");   
            fprintf (stderr, "-c <int> : set chunksize\n");
               return 1;
         default:
            fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
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
   
   dataStreamLPF = (DECIMAL_TYPE*) malloc(sizeof(DECIMAL_TYPE) * chunkSize*DSP_LPF_INTERP);
   
   dataStreamLPFTime = (DECIMAL_TYPE*) malloc(sizeof(DECIMAL_TYPE) * chunkSize*DSP_LPF_INTERP);   
   dataStreamSymbols = (DECIMAL_TYPE*) malloc(sizeof(DECIMAL_TYPE) * chunkSize);   
   dataStreamBits = (unsigned char*) malloc(sizeof(unsigned char) * chunkSize);
   
   //verify memory allocation worked right
   if (filterCoeffs == NULL || 
      waveData == NULL ||
      waveFrame == NULL ||
      waveDataTime == NULL ||
      dataStreamReal == NULL ||
      lockSignalStream == NULL ||
      dataStreamLPF == NULL ||
      dataStreamLPFTime == NULL ||
      dataStreamSymbols == NULL ||
      dataStreamBits == NULL) 
      {
      printf("Error in malloc\n");
      exit(1);
      }     
      
   //Initiate audio stream
   err = Pa_Initialize();
   
   //Check audio stream for error and jump to error if neccesary
   if( err != paNoError ) goto error;
   
   //TODO: add an option to list and specify manually an audio device
   inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
   
   if (inputParameters.device == paNoDevice) 
      {
      fprintf(stderr,"Error: No default input device.\n");
      goto error;
      }
      
   inputParameters.channelCount = NUM_CHANNELS;
   inputParameters.sampleFormat = PA_SAMPLE_TYPE;
   inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency ;
   inputParameters.hostApiSpecificStreamInfo = NULL;
   
   err = Pa_IsFormatSupported( &inputParameters, NULL, SAMPLE_RATE );
   if( err == paFormatIsSupported )
      printf( "Sample Rate %dKHz supported\n",SAMPLE_RATE/1000);
   else
   {
      printf( "ERROR: Sample Rate %dKHz not supported!\n",SAMPLE_RATE/1000);
      goto error;
   }

   
   err = Pa_OpenStream(
            &stream,
            &inputParameters,
            NULL, //no output
            SAMPLE_RATE,
            chunkSize,
            0, /* paClipOff, */  // we won't output out of range samples so don't bother clipping them 
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
   //if you open the file here, you might never have anything to write in it, maybe open it when you get your first minorframe?
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
   //since the filter coefficients are dependant on sample rate, we need to synthesize them each time
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
      
      #ifdef RAW_OUTPUT_FILES
         fwrite(waveData, sizeof(DECIMAL_TYPE complex), chunkSize,rawOutFilePtr);
      #endif
      
      //void NormalizingAGC(DECIMAL_TYPE *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE initial, DECIMAL_TYPE attack_rate, DECIMAL_TYPE decay_rate)
      //NormalizingAGC(waveData, chunkSize, normFactor, DSP_AGCC_GAIN);
      //NormalizingAGCC(waveData, chunkSize, normFactor, DSP_AGCC_GAIN); //this is no longer neccesary because the PLL is now normalizing
      
      //CheckSum0 += CheckSum((unsigned char*)waveData, sizeof(*waveData) * nSamples);
      averagePhase = CarrierTrackPLL(waveData, dataStreamReal, lockSignalStream, chunkSize, Fs, DSP_MAX_CARRIER_DEVIATION, DSP_PLL_LOCK_THRESH, DSP_PLL_LOCK_ALPHA*(2.0*M_PI/Fs), DSP_PLL_ACQ_GAIN*(2.0*M_PI/Fs), DSP_PLL_TRCK_GAIN*(2.0*M_PI/Fs)); //DSP_PLL_LOCK_ALPHA*(2.0*M_PI/Fs), DSP_PLL_ACQ_GAIN*(2.0*M_PI/Fs), DSP_PLL_TRCK_GAIN*(2.0*M_PI/Fs));      
      //CheckSum1 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      
      Squelch(dataStreamReal, lockSignalStream, chunkSize, DSP_SQLCH_THRESH); //Squelch if the PLL lock signal is too low, to avoid creeping random sync noise      
      
      //LowPassFilter(dataStreamReal, chunkSize, filterCoeffs, DSP_LPF_ORDER);
      LowPassFilterInterp(waveDataTime, dataStreamReal, dataStreamLPF, dataStreamLPFTime, chunkSize, filterCoeffs, DSP_LPF_ORDER, DSP_LPF_INTERP); 
      
      //CheckSum2 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      //SNR = 10.0 * log10(pow(signalBefore,2) / pow((signalBefore - signalAfter) ,2 )) - SNROffset;
      
      //void NormalizingAGC(DECIMAL_TYPE *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE initial, DECIMAL_TYPE attack_rate, DECIMAL_TYPE decay_rate)
      //NormalizingAGC(dataStreamReal, nSamples, 0.00025);
      //NormalizingAGC(dataStreamLPF, chunkSize*DSP_LPF_INTERP, normFactor, DSP_AGC_ATCK_RATE, DSP_AGC_DCY_RATE);
      NormalizingAGC(dataStreamLPF, chunkSize*DSP_LPF_INTERP, normFactor, DSP_AGC_ATCK_RATE*(2.0*M_PI/(Fs*DSP_LPF_INTERP)), DSP_AGC_DCY_RATE*(2.0*M_PI/(Fs*DSP_LPF_INTERP)) ); //TODO *(2.0*M_PI/(Fs*DSP_LPF_INTERP))
      //CheckSum3 += CheckSum((unsigned char *)dataStreamReal, sizeof(*dataStreamReal) * nSamples);
      
      //#ifdef RAW_OUTPUT_FILES
      //   fwrite(dataStreamLPF, sizeof(DECIMAL_TYPE), chunkSize*DSP_LPF_INTERP,rawOutFilePtr);
      //#endif
      
      //nSymbols = MMClockRecovery(dataStreamLPF, dataStreamLPFTime, chunkSize*DSP_LPF_INTERP, dataStreamSymbols, Fs*DSP_LPF_INTERP, DSP_MM_MAX_DEVIATION, DSP_MM_GAIN);
      nSymbols = GardenerClockRecovery(dataStreamLPF, dataStreamLPFTime, chunkSize*DSP_LPF_INTERP, dataStreamSymbols, Fs*DSP_LPF_INTERP, DSP_BAUD, DSP_GDNR_ERR_LIM, DSP_GDNR_GAIN);
      
      //#ifdef RAW_OUTPUT_FILES
      //   fwrite(dataStreamSymbols, sizeof(DECIMAL_TYPE), nSymbols,rawOutFilePtr);
      //#endif
      
            //ManchesterDecode(DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamInTime, unsigned long nSymbols, unsigned char *bitStream, DECIMAL_TYPE resyncThreshold)
      nBits = ManchesterDecode(dataStreamSymbols, dataStreamLPFTime, nSymbols, dataStreamBits, DSP_MCHSTR_RESYNC_LVL);
      
      #ifdef RAW_OUTPUT_FILES
         fwrite(dataStreamBits, sizeof(char), nBits,rawOutFilePtr2);
      #endif
      
      //it would be nice if this could output some quality information using averagePhase <3
      //this function is a sink, and handles fileIO as well. I don't like that. Maybe it should return complete any minorframes instead??
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
   
   err = Pa_StopStream( stream );
   if( err != paNoError ) goto error;
   
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
   free(dataStreamLPF);
   free(waveData);
   free(dataStreamBits);
   free(waveDataTime);
   free(dataStreamLPFTime);
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
   