#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>
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
#define DEFAULT_CHUNKSIZE  (1000)

#define DSP_MAX_CARRIER_DEVIATION   (250.0) //was (550.0)

#define DSP_PLL_LOCK_THRESH         (0.5) //was (1.00) //AR0.1) (doesn't really do anything if the track and acquire gains are the same)
#define DSP_PLL_LOCK_ALPHA          (0.004)
#define DSP_PLL_ACQ_GAIN            (0.015) //0.0015 //was (0.015)
#define DSP_PLL_TRCK_GAIN           (0.015) // 0.0015 //was (0.015)
#define DSP_SQLCH_THRESH            (0.15) //was (0.25) //was (0.15) //AR0.3

#define DSP_GDNR_ERR_LIM            (0.1) //was 0.1
#define DSP_GDNR_GAIN               (3.0) //was 2.5
#define DSP_BAUD                    (200*2) //was (400*2) // ARGOS 4 PTT-VLD-A4 (A4-SS-TER-SP-0079-CNES) is 200bps
#define DSP_MCHSTR_RESYNC_LVL       (1.0) //was (0.5)
#define DSP_AGC_ATCK_RATE           (1e-1) //was (1e-1) //(0.5e-1) //was (1e-1) //attack is when gain is INCREASING (weak signal)
#define DSP_AGC_DCY_RATE            (2e-1) // was (2e-1) //was (1e-1) //decay is when the gain is REDUCING (strong signal)

#define DSP_AGCC_GAIN               (0.001) //(0.0005) //was (0.001) //AR0.0015
#define DSP_LPF_FC                  (350) //(was (700)
#define DSP_LPF_ORDER               (50) //was (50)

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
  //Wave variable
  HEADER header;

  //Files we will use
  FILE *inFilePtr=NULL;
  FILE *minorFrameFile=NULL;
  FILE *rawOutFilePtr=NULL;

  unsigned long chunkSize = DEFAULT_CHUNKSIZE, nSamples, i=0, idx, nSymbols, nBits, totalSymbols=0, totalBits=0, totalSamples=0;

  double *dataStreamReal=NULL, *dataStreamSymbols=NULL, *lockSignalStream=NULL;
  double Fs;
  double LPF_Fc;
  double averagePhase, percentComplete=0;
  double normFactor=0;

  double *filterCoeffs=NULL, *waveDataTime=NULL;
  double complex *waveData=NULL;

  double dspManchesterResyncLevel = DSP_MCHSTR_RESYNC_LVL;
  double dspMaxCarrierDeviation = DSP_MAX_CARRIER_DEVIATION;

  //unsigned int CheckSum1=0, CheckSum2=0, CheckSum3=0;
  int nFrames=0, totalFrames=0,c;

  unsigned char *dataStreamBits=NULL;
  char *inFileName=NULL;
  char outFileName[100];
  char outputRawFiles=0;

  //const char *build_date = __DATE__;
  printf("Project Desert Tortoise: Realtime ARGOS Demodulator by Nebarnix.\nBuild date: %s\n",__DATE__);

  while ((c = getopt (argc, argv, "rd:m:n:c:")) != -1)
    {
    switch (c)
      {
      case 'r':
      outputRawFiles = 1;
      printf("Outputting Debugging Raw Files\n");
      break;
      case 'd':
      if(optarg == NULL)
        {
        printf("DSP Maximum Carrier Deviation unspecified");
        return 1;
        }
      dspMaxCarrierDeviation = atof(optarg);
      printf("DSP Maximum Carrier Deviation Override %f\n",dspMaxCarrierDeviation);
      break;
      case 'm':
      if(optarg == NULL)
        {
        printf("DSP Manchester Resync Level unspecified");
        return 1;
        }
      dspManchesterResyncLevel = atof(optarg);
      printf("DSP Manchester Resync Level Override %f\n",dspManchesterResyncLevel);
      break;
      case 'n':
      if(optarg == NULL)
        {
        printf("Static Gain unspecified");
        return 1;
        }
      normFactor = atof(optarg);
      printf("Static Gain Override %f\n",normFactor);
      break;
      case 'c':
      if(optarg == NULL)
        {
        printf("Chunksize unspecified");
        return 1;
        }
      chunkSize = atoi(optarg);
      printf("Using %ld chunkSize\n",chunkSize);
      break;
      case '?':
      if (optopt == 'c' || optopt == 'n' || optopt == 'd' || optopt == 'm')
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
  //printf("Using %ld chunkSize\n",chunkSize);

  //Allocate the memory we will need
  filterCoeffs = malloc(sizeof(double) * DSP_LPF_ORDER);
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
    waveData == NULL ||
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

    // get inFileName from command line
    if (argc < 2)
      {
      printf("No wave file specified\n");
      return 1;
      }
    strcpy(inFileName, argv[optind]);
    printf("%s\n", inFileName);
    }

  // open files
  printf("Opening IO files..\n");
  inFilePtr = fopen(inFileName, "rb");

  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  //printf("now: %d-%d-%d %d:%d:%d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  snprintf(outFileName, 100,"packets_%4d%02d%02d_%02d%02d%02d.txt",tm.tm_year + 1900,tm.tm_mon + 1,tm.tm_mday,tm.tm_hour, tm.tm_min, tm.tm_sec);
  minorFrameFile = fopen(outFileName, "w");

  if (minorFrameFile == NULL ||
    inFilePtr == NULL)
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

  //printHeaderInfo(header);

  MakeLPFIR(filterCoeffs, DSP_LPF_ORDER, DSP_LPF_FC, Fs, 1);

  long filePos = ftell(inFilePtr); // Store the file position

  if (normFactor == 0) // If the user has not overridden the static gain (normalization factor)
  {
    // Calculate the normalization factor based on the whole file - not just the first chunk (in case it is quiet)
    normFactor = 1000000;
    while(!feof(inFilePtr))
      {
      nSamples = GetComplexWaveChunk(inFilePtr, header, waveData, waveDataTime, chunkSize);
      double chunkNormFactor = StaticGain(waveData, nSamples, 1.0);
      if(chunkNormFactor < normFactor)
        {
        normFactor = chunkNormFactor;
        }
      }
    printf("Normalization Factor: %f\n",normFactor);
  }

  fseek(inFilePtr, filePos, SEEK_SET); // Rewind the file

  while(!feof(inFilePtr))
    {
    nSamples = GetComplexWaveChunk(inFilePtr, header, waveData, waveDataTime, chunkSize);

    //if(outputRawFiles == 1) fwrite(waveData, sizeof(double), nSamples, rawOutFilePtr);
    //if(outputRawFiles == 1) fwrite(waveDataTime, sizeof(double), nSamples, rawOutFilePtr);

    if(i == 0 && normFactor == 0)
      {
      // Calculate the static gain (normalization factor) based on the first chunk
      normFactor = StaticGain(waveData, nSamples, 1.0);
      //normFactor = 1;
      printf("Normalization Factor: %f\n",normFactor);
      }

    i+=nSamples;

    NormalizingAGCC(waveData, nSamples, normFactor, DSP_AGCC_GAIN);

    //if(outputRawFiles == 1) fwrite(waveData, sizeof(double), nSamples, rawOutFilePtr);

    averagePhase = CarrierTrackPLL(waveData, dataStreamReal, lockSignalStream, nSamples, Fs, dspMaxCarrierDeviation, DSP_PLL_LOCK_THRESH, DSP_PLL_LOCK_ALPHA, DSP_PLL_ACQ_GAIN, DSP_PLL_TRCK_GAIN);

    //if(outputRawFiles == 1) fwrite(dataStreamReal, sizeof(double), nSamples, rawOutFilePtr);
    //if(outputRawFiles == 1) fwrite(lockSignalStream, sizeof(double), nSamples, rawOutFilePtr);

    (void)averagePhase;

    LowPassFilter(dataStreamReal, nSamples, filterCoeffs, DSP_LPF_ORDER);

    //if(outputRawFiles == 1) fwrite(dataStreamReal, sizeof(double), nSamples, rawOutFilePtr);

    NormalizingAGC(dataStreamReal, nSamples, DSP_AGC_ATCK_RATE, DSP_AGC_DCY_RATE);

    //if(outputRawFiles == 1) fwrite(dataStreamReal, sizeof(double), nSamples, rawOutFilePtr);

    Squelch(dataStreamReal, lockSignalStream, nSamples, DSP_SQLCH_THRESH);

    //if(outputRawFiles == 1) fwrite(dataStreamReal, sizeof(double), nSamples, rawOutFilePtr);

    //nSymbols = MMClockRecovery(dataStreamReal, waveDataTime, nSamples, dataStreamSymbols, Fs, 3, 0.15);
    nSymbols = GardenerClockRecovery(dataStreamReal, waveDataTime, nSamples, dataStreamSymbols, Fs, DSP_BAUD, DSP_GDNR_ERR_LIM, DSP_GDNR_GAIN);

    if(outputRawFiles == 1)
      {
      for(idx=0; idx < nSymbols; idx++)
        {
        fwrite(&dataStreamSymbols[idx],sizeof(double),1,rawOutFilePtr);
        }
      }

    //nBits = ManchesterDecode(dataStreamSymbols, waveDataTime, nSymbols, dataStreamBits, DSP_MCHSTR_RESYNC_LVL);
    nBits = ManchesterDecode(dataStreamSymbols, waveDataTime, nSymbols, dataStreamBits, dspManchesterResyncLevel);

    // if(outputRawFiles == 1)
    //   {
    //   for(idx=0; idx < nBits; idx++)
    //     {
    //     fwrite(&dataStreamBits[idx],sizeof(char),1,rawOutFilePtr);
    //     }
    //   }

    // There's no need to invert the data. FindSyncWords automatically checks for the sync word and its inverse.

    // ARGOS 4 PTT-VLD-A4 A4-SS-TER-SP-0079-CNES defines a sync pattern of 0xAC5353
    nFrames = FindSyncWords(dataStreamBits, waveDataTime, nBits, "101011000101001101010011", 24, minorFrameFile);

    totalBits += nBits;
    totalFrames += nFrames;
    totalSymbols += nSymbols;
    totalSamples += nSamples;
    if((((double)( i) / num_samples)*100.0 - percentComplete > 0.15) || feof(inFilePtr))
      {
      percentComplete = ((double)( i) / num_samples)*100.0;
      printf("\r");
      printf("\n");
      printf("%0.1f%% %0.3f Ks : %0.1f Sec: %ld Sym : %ld Bits : %d Packets  ", ((double)( i) / num_samples)*100.0, (totalSamples)/1000.0, waveDataTime[0], totalSymbols, totalBits, totalFrames);
      }

    }

  //printf("\nChecksum1=%X Checksum2=%X Checksum3=%X", CheckSum1,CheckSum2,CheckSum3);
  printf("\nAll done! Closing files and exiting.\nENJOY YOUR BITS AND HAVE A NICE DAY\n");
  printf("\nARGOS 4 PTT-VLD-A4 data is convolutionally encoded, punctured and interleaved.\nOnly the first two bits (message length) are human-readable.\n");

  if(outputRawFiles == 1)
    fclose(rawOutFilePtr);

  if (fclose(inFilePtr)) { printf("error closing file."); exit(-1); }
  if (fclose(minorFrameFile)) { printf("error closing file."); exit(-1); }


  // cleanup before quitting
  free(inFileName);
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
  }
