#ifndef WAVE_H
#define WAVE_H

#define TRUE 1
#define FALSE 0

// WAVE file header format
typedef struct HEADER
   {
   unsigned char riff[4];                      // RIFF string
   unsigned int overall_size   ;               // overall size of file in bytes
   unsigned char wave[4];                      // WAVE string
   unsigned char fmt_chunk_marker[4];          // fmt string with trailing null char
   unsigned int length_of_fmt;                 // length of the format data
   unsigned int format_type;                   // format type. 1-PCM, 3- IEEE double, 6 - 8bit A law, 7 - 8bit mu law
   unsigned int channels;                      // no.of channels
   unsigned int sample_rate;                   // sampling rate (blocks per second)
   unsigned int byterate;                      // SampleRate * NumChannels * BitsPerSample/8
   unsigned int block_align;                   // NumChannels * BitsPerSample/8
   unsigned int bits_per_sample;               // bits per sample, 8- 8bits, 16- 16 bits etc
   unsigned char data_chunk_header [4];        // DATA string or FLLR string
   unsigned int data_size;                     // NumSamples * NumChannels * BitsPerSample/8 - size of the next chunk that will be read
   unsigned char type;                         //0-wav, 1-raw
   } HEADER;

void printHeaderInfo(HEADER header);
HEADER ReadWavHeader(FILE *waveFilePtr);
unsigned long int GetComplexRawChunk(FILE *waveFilePtr, HEADER header, DECIMAL_TYPE complex* waveData, DECIMAL_TYPE *waveDataTime, unsigned long int nSamples); //I know, I know, its a bastardization of a wave library but hey....
unsigned long int GetComplexWaveChunk(FILE *waveFilePtr, HEADER header, DECIMAL_TYPE complex* waveData, DECIMAL_TYPE *waveDataTime, unsigned long int nSamples);
unsigned long int GetIQWaveChunk(FILE *waveFilePtr, HEADER header, DECIMAL_TYPE *IData, DECIMAL_TYPE *QData, DECIMAL_TYPE *waveDataTime, unsigned long int nSamples);
char* seconds_to_time(DECIMAL_TYPE raw_seconds);
#endif