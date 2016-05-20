#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "wave.h"

void printHeaderInfo(HEADER header)
   {
   printf("(1-4): %.4s \n", header.riff);
   //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
   printf("(5-8) Overall size: bytes:%u, Kb:%u \n", header.overall_size, header.overall_size/1024);
   printf("(9-12) Wave marker: %.7s\n", header.wave);
   printf("(13-16) Fmt marker: %.3s\n", header.fmt_chunk_marker);
   //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
   printf("(17-20) Length of Fmt header: %u \n", header.length_of_fmt);
   char format_name[10] = "";
   
   if (header.format_type == 1)
      strcpy(format_name,"PCM");
   else if (header.format_type == 6)
      strcpy(format_name, "A-law");
   else if (header.format_type == 7)
      strcpy(format_name, "Mu-law");
   
   printf("(21-22) Format type: %u %s \n", header.format_type, format_name);
   //printf("%u %u \n", buffer2[0], buffer2[1]);
   printf("(23-24) Channels: %u \n", header.channels);
   //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
   printf("(25-28) Sample rate: %u\n", header.sample_rate);
   //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
   printf("(29-32) Byte Rate: %u , Bit Rate:%u\n", header.byterate, header.byterate*8);
   //printf("%u %u \n", buffer2[0], buffer2[1]);
   printf("(33-34) Block Alignment: %u \n", header.block_align);
   //printf("%u %u \n", buffer2[0], buffer2[1]);
   printf("(35-36) Bits per sample: %u \n", header.bits_per_sample);
   printf("(37-40) Data Marker: %.4s \n", header.data_chunk_header);
   //printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
   printf("(41-44) Size of data chunk: %u \n", header.data_size);
   // calculate no.of samples
   long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);   
   printf("Number of samples:%lu \n", num_samples);
   
   long size_of_each_sample = (header.channels * header.bits_per_sample) / 8;
   printf("Size of each sample:%ld bytes\n", size_of_each_sample);
   
   long bytes_in_each_channel = (size_of_each_sample / header.channels);
   printf("Size of each channel:%ld bytes\n", bytes_in_each_channel);
   
   // calculate duration of file
   double duration_in_seconds = (double) header.overall_size / header.byterate;
   printf("Approx.Duration in seconds=%f\n", duration_in_seconds);
   printf("Approx.Duration in h:m:s=%s\n", seconds_to_time(duration_in_seconds));
   
   }

//Read nSamples of complex data from the wave file and store in array. Actual number of returned bytes is returned as an integer.
//Advances file pointer, ready to read more new samples
int GetComplexWaveChunk(FILE *waveFilePtr, HEADER header, double complex* waveData, double *waveDataTime, int nSamples)
   {
   //Check to make sure file is actually ready for reading
   if (waveFilePtr == NULL)
      {
      printf("Error opening file\n");
      exit(1);
      }
   if ( waveData == NULL || waveDataTime == NULL)
      {
      printf("Dude, allocate your fracking memory already. UGH. \n");
      exit(1);
      }  
   if(header.channels != 2)
      {
      printf("Complex read requires 2 channels (I and Q)\n");
      exit(1);
      }
   // read each sample from data chunk if PCM
   if (header.format_type != 1)
      {
      printf("Only PCM is currently supported :(\n");
      exit(1);
      } 
   
   long i =0;
   long size_of_each_sample = (header.channels * header.bits_per_sample) / 8;
   unsigned char data_buffer[size_of_each_sample];
   int  size_is_correct = TRUE;
   int read;
   static double time=0, Ts=0;
   double realVal, imagVal;
   int16_t data_in_channel = 0;
   
   if(Ts == 0)
      Ts = 1.0/(double)header.sample_rate;
   
   // make sure that the bytes-per-sample is completely divisible by num.of channels
   long bytes_in_each_channel = (size_of_each_sample / header.channels);
   if ((bytes_in_each_channel  * header.channels) != size_of_each_sample)
      {
      printf("Error: %ld x %ud <> %ld\n", bytes_in_each_channel, header.channels, size_of_each_sample);
      size_is_correct = FALSE;
      }
   
   if (size_is_correct) 
      {
      /*switch (header.bits_per_sample)
            {
            case 8:
               maxsize = 128;
               break;
               
            case 16:
               
               maxsize = 32768;
               break;
               
            case 32:
               
               maxsize = 2147483648;
               break;
            }*/
            
      for (i =0; i < nSamples; i++) 
         {
         //printf("==========Sample %ld / %ld=============\n", i, nSamples   );
         read = fread(data_buffer, sizeof(data_buffer), 1, waveFilePtr);
         if (read == 1) 
            {            
            // dump the data read
            unsigned int  xchannels = 0;
            
            
            for (xchannels = 0; xchannels < header.channels; xchannels ++ ) 
               {
               //printf("Channel#%d : ", (xchannels+1));
               // convert data from little endian to big endian based on bytes in each channel sample
               if (bytes_in_each_channel == 4) 
                  {
                  data_in_channel =   data_buffer[0] | (data_buffer[1]<<8) | (data_buffer[2]<<16) | (data_buffer[3]<<24);
                  }                                               
               else if (bytes_in_each_channel == 2) 
                  {
                  data_in_channel = data_buffer[xchannels*bytes_in_each_channel+0] | (data_buffer[xchannels*bytes_in_each_channel+1] << 8);
                  }
               else if (bytes_in_each_channel == 1) 
                  {
                  data_in_channel = data_buffer[0];
                  }               
               
               //return normalized complex data
               if(xchannels == 0) //Real channel is first
                  {
                  realVal =  data_in_channel/32768.0;                  
                  }
               else //Imaginary channel is second
                  {
                  imagVal = data_in_channel/32768.0;                  
                  waveData[i] = realVal + imagVal * I;
                  time += Ts;
                  waveDataTime[i] = time;
                  //waveData[i] = imagVal + realVal * I;
                  }
               }
            
            //printf("\n");
            }
         else 
            {
            //EOF!           
            return i;
            }            
         } //    for (i =1; i <= nSamples; i++) {         
      } //    if (size_is_correct) {         
   return nSamples;
   }

//READ HEADER FILE FUNCTION
//Inputs file pointer to freshly opened file
//Populates header struct with details from file
//Leaves file pointer at the end of the header, ready for reading samples
HEADER ReadWavHeader(FILE *waveFilePtr)
   {
   HEADER header;
   
   unsigned char buffer4[4];
   unsigned char buffer2[2];
   //Check to make sure file is actually ready for reading
   if (waveFilePtr == NULL)
      {
      printf("Error opening file\n");
      exit(1);
      }
   
   //clear header struct to all zeros
   memset(&header, 0, sizeof header);
         
   //printf("Filepointer position: %x\n", ftell(waveFilePtr));
   //printf("Filepointer reference: %x\n",(unsigned int)waveFilePtr);
   //int read;
   
   // read header parts
   //unsigned char test[5];
   //fread(test, 4, 1, waveFilePtr);
   //test[4] = '\0';   
   //printf("(1-4): %s \n", test);
   
   fread(header.riff, sizeof(header.riff), 1, waveFilePtr);
   
   fread(buffer4, sizeof(buffer4), 1, waveFilePtr);
   
   // convert little endian to big endian 4 byte int
   header.overall_size  = buffer4[0] | (buffer4[1]<<8) | (buffer4[2]<<16) | (buffer4[3]<<24);

   fread(header.wave, sizeof(header.wave), 1, waveFilePtr);
   
   fread(header.fmt_chunk_marker, sizeof(header.fmt_chunk_marker), 1, waveFilePtr);
   
   fread(buffer4, sizeof(buffer4), 1, waveFilePtr);
   
   // convert little endian to big endian 4 byte integer
   header.length_of_fmt = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);
   
   fread(buffer2, sizeof(buffer2), 1, waveFilePtr); 
   //printf("%u %u \n", buffer2[0], buffer2[1]);
   
   header.format_type = buffer2[0] | (buffer2[1] << 8);
   
  
   fread(buffer2, sizeof(buffer2), 1, waveFilePtr);
  
   header.channels = buffer2[0] | (buffer2[1] << 8);
   
   fread(buffer4, sizeof(buffer4), 1, waveFilePtr);
   
   header.sample_rate = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);
  
   fread(buffer4, sizeof(buffer4), 1, waveFilePtr);
   
   header.byterate  = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);
   
   fread(buffer2, sizeof(buffer2), 1, waveFilePtr);
   
   header.block_align = buffer2[0] | (buffer2[1] << 8);
   
   fread(buffer2, sizeof(buffer2), 1, waveFilePtr);
  
   header.bits_per_sample = buffer2[0] | (buffer2[1] << 8);
   
   fread(header.data_chunk_header, sizeof(header.data_chunk_header), 1, waveFilePtr);
   
   fread(buffer4, sizeof(buffer4), 1, waveFilePtr);
  
   header.data_size = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24 );
     
   return header;
   }
   
   /**
 * Convert seconds into hh:mm:ss format
 * Params:
 *  seconds - seconds value
 * Returns: hms - formatted string
 **/
char* seconds_to_time(double raw_seconds)
   {
   char *hms;
   int hours, hours_residue, minutes, seconds, milliseconds;
   hms = (char*) malloc(100);
   
   sprintf(hms, "%f", raw_seconds);
   
   hours = (int) raw_seconds/3600;
   hours_residue = (int) raw_seconds % 3600;
   minutes = hours_residue/60;
   seconds = hours_residue % 60;
   milliseconds = 0;
   
   // get the decimal part of raw_seconds to get milliseconds
   char *pos;
   pos = strchr(hms, '.');
   int ipos = (int) (pos - hms);
   char decimalpart[15];
   memset(decimalpart, ' ', sizeof(decimalpart));
   strncpy(decimalpart, &hms[ipos+1], 3);
   milliseconds = atoi(decimalpart);
   
   
   sprintf(hms, "%d:%d:%d.%d", hours, minutes, seconds, milliseconds);
   return hms;
   }