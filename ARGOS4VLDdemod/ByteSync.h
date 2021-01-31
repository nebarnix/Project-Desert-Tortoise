#ifndef BYTESYNC_H
#define BYTESYNC_H
int FindSyncWords(unsigned char *bitStreamIn, double *bitStreamInTime, unsigned long nSamples,  char *syncWord, unsigned int syncWordLength, FILE *minorFrameFile);
#endif