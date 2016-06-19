#ifndef BYTESYNC_H
#define BYTESYNC_H
int ByteSyncOnSyncword(unsigned char *bitStreamIn, double *bitStreamInTime, unsigned long nSamples,  char *syncWord, unsigned int syncWordLength, int frameLength, int startBit, FILE *minorFrameFile);
#endif