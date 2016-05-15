#ifndef MMCLOCKRECOVERY_H
#define MMCLOCKRECOVERY_H
unsigned long MMClockRecovery(float *dataStreamIn, unsigned long numSamples, float *dataStreamOut, int Fs, float stepRange, float kp);
int sign(float x);
#endif