#ifndef MMCLOCKRECOVERY_H
#define MMCLOCKRECOVERY_H
unsigned long MMClockRecovery(DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamInTime, unsigned long numSamples, DECIMAL_TYPE *dataStreamOut, int Fs, DECIMAL_TYPE baud, DECIMAL_TYPE stepRange, DECIMAL_TYPE kp);
int sign(DECIMAL_TYPE x);
#endif