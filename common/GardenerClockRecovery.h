#ifndef GARDENERCLOCKRECOVERY_H
#define GARDENERCLOCKRECOVERY_H
unsigned long GardenerClockRecovery(DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamInTime, unsigned long numSamples, DECIMAL_TYPE *dataStreamOut, int Fs, DECIMAL_TYPE baud, DECIMAL_TYPE stepRange, DECIMAL_TYPE kp);

#endif