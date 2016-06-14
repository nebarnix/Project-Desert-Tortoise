#ifndef GARDENERCLOCKRECOVERY_H
#define GARDENERCLOCKRECOVERY_H
unsigned long GardenerClockRecovery(double *dataStreamIn,double *dataStreamInTime, unsigned long numSamples, double *dataStreamOut, int Fs, double stepRange, double kp);

#endif