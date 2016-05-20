#ifndef MMCLOCKRECOVERY_H
#define MMCLOCKRECOVERY_H
unsigned long MMClockRecovery(double *dataStreamIn,double *dataStreamInTime, unsigned long numSamples, double *dataStreamOut, int Fs, double stepRange, double kp);
int sign(double x);
#endif