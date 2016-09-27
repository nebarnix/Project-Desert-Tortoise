#ifndef LOWPASSFILTER_H
#define LOWPASSFILTER_H

void LowPassFilter(double *dataStream, unsigned long nSamples, double *filterCoeffs, int N);
void LowPassFilterInterp(double *dataStreamInTime, double *dataStreamIn, double *dataStreamOut, double *dataStreamOutTime, unsigned long nSamples, double *filterCoeffs, int N, int interpFactor);
int MakeLPFIR(double *h, int N, double Fc, double Fs, int interpFactor);
#endif