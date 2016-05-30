#ifndef LOWPASSFILTER_H
#define LOWPASSFILTER_H

void LowPassFilter(double *dataStream, unsigned long nSamples, double *filterCoeffs, int N);
int MakeLPFIR(double *h, int N, double Fc, double Fs);
#endif