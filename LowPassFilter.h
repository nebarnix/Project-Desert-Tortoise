#ifndef LOWPASSFILTER_H
#define LOWPASSFILTER_H

void LowPassFilter(float *dataStream, unsigned long nSamples, double *filterCoeffs, int N);
#endif