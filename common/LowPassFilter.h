#ifndef LOWPASSFILTER_H
#define LOWPASSFILTER_H

void LowPassFilter(DECIMAL_TYPE *dataStream, unsigned long nSamples, DECIMAL_TYPE *filterCoeffs, int N);
void LowPassFilterInterp(DECIMAL_TYPE *dataStreamInTime, DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamOut, DECIMAL_TYPE *dataStreamOutTime, unsigned long nSamples, DECIMAL_TYPE *filterCoeffs, int N, int interpFactor);
int MakeLPFIR(DECIMAL_TYPE *h, int N, DECIMAL_TYPE Fc, DECIMAL_TYPE Fs, int interpFactor);
#endif