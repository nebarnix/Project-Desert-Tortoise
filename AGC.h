#ifndef AGC_H
#define AGC_H

float StaticGain(double complex *complexData,unsigned int nSamples,float desiredLevel);
void NormalizingAGC(float *dataStreamIn, unsigned long nSamples, double AGC_loop_gain);
#endif