#ifndef AGC_H
#define AGC_H

double StaticGain(double complex *complexData,unsigned int nSamples,double desiredLevel);
void NormalizingAGC(double *dataStreamIn, unsigned long nSamples, double AGC_loop_gain);
#endif