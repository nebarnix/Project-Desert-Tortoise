#ifndef AGC_H
#define AGC_H
void Squelch(double *dataStream, double *squelchStreamIn, unsigned long nSamples, double squelchThreshold);
double StaticGain(double complex *complexData,unsigned int nSamples,double desiredLevel);
void NormalizingAGC(double *dataStreamIn, unsigned long nSamples, double attack_rate, double decay_rate);
void NormalizingAGCC(double complex *dataStreamIn, unsigned long nSamples, double initial, double AGC_loop_gain);
#endif