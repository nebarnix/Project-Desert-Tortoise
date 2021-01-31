#ifndef AGC_H
#define AGC_H

DECIMAL_TYPE FindSignalAmplitude(DECIMAL_TYPE *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE alpha);
void Squelch(DECIMAL_TYPE *dataStream, DECIMAL_TYPE *squelchStreamIn, unsigned long nSamples, DECIMAL_TYPE squelchThreshold);
DECIMAL_TYPE StaticGain(DECIMAL_TYPE complex *complexData,unsigned int nSamples, DECIMAL_TYPE desiredLevel);
void NormalizingAGC(DECIMAL_TYPE *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE initial, DECIMAL_TYPE attack_rate, DECIMAL_TYPE decay_rate);
void NormalizingAGCC(DECIMAL_TYPE complex *dataStreamIn, unsigned long nSamples, DECIMAL_TYPE initial, DECIMAL_TYPE AGC_loop_gain);
#endif