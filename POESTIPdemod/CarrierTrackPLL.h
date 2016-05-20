#ifndef CARRIERTRACKPLL_H
#define CARRIERTRACKPLL_H

float CarrierTrackPLL(double complex *complexDataIn, float *realDataOut, unsigned int nSamples, float Fs, float freqRange, float d_lock_threshold, float loopbw_acq, float loopbw_track);

#endif
