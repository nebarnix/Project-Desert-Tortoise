#ifndef CARRIERTRACKPLL_H
#define CARRIERTRACKPLL_H

double CarrierTrackPLL(double complex *complexDataIn, double *realDataOut, unsigned int nSamples, double Fs, double freqRange, double d_lock_threshold, double loopbw_acq, double loopbw_track);

#endif