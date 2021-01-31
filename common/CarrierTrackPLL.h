#ifndef CARRIERTRACKPLL_H
#define CARRIERTRACKPLL_H

//#define ZeroPhase 1.5708
//#define CONST_CNTR_ANGLE 1.5708 //this is just pi/2 ?!
#define CONST_CNTR_ANGLE (M_PI/2.0)
//we are looking for 60.7 degrees
//#define CONST_PERFECT 1.16937


DECIMAL_TYPE CarrierTrackPLL(DECIMAL_TYPE complex *complexDataIn, DECIMAL_TYPE *realDataOut, DECIMAL_TYPE *lockSignalStreamOut, unsigned int nSamples, DECIMAL_TYPE Fs, DECIMAL_TYPE freqRange, DECIMAL_TYPE d_lock_threshold, DECIMAL_TYPE lockSigAlpha, DECIMAL_TYPE loopbw_acq, DECIMAL_TYPE loopbw_track);

#endif