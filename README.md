# Project Desert Tortoise - Demodulator
A demodulator written in C to decode IQ data from the 137.35Mhz/137.77Mhz telemetry beacons from NOAA-15, NOAA-18, and NOAA-19

See the background of this project here
http://wiki.nebarnix.com/wiki/NOAA_POES_TIP_Demodulation

Attribution:
DSP code for the M&M Clock Recovery and Carrier Tracking PLL adapted from GNURADIO source.
Wave file format parser adapted from http://truelogic.org/wordpress/2015/09/04/parsing-a-wav-file-in-c/ 

Todo in decreasing order of importance:
- Add dynamic lowpass filter coefficient generator to compensate for varying input sample rates
- Add in 8x interpolation for the M&M clock recovery routine for better performance. 
- Add in ability to enable/disable/override things through command line options 
- Add soundcard input capability 
- Add the rest of the processing chain starting with manchester decoding and sync detection
- Add a gui to make pretty things happen in realtime
