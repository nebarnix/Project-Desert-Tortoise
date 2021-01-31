# Project Desert Tortoise - Demodulator
## And a tiny app to decode the ground transmitters as well
A demodulator written in C to decode IQ data from the 137.35Mhz/137.77Mhz telemetry beacons from NOAA-15, NOAA-18, and NOAA-19

See the background of this project here
http://wiki.nebarnix.com/wiki/NOAA_POES_TIP_Demodulation

Attribution:
DSP code for the Gardner Clock Recovery based on M&M Clock Recovery. That and the Carrier Tracking PLL adapted from GNURADIO source.
Wave file format parser adapted from http://truelogic.org/wordpress/2015/09/04/parsing-a-wav-file-in-c/ 

Todo in decreasing order of importance:
- [X] Add dynamic lowpass filter coefficient generator to compensate for varying input sample rates
- [X] Add sync word detection for metrics
- [X] Add sync word byte conversion to form minor frame data. 
- [X] Add in 8x interpolation for the clock recovery routine for better performance. Need to optimize interpolator
- [X] Keep track of local recording time and process it in parallel with the data. We will need this through ALL steps of processing
- [X] Add in ability to enable/disable/override things through command line options 
- [X] Add in ability to gauge realtime signal quality
- [X] Add support for reading and writing RAW data files. Need to write read support
- [ ] Add bytesync back into common libs with separate functions for POES and ARGOS
- [X] Add the rest of the processing chain -- See the https://github.com/nebarnix/PDT-TelemetryExplorer repo!!
- [ ] Add a gui to make pretty things happen in realtime (constellation, quality bar graph, pretty waterfall, PLL controls etc)
- [ ] Add hilbert transform to allow for real data input (large bandwidth USB or LSB recording?) This is now a priority for SATNOGS recordings

# POES demodulator 
- Spits out data to minorframes_{datetime}.txt. 
- Takes in an IQ wav file (>50Ksps if the center frequency was not tracked during recording)
- Works better than using an audio pipe

# POES demodulator using port audio
- Spits out data to minorframes_{datetime}.txt. 
- Takes in RAW I/Q audio from sdr# (32khz, unity gain, RAW mode) or similar. 
- Not as good as processing wav files, but can be done in realtime

# ARGOS demodulator
- Demodulates 401.65Mhz transmissions from ground transmitters to POES satellites
- Spits out packets to packets.txt. 
- Takes in .wav files. 
- TODO
- [ ] Add support for real data input
- [ ] Append Quality/SNR for triangulation

# ARGOS demodulator using port audio
- Demodulates 401.65Mhz transmissions from ground transmitters to POES satellites
- Spits out packets to packets.txt. 
- Takes in RAW I/Q audio from sdr# or similar (32khz, unity gain, RAW mode). 
- TODO
- [ ] Add UTC timestamps for ARGOS realtime packets
- [ ] Add support for real data input
- [ ] Append Quality/SNR for triangulation
