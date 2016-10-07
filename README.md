# Project Desert Tortoise - Demodulator
## And a tiny app to decode the ground transmitters as well
A demodulator written in C to decode IQ data from the 137.35Mhz/137.77Mhz telemetry beacons from NOAA-15, NOAA-18, and NOAA-19

See the background of this project here
http://wiki.nebarnix.com/wiki/NOAA_POES_TIP_Demodulation

Attribution:
DSP code for the M&M Clock Recovery and Carrier Tracking PLL adapted from GNURADIO source.
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
- [X] Add a gui to make pretty things happen in realtime -- See above except its not really in realtime so nevermind but whatever.
- [ ] Add hilbert transform to allow for real data input (large bandwidth USB or LSB recording?)

# POES demodulator using port audio
- Spits out data to minorframes_{datetime}.txt. 
- Takes in RAW I/Q audio from sdr# (32khz, unity gain, RAW mode) or similar. 
- Works! [NEW!]

# ARGOS demodulator
- Demodulates 401.65Mhz transmissions from ground transmitters to POES satellites
- Spits out packets to packets.txt. 
- Takes in wave files. 
- Mostly works! (NEW!)
- TODO
- [ ] Add support for RAW data files
- [ ] Add support for real data input
- [ ] Figure out how to calculate Signal to Noise Ratio which is important for triangulation

# ARGOS demodulator using port audio
- Demodulates 401.65Mhz transmissions from ground transmitters to POES satellites
- Spits out packets to packets.txt. 
- Takes in RAW I/Q audio from sdr# or similar (32khz, unity gain, RAW mode) . 
- Mostly works! (NEW!)
- TODO
- [ ] Add UTC timestamps for ARGOS realtime packets
- [ ] Add support for real data input
- [ ] Figure out how to calculate Signal to Noise Ratio which is important for triangulation
