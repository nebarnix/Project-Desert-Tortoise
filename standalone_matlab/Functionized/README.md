# Project Desert Tortoise - Matlab Demodulator
A demodulator written in Matlab to decode IQ data from the 137.35Mhz/137.77Mhz telemetry beacons from NOAA-15, NOAA-18, and NOAA-19

## What is in this folder ##
This folder contains the pathfinder matlab project in TWO parts. 
- Part 1 The Pre-processor
* Main File -> preprocess.m
* This takes the IQ wave file and processes it all the way to minor frame data, spacecraft and time data, and parity check

- Part 2 The data processor
*  Main File -> POES.m
* This takes in data from the bitsync in raw, wav, or internal memory from the preprocessor and dissects it into the various instruments and subcommutated streams. 
* The DCS routine is of special interest to Project Desert Tortoise
* You have to make a folder called DCS and make it active or risk never being able to find anything in your active directory ever again. 

See the background of this project here
http://wiki.nebarnix.com/wiki/NOAA_POES_TIP_Demodulation

Attribution:
DSP code for the M&M Clock Recovery and Carrier Tracking PLL adapted from GNURADIO source.

