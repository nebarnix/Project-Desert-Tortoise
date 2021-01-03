#ifndef MANCHESTERDECODE_H
#define MANCHESTERDECODE_H
unsigned long ManchesterDecode(DECIMAL_TYPE *dataStreamIn, DECIMAL_TYPE *dataStreamInTime, unsigned long nSymbols, unsigned char *bitStream, DECIMAL_TYPE resyncThreshold);
#endif