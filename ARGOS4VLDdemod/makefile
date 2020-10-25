CC=gcc
CFLAGS=-Wall -Wno-unused-parameter -O3 -pipe -fomit-frame-pointer -march=i686
LIBS = -lm

vpath %.c ../common

SRC = main.c AGC.c ByteSync.c CarrierTrackingPLL.c GardenerClockRecovery.c LowPassFilter.c ManchesterDecode.c MMClockRecovery.c wave.c 
OBJ = $(SRC:.c=.o)

demodARGOS: $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(OBJ) demodARGOS