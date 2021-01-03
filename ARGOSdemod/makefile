CC=gcc
CFLAGS= -include config.h -Wall -Wno-unused-parameter -Ofast -pipe -fomit-frame-pointer -march=i686
LIBS = -lm

vpath %.c ../common

SRC = main.c AGC.c ByteSync.c CarrierTrackingPLL.c GardenerClockRecovery.c LowPassFilter.c ManchesterDecode.c wave.c MMClockRecovery.c 
OBJ = $(SRC:.c=.o)

demodARGOS: $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(OBJ) demodARGOS