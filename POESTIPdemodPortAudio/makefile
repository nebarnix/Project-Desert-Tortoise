CC=gcc
//CFLAGS=-Wall -Wno-unused-parameter -Ofast -pipe -fomit-frame-pointer -march=i686
CFLAGS = -include config.h -Wall -Wno-unused-parameter -Ofast -pipe -fomit-frame-pointer -march=i686
LIBS = -lm -lportaudio

vpath %.c ../common

SRC = main.c AGC.c ByteSync.c CarrierTrackingPLL.c GardenerClockRecovery.c LowPassFilter.c ManchesterDecode.c MMClockRecovery.c
OBJ = $(SRC:.c=.o)

demodPOESP: $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(OBJ) demodPOESP