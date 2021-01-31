CC=gcc
CFLAGS = -include config.h -Wall -Wno-unused-parameter -Ofast -pipe -fomit-frame-pointer -march=i686
#CFLAGS = -Wall -Wno-unused-parameter -Ofast -pipe -fomit-frame-pointer -march=corei7 -msse4 -flto
#CFLAGS = -Wall -Wno-unused-parameter -Ofast -pipe -fomit-frame-pointer -march=native -msse2
#CFLAGS = -g 

LIBS = -lm

vpath %.c ../common

SRC = main.c AGC.c ByteSync.c CarrierTrackingPLL.c GardenerClockRecovery.c LowPassFilter.c ManchesterDecode.c MMClockRecovery.c wave.c #atan2_approx.c
OBJ = $(SRC:.c=.o)

demodPOES: $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(OBJ) demodPOES