CC=gcc
CFLAGS=-Wall -Wno-unused-parameter -O3 -pipe -fomit-frame-pointer -march=i686
LIBS = -lm -lportaudio

src = $(wildcard *.c)
obj = $(src:.c=.o)

demodPOESTIP: $(obj)
	 $(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(obj) demodPOESTIP