CC=gcc
CFLAGS=-Wall -Wno-unused-parameter -O3 -pipe -fomit-frame-pointer -march=i686
LIBS = -lm

src = $(wildcard *.c)
obj = $(src:.c=.o)

demodARGOS: $(obj)
	 $(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(obj) demodARGOS