CC=gcc
CFLAGS=-Wall  -O1 -pipe -fomit-frame-pointer -march=i686
LIBS =

src = $(wildcard *.c)
obj = $(src:.c=.o)

demodPOESTIP: $(obj)
	 $(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(obj) demodPOESTIP