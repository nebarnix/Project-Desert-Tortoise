CC=gcc
CFLAGS=-Wall  -O1 -pipe -fomit-frame-pointer 
LIBS = -lm

src = $(wildcard *.c)
obj = $(src:.c=.o)

demodPOESTIP: $(obj)
	 $(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(obj) demodPOESTIP
