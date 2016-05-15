CC=gcc
CFLAGS=-Wall -O2
LIBS = -lkoolplot -lgdi32 -lcomdlg32 -luuid -loleaut32 -lole32 -lstdc++ -lsupc++

src = $(wildcard *.c)
obj = $(src:.c=.o)

demodPOESTIP: $(obj)
	 $(CC) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	del -f $(obj) demodPOESTIP