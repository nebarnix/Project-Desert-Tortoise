# This code plots the 32-bit double data in the RAW
# file produced by Nebarnix's ARGOS Demodulator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import struct

y = []
x = []
count = 0
numBytes = 8

if (len(sys.argv) >= 3):
    numBytes = int(sys.argv[2])

if numBytes == 8:
    dataType = 'd'
if numBytes == 4:
    dataType = 'f'

with open(sys.argv[1], "rb") as reader:
    
    reader.seek(0,2) #Jumps to the end
    endLocation = reader.tell() #Give you the end location (characters from start)
    reader.seek(0)   #Jump to the beginning of the file again
    
    while True:
        if numBytes == 1:
            val = int.from_bytes(reader.read(numBytes), byteorder='little')
        else:
            val = struct.unpack(dataType,reader.read(numBytes))
        y.append(val)
        x.append(count)
        count += 1
        if reader.tell() == endLocation:
            break

fig, ax = plt.subplots()
ax.plot(x, y)

ax.set(xlabel='sample', ylabel='signal',
       title=sys.argv[1])
ax.grid()

#fig.savefig("test.png")
plt.show()
