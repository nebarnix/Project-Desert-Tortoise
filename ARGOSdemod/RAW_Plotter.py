# This code plots the 32-bit signed (?) data in the RAW
# file produced by Nebarnix's ARGOS Demodulator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

y = []
x = []
count = 0

with open("output.raw", "rb") as reader:
    while True:
        val = int.from_bytes(reader.read(8), byteorder='little')
        if not val:
            break
        y.append(val)
        x.append(count)
        count += 1

fig, ax = plt.subplots()
ax.plot(x, y)

ax.set(xlabel='sample', ylabel='signal',
       title='output.raw')
ax.grid()

#fig.savefig("test.png")
plt.show()
