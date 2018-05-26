from tastools import load_ill_ascii
import numpy as np
import pylab as plt

class tasdata:
    def __init__(self):
        self.fields = ['signal', 'monitor']
        self.signal = None
        self.monitor = None

    def __str__(self):
        if self.signal is None:
            return 'empty tasdata object'

        s = 'tasdata object with the following fields of length ' + str(len(self.signal)) + ':'
        for f in self.fields:
            s += '\n' + str(f) + ': ' + str(getattr(self, f))        
        return s
        
    def add_axis(self, axis):
        """ add an axis to the dataset """
        self.fields.append(axis)

    def load(self, filename, format='ILL', prefix=''):
        """ load data and asign signal, monitor and any added axes """
        d = load_ill_ascii(prefix + filename)
        self.signal = d['CNTS']
        self.monitor = d['M1']

        for f in self.fields[2::]:
            setattr(self, f, d[f])

    def plot(self, axis):
        plt.errorbar(getattr(self, axis), self.signal, np.sqrt(self.signal), fmt='o')
        plt.show()

## example

data = 'data/rawdata/018530'

d = tasdata()
d.add_axis('QH')
d.load(data)
print(d)
d.plot('QH')