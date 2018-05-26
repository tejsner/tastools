from tastools import load_ill_ascii
import numpy as np
import pylab as plt

class tasdata:
    def __init__(self, filenumbers=None, prefix='', suffix='', axes=[]):
        self.measurements = ['signal', 'monitor']
        self.axes = []
        self.signal = None
        self.monitor = None

        if filenumbers is not None:
            filenames = []
            for f in filenumbers:
                filenames.append(prefix + str(f) + suffix)

            self.add_axis(axes)
            self.load(filenames)
            self.update()

    def __str__(self):
        s = 'monitors: ' + str(self.measurements) + '\n'
        s += 'axes: ' + str(self.axes)
        return s
        
    def add_axis(self, axes):
        """ add an axis to the dataset """
        if isinstance(axes, str):
            axes = [axes]           
            
        for ax in axes:
            self.axes.append(ax)        

    def update(self):
        """ update intensity and error """
        self.I = self.signal/self.monitor
        self.err = np.sqrt(self.signal)/self.monitor

    def load(self, filenames, format='ILL'):
        """ load data and assign signal, monitor and any added axes """
        self.measurements.append('M2')
        self.measurements.append('TIME')
                
        if isinstance(filenames, str):
            filenames = [filenames]

        for meas in self.measurements:
            setattr(self, meas, np.array([]))

        for ax in self.axes:
            setattr(self, ax, np.array([]))
        
        for f in filenames:        
            d = load_ill_ascii(f)
            self.signal = np.concatenate((self.signal, d['CNTS']))
            self.monitor = np.concatenate((self.monitor, d['M1']))
            self.M2 = np.concatenate((self.M2, d['M2']))
            self.TIME = np.concatenate((self.TIME, d['TIME']))

            for ax in self.axes:                
                setattr(self, ax, np.concatenate((getattr(self, ax), d[ax])))

        self.update()

    def plot(self, axis):
        plt.errorbar(getattr(self, axis), self.I, self.err, fmt='o')
        plt.show()

# example

d = tasdata([18530,18531,18532], prefix='data/rawdata/0', axes=['QH', 'QK'])
print(d)
d.plot('QH')
