from tastools import load_ill_ascii
import numpy as np
import pylab as plt

class tasdata:
    def __init__(self, filenumbers=None, prefix='', suffix='', axes=[]):
        self.measurements = ['signal', 'monitor']
        self.axes = []
        self.signal = None
        self.monitor = None
        
        # if given axes, asign the first one as our x-axis
        if axes:
            self.xaxis = axes[0]
        else:
            self.xaxis = False

        # if given filenumbers, load data right away
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
        if self.xaxis:
            self.x = getattr(self, self.xaxis)

    def load_ill(self, filenames):
        """ loads data from ILL TAS """
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

            # if we have not defined an x-axis, take the scan axis from the first file
            if not self.axes:
                self.axes.append(d.dtype.names[1])
                setattr(self, self.axes[0], np.array([]))
                self.xaxis = self.axes[0]
            
            for ax in self.axes:
                setattr(self, ax, np.concatenate((getattr(self, ax), d[ax])))


    def load(self, filenames, fileformat='ILL'):
        """ load data and assign signal, monitor and any added axes """
        if fileformat == 'ILL':
            self.load_ill(filenames)

        self.update()

    def bin(self, bins, axis=None, method='avg'):
        """ Bins the data along a certain axis 
            Either specify a binsize (float) or a vector (array or list) with the actual bins 
            methods: 'avg' takes the average of x-values in each bin.
                     'center' uses the bin center as x-value

            TODO: bin center, save binned x-axis in a seperate variabe
        """
        # use the scan-axis if not specified
        if axis is None:
            axis = self.xaxis

        if isinstance(bins, float):
            xmin = min(getattr(self, axis))
            xmax = max(getattr(self, axis))
            bins = np.arange(xmin-bins/2, xmax+bins, bins)
        
        inds = np.digitize(getattr(self, axis), bins)

        # make lists of temporary axes (x) and measurements (y) with each element being initialized as an array
        x = []
        y = []
        for ax in self.axes:
            x.append(np.array([]))

        for meas in self.measurements:
            y.append(np.array([]))
        
        # fill up the temporary axes and measurments
        for i in np.unique(inds):
            # for axes we take the mean value in each bin
            for j, ax in enumerate(self.axes):
                x[j] = np.append(x[j], getattr(self, ax)[inds == i].mean())

            # for measurements we add up the values
            for k, meas in enumerate(self.measurements):
                y[k] = np.append(y[k], getattr(self, meas)[inds == i].sum())

        # update axes and mreasurements
        for i, ax in enumerate(self.axes):
            setattr(self, ax, x[i])

        for j, meas in enumerate(self.measurements):
            setattr(self, meas, y[j])

        self.update()

    def plot(self):
        plt.errorbar(self.x, self.I, self.err, fmt='o')

    def __sub__(self, other):
        """ subtraction given that the two objects have intersecting bins """
        pass

    def __add__(self, other):
        """ addition given that the two objects have intersecting bins """
        pass

# example
d_40K = tasdata([18543,18578], prefix='data/rawdata/0')
d_80K = tasdata([18548], prefix='data/rawdata/0')
d_10K = tasdata([18536], prefix='data/rawdata/0')
d_40K.bin(0.007, 'QH')
d_80K.bin(0.007, 'QH')
d_10K.bin(0.007, 'QH')

# d_10K.plot()
# d_40K.plot()
# d_80K.plot()
plt.plot(d_10K.x, d_10K.I-d_40K.I+0.005, 'o')
plt.plot(d_10K.x, len(d_10K.x)*[0.005], '-')
plt.show()