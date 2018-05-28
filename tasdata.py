from tastools import load_ill_ascii
import numpy as np
import pylab as plt

class tasdata:
    def __init__(self, filenumbers=None, prefix='', suffix='', axes=[]):
        self.monitors = ['signal', 'monitor'] # we allways have at least one signal and monitor
        self.axes = []
        self.binned = False
        
        # if given axes, asign the first one as our x-axis
        if axes:
            self.xlabel = axes[0]
        else:
            self.xlabel = False

        # if given filenumbers, load data right away
        if filenumbers is not None:
            filenames = []
            for f in filenumbers:
                filenames.append(prefix + str(f) + suffix)

            self.add_axis(axes)
            self.load(filenames)
            self.update()

    def __str__(self):
        s = 'monitors: ' + str(self.monitors) + '\n'
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
        if self.xlabel:
            self.x = getattr(self, self.xlabel)

    def load_ill(self, filenames):
        """ loads data from ILL TAS """
        # we allways have the following additional monitors when doing ILL TAS experiments
        self.monitors.append('M2')
        self.monitors.append('TIME')
                
        if isinstance(filenames, str):
            filenames = [filenames]

        for mon in self.monitors:
            setattr(self, mon, np.array([]))

        for ax in self.axes:
            setattr(self, ax, np.array([]))
        
        for f in filenames:        
            d = load_ill_ascii(f)            
            self.signal = np.concatenate((self.signal, d['CNTS']))
            self.monitor = np.concatenate((self.monitor, d['M1']))
            
            for mon in self.monitors[2:]:
                setattr(self, mon, np.concatenate((getattr(self, mon), d[mon])))

            # if we have not defined an x-axis, take the scan axis from the first file
            if not self.axes:
                self.axes.append(d.dtype.names[1])
                setattr(self, self.axes[0], np.array([]))
                self.xlabel = self.axes[0]
            
            for ax in self.axes:
                setattr(self, ax, np.concatenate((getattr(self, ax), d[ax])))


    def load(self, filenames, fileformat='ILL'):
        """ load data and assign signal, monitor and any added axes """
        if fileformat == 'ILL':
            self.load_ill(filenames)

        self.update()

    def hist(self, bins='auto'):
        """ bin the data using the numpy histogram function
            bins are either an integer definining the number of bins or a list of edges
            if bins == 'auto', use self.get_bin_edges()
            bin along axis self.xlabel
            all monitors in self.monitors are binned in the same way """
        
        x = getattr(self, self.xlabel)

        if isinstance(bins, str):
            if bins == 'auto':
                self.calc_bin_edges()
                bins = self.bin_edges

        for mon in self.monitors:
            s, b = np.histogram(x, weights=getattr(self, mon), bins=bins)
            setattr(self, mon, s)

        self.bin_edges = b
        self.bin_centers = np.diff(self.bin_edges)/2 + self.bin_edges[:-1]
        self.binned = True
        
        setattr(self, self.xlabel, self.bin_centers)

        self.update()

    def calc_bin_edges(self):
        """ finds the bin edges assuming discrete measurements """
        diff = np.diff(self.x)
        self.bin_edges = self.x[:-1] + diff/2
        self.bin_edges = np.insert(self.bin_edges, 0, self.x[0] - diff[0]/2)
        self.bin_edges = np.append(self.bin_edges, self.x[-1] + diff[-1]/2)

        self.bin_edges = np.sort(self.bin_edges)

    def bin(self, bins, axis=None):
        """ Bins the data along a certain axis with a specified binsize.
            Usefull for combining similar datasets
            x-values becomes averages rather than bin centers
        """
        # use the scan-axis if not specified
        if axis is None:
            axis = self.xlabel

        xmin = min(getattr(self, axis))
        xmax = max(getattr(self, axis))
        bins = np.arange(xmin-bins/2, xmax+bins, bins)       
        inds = np.digitize(getattr(self, axis), bins)

        # make lists of temporary axes (x) and measurements (y) with each element being initialized as an array
        x = []
        y = []
        for ax in self.axes:
            x.append(np.array([]))

        for meas in self.monitors:
            y.append(np.array([]))
        
        # fill up the temporary axes and measurments
        for i in np.unique(inds):
            # for axes we take the mean value in each bin
            for j, ax in enumerate(self.axes):
                x[j] = np.append(x[j], getattr(self, ax)[inds == i].mean())

            # for measurements we add up the values
            for k, meas in enumerate(self.monitors):
                y[k] = np.append(y[k], getattr(self, meas)[inds == i].sum())

        # update axes and mreasurements
        for i, ax in enumerate(self.axes):
            setattr(self, ax, x[i])

        for j, meas in enumerate(self.monitors):
            setattr(self, meas, y[j])
       
        self.update()

    def plot(self):
        plt.errorbar(self.x, self.I, self.err, fmt='o')

    def plot_bins(self):
        for b in self.bin_edges:
            plt.axvline(x=b, color='k')            

    def __sub__(self, other):
        """ subtraction given that one object is binned 
            returns an object with x, I, err rather than counts and monitor
            err is added in quadrature """
        d = tasdata()
        if self.binned:
            d.x = self.x            
            other.hist(self.bin_edges)            
            d.I = self.I - other.I
            d.err = np.sqrt(self.err**2 + other.err**2)
        elif other.binned:
            d.x = other.x            
            self.hist(other.bin_edges)            
            d.I = self.I - other.I
            d.err = np.sqrt(self.err**2 + other.err**2)

        return d

    def __add__(self, other):
        """ addition given that one object is binned """
        d = tasdata()
        if self.binned:
            d.x = self.x            
            other.hist(self.bin_edges)            
            d.signal = other.signal + self.signal
            d.monitor = other.monitor + self.monitor
        elif other.binned:
            d.x = other.x            
            self.hist(other.bin_edges)            
            d.signal = other.signal + self.signal
            d.monitor = other.monitor + self.monitor

        d.update()
        return d
