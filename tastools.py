""" tastools module """

import numpy as np
import pylab as plt
from numpy.lib.recfunctions import append_fields
from scipy.special import wofz

def load_ill_ascii(filename):
    """ 
    Loads ILL ASCII data
    Input: filename (ASCII)
    Output: Structured numpy array indexed by column names
    """
    fh = open(filename)
    with fh as f:
        for linenum, line in enumerate(f):
            if 'DATA_:' in line:
                sh = linenum+1
    fh.close()

    return np.genfromtxt(filename, skip_header=sh, names=True)

def load_ill_data(filenumbers, prefix, monitor='M1'):
    """
    Loads one or several ILL data files and returns a single structured array
    """
    if type(filenumbers) is int:
        d = load_ill_ascii(prefix + str(filenumbers))
    else:
        d = load_ill_ascii(prefix + str(filenumbers[0]))
        for f in filenumbers[1::]:
            d = np.append(d, load_ill_ascii(prefix + str(f)))
    
    I = d['CNTS']/d['M1']
    err = np.sqrt(d['CNTS'])/d['M1']

    d = append_fields(d, ['I', 'err'], [I, err])
    return d

def get_edges(grid):
    """
    Get the bin edges of regular 1d bins
    Input: [xmin, xmax, dx] list defining the 1d bins
    Output: list of bind edges
    """
    return np.arange(grid[0]-grid[2]/2, grid[1]+grid[2], grid[2])

def colorplot_2d(x, y, I, xgrid=None, ygrid=None, x_edges=np.array([]), y_edges=np.array([])):
    """
    Create 2d mesh dataset for a colorplot
    Input:
        x: list of x-values
        y: list of y-values
        I: list of intensities
        xgrid: 1d bin centers [xmin, xmax, dx]
        ygrid: 1d bin centers [ymin, ymax, dy]
        x_edges: manually specified edges of x-bins. Ignores xgrid and ygrid if set
        y_edges: manually specified edges of y-bins. Ignores xgrid and ygrid if set
        NOTE: Intensities are calculated as a running average, so somewhat approximated.
    """
    # calculate the edges for our binning unless specified
    if not (x_edges.any() or y_edges.any()):
        x_edges = get_edges(xgrid)
        y_edges = get_edges(ygrid)

    # generate an empty (nan) array based on our binning
    R = np.ones((len(x_edges), len(y_edges)))*np.nan

    # fill the array
    for i in range(len(x)):
        x_index, y_index = np.digitize(x[i], x_edges) - 1, np.digitize(y[i], y_edges) - 1
        if np.isnan(R[x_index, y_index]):
            R[x_index, y_index] = I[i]
        else:
            R[x_index, y_index] = R[x_index, y_index]/2 + I[i]/2

    # generate meshes for X and Y
    X, Y = np.meshgrid(x_edges, y_edges)

    return X, Y, R

def sine_dispersion(qrange, E0=70, width=15, fmt='k--', xgrid=100):
    """
    Simple sine dispersion plot of LO phonon with maximal energy at the zone center.
    hrange: min and max value of q (rlu)
    E0: Minimum energy
    width: width of the dispersion, meaning that Emax = E0 + width
    xgrid: number of points to evaluate
    RETURNS: dispersion as q, E points
    """
    x = np.linspace(qrange[0], qrange[1], xgrid)
    y = width/2*np.cos(x*2*np.pi) + E0 + width/2
    return x, y

def gauss(x, x0, amp, sig):
    norm = 1/sig/np.sqrt(2*np.pi)
    return amp*norm*np.exp(-0.5*((x-x0)/sig)**2)

def lorz(x, x0, amp, gam):
    norm = 1/np.pi/gam
    return amp*norm*gam**2/((x-x0)**2 + gam**2)

def DHO(x, x0, amp, gam):
    return lorz(x, x0, amp, gam)/gam

def voigt(x, x0, amp, sig, gam):
    alpha = sig/np.sqrt(2 * np.log(2))
    return amp*np.real(wofz(((x-x0) + 1j*gam)/alpha/np.sqrt(2)))/alpha/np.sqrt(2*np.pi)
