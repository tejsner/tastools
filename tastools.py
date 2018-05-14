""" tastools module """

import numpy as np
import pylab as plt

def load_ill_data(filename):
    """ 
    Loads ILL data
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

def get_ill_data(filenumbers, prefix='data/rawdata/0', axis='EN'):
    """
    Loads ILL data file(s) along a certain axis and return x, y, err
    Input: filenumber(s) as list or int
           prefix: path where the files are stored
           axis: x-axis to evaluate
    Output: x, I, err: 3 numpy arrays
    """
    x = np.array([])
    I = np.array([])
    err = np.array([])

    if type(filenumbers) is int:
        filenumbers = [filenumbers]

    for f in filenumbers:
        d = load_ill_data(prefix + str(f))
        x = np.append(x, d[axis])
        I = np.append(I, d['CNTS']/d['M1'])
        err = np.append(err, np.sqrt(d['CNTS'])/d['M1'])

    return x, I, err

def get_ill_data_2d(filenumbers, prefix='data/rawdata/0', axis1='EN', axis2='QH'):
    """
    Loads ILL data files along two seperate axes
    Input: filenumbers as list
           prefix: path
           axis1: first axis
           axis2: second axis
    Output:: axis1, axis2, I: 3 numpy arrays
    """
    ax1 = np.array([])
    ax2 = np.array([])
    I = np.array([])

    for f in filenumbers:
        d = load_ill_data(prefix + str(f))
        ax1 = np.append(ax1, d[axis1])
        ax2 = np.append(ax2, d[axis2])
        I = np.append(I, d['CNTS'])

    return ax1, ax2, I
