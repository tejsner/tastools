import scipy.constants as const
from scipy.constants import milli, angstrom
import numpy as np
import pylab as plt
from numpy.linalg import norm

kb = const.value('Boltzmann constant')
m_e = const.value('neutron mass')
h = const.value('Planck constant')
hbar = const.value('Planck constant over 2 pi')

def gauss(x, sig, x0=0, amp=1):
    norm = (np.sqrt(2.*np.pi) * sig)
    return amp * np.exp(-0.5*((x-x0)/sig)**2.) / norm

def neutron_convert(x, unit='E'): 
    if unit == 'E':
        E = const.eV*milli*x
    elif unit == 'T':
        E = kb*x
    elif unit == 'v':
        E = 0.5*m_e*x**2
    elif unit == 'wl':
        E = h**2/2/m_e/(x*angstrom)**2
    elif unit == 'k':
        E = hbar**2*(x/angstrom)**2/2/m_e
        
    res = dict([])
    res['E'] = E/const.eV/milli
    res['T'] = E/kb
    res['v'] = np.sqrt(2*E/m_e)
    res['wl'] = np.sqrt(h**2/2/m_e/E)/angstrom
    res['k'] = np.sqrt(E*2*m_e/hbar**2)*angstrom
    
    return res

# get scattering triangle vectors in rlu. only in 2d systems with orthogonal axes right now.
def get_triangle(kf=2.661, dE=78.664, h=4.7, k=4.7, a=5.34, b=5.356):
    # get |ki|
    Ei = neutron_convert(kf, 'k')['E'] + dE
    ki = neutron_convert(Ei, 'E')['k']
    
    # get |q|
    a_star = 2*np.pi/a
    b_star = 2*np.pi/b
    qh = h*a_star
    qk = k*b_star
    q = np.sqrt(qh**2 + qk**2)

    # place q along x
    Q = np.array([q, 0])

    # find ki as a vector
    ki_theta = np.arccos((q**2 + ki**2 - kf**2)/(2*q*ki))
    KI = np.array([np.cos(ki_theta), np.sin(ki_theta)])*ki

    # rotate q and ki
    theta = np.arccos(qh/q)
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    Q = np.dot(R, Q)
    KI = np.dot(R, KI)

    # scale to rlu
    Q[0] = Q[0]/a_star
    Q[1] = Q[1]/b_star
    KI[0] = KI[0]/a_star
    KI[1] = KI[1]/b_star

    # trivially find KF
    KF = Q - KI

    return KI, KF, Q

# extend a triangle given by ki and q (in rlu) to kf=ki (neutron loss) in order to look for spurions
def extend_triangle(KI, Q):
    KF = Q - KI

    kf = np.sqrt(sum(KF**2))
    ki = np.sqrt(sum(KI**2))
    KF = (ki/kf)*KF

    Q = KF + KI

    return KI, KF, Q

def shrink_triangle(KI, Q):
    KF = Q - KI

    kf = np.sqrt(sum(KF**2))
    ki = np.sqrt(sum(KI**2))
    KI = (kf/ki)*KI

    Q = KF + KI

    return KI, KF, Q

# plot a scattering triangle
def plot_triangle(KI, KF, Q, fmt='-'):
    plt.plot([0, Q[0]], [0, Q[1]], 'r' + fmt, label='q')
    plt.plot([0, KI[0]], [0, KI[1]], 'b' + fmt, label='ki')
    plt.plot([KI[0], KI[0]+KF[0]], [KI[1], KI[1]+KF[1]], 'k' + fmt, label='kf')

# plot allowed Bmab reflections
def plot_bmab(hmin, hmax, kmin, kmax, fmt='g.'):
    for h in range(hmin, hmax+1):
        for k in range(kmin, kmax+1):
            if h % 2 == 0 and k % 2 == 0:
                plt.plot(h, k, fmt)

# check configuration
def check_config(kf=2.661, dE=78.664, h=4.5, k=4.5, a=5.34, b=5.356, hmin=0, hmax=10, kmin=0, kmax=8, mode='A'):
    KI, KF, Q = get_triangle(kf=kf, dE=dE, h=h, k=k, a=a, b=b)
    plot_triangle(KI, KF, Q, fmt='-')
    
    if mode == 'A':
        KI, KF, Q = extend_triangle(KI, Q)
    elif mode == 'M':
        KI, KF, Q = shrink_triangle(KI, Q)

    plot_triangle(KI, KF, Q, fmt='--')
    plot_bmab(hmin, hmax, kmin, kmax)
    plt.title('({} {} 0) hw={}, kf={}'.format(h, k, dE, kf))
    plt.legend(loc=1)
    plt.xlabel('h [rlu]')
    plt.ylabel('k [rlu]')

# vizualize scan in q
def check_scan(scan, kf=2.661, dE=78.664, a=5.34, b=5.356, mode='A'):
    Q_x = np.array([])
    Q_y = np.array([])
    for hkE in scan:
        KI, KF, Q = get_triangle(kf=kf, h=hkE[0], k=hkE[1], dE=hkE[2], a=a, b=b)
        
        if mode == 'A':
            KI, KF, Q = extend_triangle(KI, Q)
        elif mode == 'M':
            KI, KF, Q = shrink_triangle(KI, Q)
        
        Q_x = np.append(Q_x, Q[0])
        Q_y = np.append(Q_y, Q[1])

    plt.plot(Q_x, Q_y, '.')
    plt.xlabel('h [rlu]')
    plt.ylabel('k [rlu]')

# get intensity of a particular spurion at h, k
def get_spur_I(h, k, dE=78.664, kf=2.661, a=5.34, b=5.356, spur_Q=np.array([8, 4]), sigma=0.1, mode='A'):
    KI, KF, Q = get_triangle(kf=kf, h=h, k=k, dE=dE, a=a, b=b)
    
    if mode == 'A':
        KI, KF, Q = extend_triangle(KI, Q)
    elif mode == 'M':
        KI, KF, Q = shrink_triangle(KI, Q)
    
    dist = norm(Q - spur_Q)
    return gauss(dist, sigma)

# simulate spurious signal on hk, E axes
def vis_spur(hmin, hmax, dh, Emin, Emax, dE, kf=2.661, a=5.34, b=5.356, spur_Q=np.array([8, 4]), sigma=0.1, cmap='viridis', mode='A'):
    gh, gE = np.mgrid[slice(hmin, hmax + dh, dh), slice(Emin, Emax + dE, dE)]
    gI = np.array([[get_spur_I(h, h, dE=E, kf=kf, a=a, b=b, mode=mode, spur_Q=spur_Q, sigma=sigma) for h in gh[:,0]]  for E in gE[0,:]]).T
    gI = gI[:-1, :-1]
    plt.pcolor(gh, gE, gI, cmap=cmap)

# check_config(kf=4.1, dE=78.664, h=5.0, k=5.0, a=5.34, b=5.356, hmin=0, hmax=10, kmin=0, kmax=6, mode='M')
# plt.show()
