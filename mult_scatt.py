""" example of how to calculate a multiple scattering candidate """

from numpy import *
from numpy.linalg import norm
from itertools import product
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D

degrees = 360/2/pi
null_vector = array([0, 0, 0])

def draw_vector(x0, x1, fmt='-', dim=2):
    plt.plot([x0[0], x1[0]], [x0[1], x1[1]], [x0[2], x1[2]], fmt)


a, b, c = 5.376, 5.376, 13.2
a_star, b_star, c_star = 2*pi/a, 2*pi/b, 2*pi/c
h, k, l = 0, 1, 0
hp, kp, lp = -2, 1, 2
gamma = pi/2

ki = 1.48

Q = array([h*a_star, k*b_star, l*c_star])
Q_prime = array([hp*a_star, kp*b_star, lp*c_star])

theta = arcsin(norm(Q)/2/ki)
phi = (pi - 2*theta)/2

KI = - array([ki*cos(phi+gamma), ki*sin(phi+gamma), 0])

KF_prime = Q_prime + KI

print(norm(KF_prime), norm(KI))


fig = plt.figure(dpi=150, figsize=(7,7))
ax = fig.add_subplot(111, projection='3d')

draw_vector(null_vector, -KI, fmt='b-')
draw_vector(-KI, Q, fmt='b-')
draw_vector(null_vector, Q, fmt='b-')

draw_vector(-KI, -KI+KF_prime, fmt='g-')
draw_vector(null_vector, Q_prime, fmt='g-')
draw_vector(Q_prime, Q, fmt='g-')

ax.set_xlim(-2.5, 0.5)
ax.set_xlabel('q_x [AA^(-1)]')
ax.set_ylim(-1, 2)
ax.set_ylabel('q_y [AA^(-1)]')
ax.set_zlim(0, 3)
ax.set_zlabel('q_z [AA^(-1)]')

plt.show()