from __future__ import print_function
import numpy as np
from scipy.special import sph_harm
from pyglib.math.coordtrans import cart2sph


with open('param.inp', 'r') as f:
    lines = f.readlines()
    lmax = int(lines[0].split()[0])
    v = map(float, lines[1].split()[:3])

r, theta, phi = cart2sph(v)

print(r, theta, phi)

with open('g_ylm.txt', 'w') as f:
    f.write('SINTH,COSTH,SINPH,COSPH:\n')
    f.write(('{:12.7f}'*4+'\n').format(np.sin(theta), np.cos(theta),
            np.sin(phi), np.cos(phi)))
    for l in range(lmax+1):
        for m in range(-l,l+1):
            ylm = sph_harm(m, l, phi, theta)
            f.write('{:2d}{:3d}{:16.10f}{:16.10f}\n'.format(
                    l, m, ylm.real, ylm.imag))

# specific examples
print(sph_harm(-1,1,0.,np.pi/2))
print(sph_harm(1,1,0.,np.pi/2))
