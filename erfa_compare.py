import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u
from astropy.time import Time
from astropy.io import ascii
from astropy.coordinates import *

ltab = ascii.read('radec_icrs.dat', names=['ra', 'dec'])
scl = SkyCoord(ltab['ra']*u.deg, ltab['dec']*u.deg)

obstime = Time('2010-1-1 0:00')
loc = EarthLocation(lat=48*u.deg, lon=-95*u.deg)
scaa = scl.transform_to(AltAz(obstime=obstime, location=loc))

derfa = np.loadtxt('altaz_erfa.dat').T
scerfa = SkyCoord(derfa[1]*u.deg, derfa[0]*u.deg, frame=scaa.frame)

sep1 = scaa.separation(scerfa)
plt.figure()
plt.axes(projection='hammer')
plt.title('ICRS RA/Dec')
plt.scatter(Angle(scl.ra.to(u.radian)).wrap_at(180*u.deg), scl.dec.to(u.radian),
            c=sep1.arcsec)
plt.colorbar().set_label('arcsec between astropy ICRS->AA and eraAtco13')

plt.show()
