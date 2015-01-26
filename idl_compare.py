import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u
from astropy.time import Time
from astropy.io import ascii
from astropy.coordinates import *

ltab = ascii.read('radec_fk5.dat', names=['ra', 'dec'])
scl = SkyCoord(ltab['ra']*u.deg, ltab['dec']*u.deg)

obstime = Time('2010-1-1 0:00')
loc = EarthLocation(lat=48*u.deg, lon=-95*u.deg)
scaa = scl.transform_to(AltAz(obstime=obstime, location=loc))

deq2hor = np.loadtxt('altaz_eq2hor.dat').T
sceq2hor = SkyCoord(deq2hor[1]*u.deg, deq2hor[0]*u.deg, frame=scaa.frame)
dhor2eq = np.loadtxt('radec_hor2eq.dat').T
schor2eq = SkyCoord(dhor2eq[0]*u.deg, dhor2eq[1]*u.deg, frame='icrs')

sep1 = scaa.separation(sceq2hor)
plt.figure()
plt.axes(projection='hammer')
plt.title('FK5 RA/Dec')
plt.scatter(Angle(scl.ra.to(u.radian)).wrap_at(180*u.deg), scl.dec.to(u.radian),
            c=sep1.arcsec)
plt.colorbar().set_label('arcsec between astropy FK5->AA and eq2hor')

plt.figure()
plt.axes(projection='hammer')
plt.title('Astropy AA')
plt.scatter(Angle(scaa.az.to(u.radian)).wrap_at(180*u.deg), scaa.alt.to(u.radian),
            c=sep1.arcsec)
plt.colorbar().set_label('arcsec between astropy FK5->AA and eq2hor')


sep2 = scl.separation(schor2eq)
plt.figure()
plt.axes(projection='hammer')
plt.title('FK5 RA/Dec input')
plt.scatter(Angle(scl.ra.to(u.radian)).wrap_at(180*u.deg), scl.dec.to(u.radian),
            c=sep2.arcsec)
plt.colorbar().set_label('arcsec between astropy FK5 and hor2eq')

plt.figure()
plt.axes(projection='hammer')
plt.title('Astropy AA')
plt.scatter(Angle(scaa.az.to(u.radian)).wrap_at(180*u.deg), scaa.alt.to(u.radian),
            c=sep2.arcsec)
plt.colorbar().set_label('arcsec between astropy FK5 and hor2eq')



plt.show()
