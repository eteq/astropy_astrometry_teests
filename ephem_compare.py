import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u
from astropy.time import Time
from astropy.io import ascii
from astropy.coordinates import *

import ephem

ltab = ascii.read('radec_fk5.dat', names=['ra', 'dec'])
scl = SkyCoord(ltab['ra']*u.deg, ltab['dec']*u.deg)

obstime = Time('2010-1-1 0:00')
loc = EarthLocation(lat=48*u.deg, lon=-95*u.deg)
scaa = scl.transform_to(AltAz(obstime=obstime, location=loc))


ephemalt = []
ephemaz = []
eobserver = ephem.Observer()
eobserver.lon = '-95:00:00'
eobserver.lat = '48:00:00'
eobserver.elevation = 0
eobserver.date = '2010/1/1 0:00'
eobserver.pressure = 0
for ra, dec in zip(scl.ra, scl.dec):
    f = ephem.FixedBody()
    f.name = 'testbody'
    f._ra = ra.radian
    f._dec = dec.radian
    f.compute(eobserver)
    ephemalt.append(float(f.alt))
    ephemaz.append(float(f.az))

scephem = SkyCoord(alt=ephemalt*u.radian, az=ephemaz*u.radian, frame=scaa.frame)
with open('altaz_ephem.dat', 'w') as f:
    for altdeg, azdeg in zip(scephem.alt.degree, scephem.az.degree):
        f.write('{0} {1}\n'.format(altdeg, azdeg))

sep1 = scaa.separation(scephem)
plt.figure()
plt.axes(projection='hammer')
plt.title('FK5 RA/Dec')
plt.scatter(Angle(scl.ra.to(u.radian)).wrap_at(180*u.deg), scl.dec.to(u.radian),
            c=sep1.arcsec)
plt.colorbar().set_label('arcsec between astropy FK5->AA and pyephem')

plt.show()
