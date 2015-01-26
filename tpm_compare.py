import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u
from astropy.time import Time
from astropy.io import ascii
from astropy.coordinates import *
from astropy.utils import iers

import pytpm

ltab = ascii.read('radec_fk5.dat', names=['ra', 'dec'])
scl = SkyCoord(ltab['ra']*u.deg, ltab['dec']*u.deg)

obstime = Time('2010-1-1 0:00')
loc = EarthLocation(lat=48*u.deg, lon=-95*u.deg)
scaa = scl.transform_to(AltAz(obstime=obstime, location=loc))

pm_x, pm_y = iers.IERS.open().pm_xy(obstime.jd)
pm_x_rad = pm_x.to(u.radian).value
pm_y_rad = pm_y.to(u.radian).value
tt = obstime.tdb.jd
v6cat = []
for rarad, decrad in zip(scl.ra.radian, scl.dec.radian):
    v6cat.extend(pytpm.convert.cat2v6(rarad, decrad))
v6cat = pytpm.convert.proper_motion(v6cat, tt, pytpm.tpm.J2000)
v6_app = pytpm.convert.convertv6(v6cat, s1=6, s2=19, utc=obstime.jd,
                                 lon=loc.longitude.value, lat=loc.latitude.value,
                                 alt=loc.height.value,
                                 P=1, xpole=pm_x_rad, ypole=pm_y_rad)

dlist = pytpm.convert.v62cat(v6_app)
tpmazrad = []
tpmzenrad = []
for d in dlist:
    tpmazrad.append(d['alpha'])
    tpmzenrad.append(d['delta'])
sctpm = SkyCoord(alt=tpmzenrad*u.radian, az=tpmazrad*u.radian, frame=scaa.frame)
with open('altaz_tpm.dat', 'w') as f:
    for altdeg, azdeg in zip(sctpm.alt.degree, sctpm.az.degree):
        f.write('{0} {1}\n'.format(altdeg, azdeg))

sep1 = scaa.separation(sctpm)
print('avg dalt:', np.mean((scaa.alt - sctpm.alt).to(u.arcsec)),
           'daz:', np.mean(scaa.az - sctpm.az).to(u.arcsec))

plt.figure()
plt.axes(projection='hammer')
plt.title('FK5 RA/Dec')
plt.scatter(Angle(scl.ra.to(u.radian)).wrap_at(180*u.deg), scl.dec.to(u.radian),
            c=sep1.arcsec)
plt.colorbar().set_label('arcsec between astropy FK5->AA and pytpm')

plotlon = sctpm.az
plotlat = sctpm.alt
plottitle = 'TPM altaz'
plt.figure()
plt.axes(projection='hammer')
plt.title(plottitle)
plt.scatter(Angle(plotlon.to(u.radian)).wrap_at(180*u.deg), plotlat.to(u.radian),
            c=sep1.arcsec)
plt.colorbar().set_label('arcsec between astropy FK5->AA and pytpm')


plt.figure()
plt.axes(projection='hammer')
plt.title(plottitle)
plt.scatter(Angle(plotlon.to(u.radian)).wrap_at(180*u.deg), plotlat.to(u.radian),
            c=(scaa.alt - sctpm.alt).to(u.arcsec).value)
plt.colorbar().set_label('dalt')


plt.figure()
plt.axes(projection='hammer')
plt.title(plottitle)
plt.scatter(Angle(plotlon.to(u.radian)).wrap_at(180*u.deg), plotlat.to(u.radian),
            c=(scaa.az - sctpm.az).to(u.arcsec).value)
plt.colorbar().set_label('daz')

plt.show()
