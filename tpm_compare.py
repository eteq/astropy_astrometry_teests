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

# the value of toffset comes from using find_optimal_offset - the default value
# is what's needed to remove the effect
def compute_pytpm_aa(scl, toffset=-1.368*u.second):
    v6cat = []
    for rarad, decrad in zip(scl.ra.radian, scl.dec.radian):
        v6cat.extend(pytpm.convert.cat2v6(rarad, decrad))
    v6cat = pytpm.convert.proper_motion(v6cat, tt, pytpm.tpm.J2000)
    v6_app = pytpm.convert.convertv6(v6cat, s1=6, s2=18, utc=obstime.jd + toffset.to(u.day).value,
                                     lon=loc.longitude.value, lat=loc.latitude.value,
                                     alt=loc.height.value,
                                     P=1, xpole=pm_x_rad, ypole=pm_y_rad)

    dlist = pytpm.convert.v62cat(v6_app)
    tpmazrad = []
    tpmzenrad = []
    for d in dlist:
        tpmazrad.append(d['alpha'])
        tpmzenrad.append(d['delta'])
    return SkyCoord(alt=tpmzenrad*u.radian, az=tpmazrad*u.radian, frame=scaa.frame)


def find_optimal_offset(tlower=-3, tupper=0, n=250):
    """
    The astropy and pytpm timescales seem to be offset by a fixed amount.  This
    computes the total offset for a variety of different offset times and
    plots them.  I used this to set the default for compute_pytpm_aa (and
    therefore main)
    """
    toffsets = np.linspace(tlower, tupper, n)*u.second
    avgsep = []
    for i, toffset in enumerate(toffsets):
        print 'Computing offset for #', i, 'of', len(toffsets)
        sctpm = compute_pytpm_aa(scl, toffset)
        sep = scaa.separation(sctpm)
        avgsep.append(np.mean(sep))

    avgsep = u.Quantity(avgsep).to(u.arcsec)
    plt.plot(toffsets, avgsep,'-o')
    plt.xlabel('time offset [sec]')
    plt.ylabel('avg separation [arcsec]')

    plt.show()


def main():
    toffset = 1*u.second
    sctpm = compute_pytpm_aa(scl)
    with open('altaz_tpm.dat', 'w') as f:
        for altdeg, azdeg in zip(sctpm.alt.degree, sctpm.az.degree):
            f.write('{0} {1}\n'.format(altdeg, azdeg))

    sep = scaa.separation(sctpm)

    plt.figure()
    plt.axes(projection='hammer')
    plt.title('FK5 RA/Dec')
    plt.scatter(Angle(scl.ra.to(u.radian)).wrap_at(180*u.deg), scl.dec.to(u.radian),
                c=sep.arcsec)
    plt.colorbar().set_label('arcsec between astropy FK5->AA and pytpm')

    plt.figure()
    plt.axes(projection='hammer')
    plt.title('TPM altaz')
    plt.scatter(Angle(sctpm.az.to(u.radian)).wrap_at(180*u.deg), sctpm.alt.to(u.radian),
                c=sep.arcsec)
    plt.colorbar().set_label('arcsec between astropy FK5->AA and pytpm')

    plt.show()

find_optimal_offset()
#main()
