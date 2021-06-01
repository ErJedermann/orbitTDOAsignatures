from astropy import units as u
from astropy import coordinates as coord
from astropy.time import Time
import skyfield.sgp4lib as sgp4lib

import numpy as np

def TEME_2_ITRS(pos: np.ndarray, vel: np.ndarray, jDay:float, jDayF: float) -> coord.ITRS:
    # pos in km and vel in km/s
    # sgp4lib: ~200 sats in 0.1 sec
    pos, vel = sgp4lib.TEME_to_ITRF(jDay + jDayF, np.asarray(pos), np.asarray(vel) * 86400)
    vel = vel / 86400
    time = Time(jDay+jDayF, format='jd')
    itrs = coord.ITRS(pos[0] * u.km, pos[1] * u.km, pos[2] * u.km,
                      vel[0] * u.km / u.s, vel[1] * u.km / u.s, vel[2] * u.km / u.s,
                      obstime=time)
    return itrs

def ITRS_2_LonLatHeight(obj: coord.ITRS) -> (float, float, float):
    location = obj.earth_location
    lon = location.geodetic.lon
    lat = location.geodetic.lat
    height = location.geodetic.height
    # lon in deg, lat in deg, height in km
    return lon.value, lat.value, height.value

def LatLonHeight_2_ITRS(lat: float, lon: float, height: float) -> coord.ITRS:
    # lat is +north-south, lon is +east-west
    # lon in deg, lat in deg, height in m
    location = coord.EarthLocation(lon=lon * u.deg, lat=lat * u.deg, height=height * u.m)
    return location.itrs

def Cartesian_2_ITRS(pos: (float, float, float), vel: (float, float, float),
                     jDay: float=0.0, jDayF: float=0.0) -> coord.ITRS:
    # receives pos in [km] and vel in [km/s]
    if jDay == 0.0 and jDayF == 0.0:
        itrs = coord.ITRS(pos[0] * u.km, pos[1] * u.km, pos[2] * u.km,
                          vel[0] * u.km / u.s, vel[1] * u.km / u.s, vel[2] * u.km / u.s)
    else:
        itrs = coord.ITRS(pos[0] * u.km, pos[1] * u.km, pos[2] * u.km,
                          vel[0]*u.km/u.s, vel[1]*u.km/u.s, vel[2]*u.km/u.s,
                          obstime=Time(jDay+jDayF, format='jd'))
    return itrs

def ITRS_2_TEME(pos_itrs: np.ndarray, jDay: float, jDayF: float = 0.0, vel_itrs: np.ndarray = None, xp: float = 0.0, yp: float = 0.0):
    if vel_itrs is None:
        vel_itrs = np.zeros(np.shape(pos_itrs))
    # go in the other direction than in the sgp4lib: ITRS -> PEF -> TEME
    if xp == 0.0 and yp == 0.0:
        pos_pef = pos_itrs
        vel_pef = vel_itrs
    else:
        W = (sgp4lib.rot_x(-yp)).dot(sgp4lib.rot_y(-xp))
        pos_pef = (W).dot(pos_itrs)
        vel_pef = (W).dot(vel_itrs)
    # now calculate the angle theta
    theta, theta_dot = sgp4lib.theta_GMST1982(jDay, jDayF)
    angular_velocity = sgp4lib.multiply.outer(sgp4lib._zero_zero_minus_one, -theta_dot)  # theta_dot -> -theta_dot
    # invert the angle to rotation back (-theta -> theta)
    R = sgp4lib.rot_z(theta)
    if len(pos_pef.shape) == 1:
        pos_teme = (R).dot(pos_pef)
        vel_pef = vel_pef + sgp4lib._cross(angular_velocity, pos_pef)
        vel_teme = (R).dot(vel_pef)
    else:
        pos_teme = sgp4lib.mxv(R, pos_pef)
        vel_pef = vel_pef + sgp4lib._cross(angular_velocity[:,None], pos_pef)
        vel_teme = sgp4lib.mxv(R, vel_pef)
    return pos_teme, vel_teme








