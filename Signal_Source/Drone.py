from Signal_Source.SignalSourceInterface import SignalSource
from TLEcalculator import TLEcalculator
from astropy.time import Time
import astropy.coordinates as coord
import ITRSconverter
import random
import numpy as np
import astropy.units as u

class DroneSource(SignalSource):
    def __init__(self, auth: bool, height_min: float, height_max: float, max_speed: float = 0, above_receiver: int = -1):
        self.receiving_angle = 80 # 80° from the zenith are receivable
        self.speed_of_light = 299792458 # meter/sec
        self.one_sec_jDay = 1 / 86400  # one second in julian day
        self.auth = auth
        self.height_min = height_min  # height in m
        self.height_max = height_max
        self.drone_speed = max_speed  # speed in m/s
        self.above_receiver = above_receiver  # -1: (default) the centroid should be used, 0-n: above receiver 0-n

    def __is_sender_visible(self, rec_pos: [(float, float, float)], sat_pos: (float, float, float)) -> bool:
        sat_pos = np.array(sat_pos)
        for i in range(len(rec_pos)):
            temp_position = np.array(rec_pos[i])
            vector_rec_sat = sat_pos - temp_position
            angle_between = np.dot(temp_position, vector_rec_sat) / (
                    np.sqrt(sum(temp_position ** 2)) * np.sqrt(np.sum(np.power(vector_rec_sat, 2))))
            angle_between = np.arccos(angle_between)
            angle_between = np.degrees(angle_between)
            if angle_between > self.receiving_angle:
                return False
        return True


    def __get_pos_in_m(self, input: coord.ITRS):
        if input.cartesian.xyz.unit is u.km:
            return input.cartesian.xyz.value * 1000
        else:
            return input.cartesian.xyz.value

    def __rec_ITRS_2_TEME(self, receivers: [coord.ITRS], jDay: float, jDayF: float) -> [(float, float, float)]:
        pos_list_TEME = []
        for i in range(len(receivers)):
            pos_ITRS = self.__get_pos_in_m(receivers[i])
            pos_TEME, vel_TEME = ITRSconverter.ITRS_2_TEME(pos_ITRS, jDay, jDayF)
            pos_list_TEME.append(pos_TEME)
        return pos_list_TEME

    def __calc_centroid(self, pos_list: [(float, float, float)]):
        pos_list = np.array(pos_list)  # shape(3, rec)
        pos_sum = np.sum(pos_list, axis=0)  # shape(3)
        return pos_sum / len(pos_list)

    def __find_min_height(self, cen_lon: float, cen_lat: float, cen_h: float, rec_ITRS: [(float, float, float)]):
        # increase until it is visible
        temp_height = cen_h  # the centroid should not be visible
        temp_pos = ITRSconverter.LatLonHeight_2_ITRS(cen_lat, cen_lon, temp_height)
        temp_pos = self.__get_pos_in_m(temp_pos)
        while not self.__is_sender_visible(rec_ITRS, temp_pos):
            temp_height = temp_height * 2
            temp_pos = ITRSconverter.LatLonHeight_2_ITRS(cen_lat, cen_lon, temp_height)
            temp_pos = self.__get_pos_in_m(temp_pos)
        # now the visibility-border is between temp_height and temp_height/2
        while self.__is_sender_visible(rec_ITRS, temp_pos):
            temp_height = temp_height * 0.95  # this will decrease in max 14 steps below temp_height/2
            temp_pos = ITRSconverter.LatLonHeight_2_ITRS(cen_lat, cen_lon, temp_height)
            temp_pos = self.__get_pos_in_m(temp_pos)
        # now it is not visible, so remove the last operation
        temp_height = temp_height / 0.95
        return temp_height
    
    def __east_north_up_shift_2_xyz_shift(self, start_pos: coord.ITRS, north_shift: float, east_shift: float, up_shift):
        # north_shift and east_shift in [m]
        # from https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
        phi = start_pos.earth_location.geodetic.lat.value
        lam = start_pos.earth_location.geodetic.lon.value
        delta_x = - np.sin(np.radians(lam)) * east_shift \
                  - np.cos(np.radians(lam)) * np.sin(np.radians(phi)) * north_shift \
                  + np.cos(np.radians(lam)) * np.cos(np.radians(phi)) * up_shift
        delta_y = + np.cos(np.radians(lam)) * east_shift \
                  - np.sin(np.radians(lam)) * np.sin(np.radians(phi)) * north_shift \
                  + np.sin(np.radians(lam)) * np.cos(np.radians(phi)) * up_shift
        delta_z = 0 + np.cos(np.radians(phi)) * north_shift + np.sin(np.radians(phi)) * up_shift
        return delta_x, delta_y, delta_z # in [m]

    def get_positions_TEME(self, receivers: [coord.ITRS], measurement_times: [Time]) -> (int, [(float, float, float)]):
        rec_start_ITRS = []
        for i in range(len(receivers)):
            rec_start_ITRS.append(self.__get_pos_in_m(receivers[i]))
        if self.above_receiver == -1:
            # get the centroid of the receivers
            centroid_pos_ITRS = self.__calc_centroid(rec_start_ITRS)
            centroid_ITRS = ITRSconverter.Cartesian_2_ITRS(centroid_pos_ITRS / 1000, (0,0,0), measurement_times[0].jd1, measurement_times[0].jd2)
            cen_lon, cen_lat, cen_h = ITRSconverter.ITRS_2_LonLatHeight(centroid_ITRS)
            cen_h = cen_h * 1000  # km in m
        else:
            # get position of receiver n
            n_rec = self.above_receiver % len(receivers)
            cen_lon, cen_lat, cen_h = ITRSconverter.ITRS_2_LonLatHeight(receivers[n_rec])
            cen_h = cen_h * 1000  # km in m
        if cen_h < 0:
            cen_h = 1
        min_height = self.__find_min_height(cen_lon, cen_lat, cen_h, rec_start_ITRS)
        lower_height = max(self.height_min, min_height)
        upper_height = min(self.height_max, min_height)
        range_height = upper_height - lower_height
        random_height = random.random() * range_height + lower_height
        drone_start_ITRS = ITRSconverter.LatLonHeight_2_ITRS(cen_lat, cen_lon, random_height)
        # get a velocity direction
        random_direction = random.random() * 360  # 0° is north, 90° is east
        north_movement = np.cos(np.radians(random_direction)) * self.drone_speed
        east_movement = np.sin(np.radians(random_direction)) * self.drone_speed
        # calculate the TEME positions
        drone_TEME_pos = []
        start_time = measurement_times[0].jd1 + measurement_times[0].jd2
        start_pos_ITRS = self.__get_pos_in_m(drone_start_ITRS)
        for i_time in range(len(measurement_times)):
            temp_jd = measurement_times[i_time]
            temp_time = temp_jd.jd1 + temp_jd.jd2
            temp_time = temp_time - start_time
            delta_ITRS = self.__east_north_up_shift_2_xyz_shift(drone_start_ITRS, north_movement * temp_time, east_movement * temp_time, 0)
            delta_ITRS = np.array(delta_ITRS)
            pos_TEME, vel_TEME = ITRSconverter.ITRS_2_TEME(start_pos_ITRS + delta_ITRS, temp_jd.jd1, temp_jd.jd2)
            drone_TEME_pos.append(pos_TEME)
        index_return = 0
        return True, index_return, drone_TEME_pos

    def set_deviation(self, range: float):
        pass


    def is_authentic(self):
        return self.auth

    def get_names(self):
        return [f'Drone.attack:{self.auth}']
