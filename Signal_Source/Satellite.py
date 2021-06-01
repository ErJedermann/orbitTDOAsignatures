from Signal_Source.SignalSourceInterface import SignalSource
from TLEcalculator import TLEcalculator
from astropy.time import Time
import astropy.coordinates as coord
import ITRSconverter
import random
import numpy as np
import astropy.units as u

class SatelliteSource(SignalSource):
    def __init__(self, auth: bool, tle_file: str = None, verbose: bool = True):
        self.auth = auth
        self.tle_file = tle_file
        self.verbose = verbose
        self.tle_calc = TLEcalculator(tle_file, verbose=False)
        self.receiving_angle = 80 # 80° from the zenith are receivable
        self.speed_of_light = 299792458 # meter/sec
        self.one_sec_jDay = 1 / 86400  # one second in julian day
        self.deviation = 0  # deviation of the signal source from the data-base position in km

    def __is_satellite_visible(self, rec_pos: [(float, float, float)], sat_pos: (float, float, float)) -> bool:
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


    def __get_satellite_positions_TEME(self, sat_index: int, measurement_times: [Time]):
        sat_pos_TEME = []
        target_sat = self.tle_calc.satList[sat_index]
        for i in range(len(measurement_times)):
            jDay = measurement_times[i].jd1
            jDayF = measurement_times[i].jd2
            err, pos, vel = self.tle_calc.calculate_position_single(target_sat, jDay, jDayF)
            pos = np.array(pos) * 1000  # km in m
            sat_pos_TEME.append(pos)
        return sat_pos_TEME

    def __get_deviation(self) -> [float]:
        # returns the deviation of the satellite in meter
        if self.deviation is None:
            return [0,0,0]
        else:
            vector = np.random.rand(3)
            norm_vector = vector / np.linalg.norm(vector)
            return norm_vector * np.random.random(1) * self.deviation * 1000

    def get_positions_TEME(self, receivers: [coord.ITRS], measurement_times: [Time]) -> (bool, int, [(float, float, float)]):
        sat_amount = len(self.tle_calc.satList)
        invisible_counter = 0
        invisible_border = sat_amount
        start_rec_pos_TEME = self.__rec_ITRS_2_TEME(receivers, measurement_times[0].jd1, measurement_times[0].jd2)
        end_rec_pos_TEME = self.__rec_ITRS_2_TEME(receivers, measurement_times[-1].jd1, measurement_times[-1].jd2)
        while True:
            temp_i = int(random.random() * sat_amount)
            if invisible_counter >= invisible_border:
                # start searching systematic
                temp_i = invisible_counter - invisible_border
                if temp_i == sat_amount:
                    print(f"ERROR: Signal_Source.Satellite: No visible satellite! Consider selecting other time-point.")
                    return False, None, None
            temp_sat = self.tle_calc.satList[temp_i]
            temp_err = np.array(self.__get_deviation())
            # pos & vel in km in TEME
            err, pos, vel = self.tle_calc.calculate_position_single(temp_sat, measurement_times[0].jd1,
                                                                    measurement_times[0].jd2)
            err2, pos2, vel2 = self.tle_calc.calculate_position_single(temp_sat, measurement_times[-1].jd1,
                                                                       measurement_times[-1].jd2)
            if err == 0 & err2 == 0:
                sat_start_pos = np.array(pos+temp_err) * 1000  # km in m
                sat_end_pos = np.array(pos2+temp_err) * 1000
                if self.__is_satellite_visible(start_rec_pos_TEME, sat_start_pos+temp_err) & self.__is_satellite_visible(
                        end_rec_pos_TEME, sat_end_pos+temp_err):
                    sat_pos_list = self.__get_satellite_positions_TEME(temp_i, measurement_times)
                    sat_pos_list = sat_pos_list + temp_err
                    return True, temp_i, sat_pos_list
                else:
                    invisible_counter += 1
            else:
                invisible_counter += 1

    def set_deviation(self, range: float):
        # position deviation of the signal source from the real data-base position in km
        self.deviation = range

    def is_authentic(self):
        return self.auth

    def get_names(self):
        return self.tle_calc.sat_names
