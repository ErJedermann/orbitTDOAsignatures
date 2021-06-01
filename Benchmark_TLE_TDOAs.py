from astropy import units as u
from astropy import coordinates as coord
from astropy.time import Time
from TLEcalculator import TLEcalculator
from Benchmark_SetupReceiver import BenchmarkSetupReceiver
from Authentication_Algorithms.Authentication_Algorithm_Interface import Authentication_Algorithm
from TimingErrors.TimingErrorInterface import TimingError

import numpy as np
import random
import ITRSconverter


# This class is the core of the benchmark. It can handle single and multiple measurements
class Benchmark_TLE_TDOAs:
    def __init__(self, TLE_file: str):
        self.rec_setup = BenchmarkSetupReceiver()
        self.tle_file = TLE_file
        self.tle_calculator = TLEcalculator(TLE_file, warnings=False)
        self.receiving_angle = 80 # 80° from the zenith are receivable
        self.speed_of_light = 299792458 # meter/sec
        self.one_sec_jDay = 1 / 86400  # one second in julian day


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

    # step1: get some receivers of your desired amount and distribution
    def get_receivers(self, amount:int, min_diameter: float = 1000, max_diameter: float = 1000) -> [coord.ITRS]:
        return self.rec_setup.get_random_receivers(amount, min_diameter, max_diameter)

    # step2: get some measurement-times
    def get_measurement_times(self, amount: int, measurement_delay: float = 1) -> [Time]:
        measurement_times = []
        epoch_min, epoch_max = self.tle_calculator.get_min_max_epoch()
        valid_days = self.tle_calculator.valid_days
        jDay_duration = amount * measurement_delay * self.one_sec_jDay
        jDay_min = self.tle_calculator.epoch_2_jDay(epoch_min[0], epoch_min[1])
        jDay_max = self.tle_calculator.epoch_2_jDay(epoch_max[0], epoch_max[1])
        allowed_range = (jDay_max[0]+jDay_max[1]) - (jDay_min[0]+jDay_min[1]) - jDay_duration + 2*valid_days
        start_time = random.random() * allowed_range - valid_days + (jDay_min[0] + jDay_min[1])
        for i2 in range(amount):
            temp_delay = i2 * measurement_delay * self.one_sec_jDay
            temp_jd = Time(start_time + temp_delay, format='jd')
            measurement_times.append(temp_jd)
        return measurement_times

    # step3: get a visible satellite, that will be used as a sender (alternatively use a prepared list of sender-pos)
    # If no extra TLEcalculator is given, a random satellite from the default TLEcalculator is used.
    def get_visible_satellite(self, receivers: [coord.ITRS], measurement_times: [Time], used_tle_calc: TLEcalculator = None, verbose: bool=False) -> (bool, int, [(float, float, float)]):
        if used_tle_calc is None:
            used_tle_calc = self.tle_calculator
        sat_amount = len(used_tle_calc.satList)
        invisible_counter = 0
        start_rec_pos_TEME = self.__rec_ITRS_2_TEME(receivers, measurement_times[0].jd1, measurement_times[0].jd2)
        end_rec_pos_TEME = self.__rec_ITRS_2_TEME(receivers, measurement_times[-1].jd1, measurement_times[-1].jd2)
        while True:
            if invisible_counter > len(used_tle_calc.satList):
                if verbose:
                    print(f"DEBUG: Ben_TLE_TDOA.get_vis_sat: Many not visible satellites. Select other time-point.")
                return False, None, None
            temp_i = int(random.random() * sat_amount)
            temp_sat = used_tle_calc.satList[temp_i]
            # pos & vel in km in TEME
            err, pos, vel = used_tle_calc.calculate_position_single(temp_sat, measurement_times[0].jd1, measurement_times[0].jd2)
            err2, pos2, vel2 = used_tle_calc.calculate_position_single(temp_sat, measurement_times[-1].jd1, measurement_times[-1].jd2)
            if err == 0 & err2 == 0:
                sat_start_pos = np.array(pos) * 1000  # km in m
                sat_end_pos = np.array(pos2) * 1000
                if self.__is_satellite_visible(start_rec_pos_TEME, sat_start_pos) & self.__is_satellite_visible(end_rec_pos_TEME, sat_end_pos):
                    sat_pos_list = self.__get_satellite_positions_TEME(temp_i, measurement_times, used_tle_calc)
                    return True, temp_i, sat_pos_list
                else:
                    invisible_counter += 1
            else:
                invisible_counter += 1

    # step4: verify that the receivers can see the sender (especially when the sender is an attacker)
    def verify_sender_visibility(self, receivers: [coord.ITRS], measurement_times: [Time], sender_pos_TEME: [(float, float, float)]) -> bool:
        for i in range(len(measurement_times)):
            temp_rec_pos_TEME = self.__rec_ITRS_2_TEME(receivers, measurement_times[i].jd1, measurement_times[i].jd2)
            temp_sender_pos_TEME = sender_pos_TEME[i]
            if not self.__is_satellite_visible(temp_rec_pos_TEME, temp_sender_pos_TEME):
                return False
        return True



    def __get_satellite_positions_TEME(self, sat_index: int, measurement_times: [Time], used_tle_calculator: TLEcalculator):
        sat_pos_TEME = []
        target_sat = used_tle_calculator.satList[sat_index]
        for i in range(len(measurement_times)):
            jDay = measurement_times[i].jd1
            jDayF = measurement_times[i].jd2
            err, pos, vel = used_tle_calculator.calculate_position_single(target_sat, jDay, jDayF)
            pos = np.array(pos) * 1000  # km in m
            sat_pos_TEME.append(pos)
        return sat_pos_TEME


    def __calculate_TOF(self, rec_pos: [(float, float, float)], sender_pos: (float, float, float)):
        sender_pos = np.array(sender_pos)
        TOF_list = []
        for i in range(len(rec_pos)):
            temp_rec = np.array(rec_pos[i])
            temp_dist = np.sqrt(np.sum(np.power(sender_pos - temp_rec, 2)))
            temp_TOF = temp_dist / self.speed_of_light
            TOF_list.append(temp_TOF)
        return TOF_list

    def __calculate_TOF2(self, rec_pos: [(float, float, float)], sender_pos: (float, float, float)):
        sender_pos = np.array(sender_pos)  # shape(3,)
        rec_pos = np.array(rec_pos)  # shape(rec, 3)
        dists = np.sqrt(np.power(rec_pos[:, 0] - sender_pos[0], 2) +
                       np.power(rec_pos[:, 1] - sender_pos[1], 2) +
                       np.power(rec_pos[:, 2] - sender_pos[2], 2))
        TOFs = dists / self.speed_of_light
        return TOFs


    def __calculate_splitted_distance(self, receivers: [coord.ITRS], index1: int, index2: int, jDay: float, jDayF: float):
        err1, pos1, vel1 = self.tle_calculator.calculate_position_single(self.tle_calculator.satList[index1], jDay, jDayF)
        sat1 = ITRSconverter.TEME_2_ITRS(pos1, vel1, jDay, jDayF)
        sat1_pos = sat1.cartesian.xyz.value  # km
        err2, pos2, vel2 = self.tle_calculator.calculate_position_single(self.tle_calculator.satList[index2], jDay, jDayF)
        sat2 = ITRSconverter.TEME_2_ITRS(pos2, vel2, jDay, jDayF)
        sat2_pos = sat2.cartesian.xyz.value  # km
        vector_total_dist = sat1_pos - sat2_pos
        dist_total_scalar = np.sqrt(vector_total_dist.dot(vector_total_dist))
        centroid = self.rec_setup.get_centroid_pos(receivers) # x,y,z in km
        vector_rec_sat = sat1_pos - centroid
        # get the 'vertical' distance by calculating the projection of vector_dist on vector_rec_sat
        two_norm = np.sqrt(sum(vector_rec_sat ** 2))
        dist_vertical = (np.dot(vector_total_dist, vector_rec_sat) / two_norm ** 2) * vector_rec_sat
        dist_vertical_scalar = np.sqrt(dist_vertical.dot(dist_vertical))
        # get the 'horizontal' distance by vector_total_dist - dist_vertical
        dist_horizontal = vector_total_dist - dist_vertical
        dist_horizontal_scalar = np.sqrt(dist_horizontal.dot(dist_horizontal))
        return dist_total_scalar, dist_horizontal_scalar, dist_vertical_scalar

    def __calculate_average_splitted_distances(self, receivers: [coord.ITRS], sat_index1: int, sat_index2: int,
                                      receiving_times:[[float]]) -> (float, float, float):
        dist_tot = 0.0
        dist_hor = 0.0
        dist_vert = 0.0
        for time_index in range(len(receiving_times)):
            temp_time = receiving_times[time_index]
            temp_tot, temp_hor, temp_vert = self.__calculate_splitted_distance(receivers, sat_index1, sat_index2, temp_time[0], temp_time[1])
            dist_tot += temp_tot
            dist_hor += temp_hor
            dist_vert += temp_vert
        return dist_tot/len(receiving_times), dist_hor/len(receiving_times), dist_vert/len(receiving_times)

    # step5: The main step. Simulate the situation and get the reaction of the algorithms.
    def perform_mesaurement(self, receiver_list: [coord.ITRS], auth_algs: [Authentication_Algorithm],
                            sender_TEME_positions: [(float, float, float)], measurement_times: [Time],
                            sender_is_valid: bool, sender_sat_name: str = None, error: TimingError = None,
                            verbose: bool = True, plot_tdof_graph: bool = False) -> ([float], [float], [float], [float], [float]):
        TDOA_measurements = []
        TOA_times = []
        # go through all measurement-times and collect the TDOA-measurements
        for i_time in range(len(measurement_times)):
            temp_time = measurement_times[i_time]
            temp_rec_pos_list = self.__rec_ITRS_2_TEME(receiver_list, temp_time.jd1, temp_time.jd2)
            temp_sender_pos = sender_TEME_positions[i_time]
            temp_TOF = self.__calculate_TOF(temp_rec_pos_list, temp_sender_pos)
            temp_TOF = np.array(temp_TOF)
            if error is not None:
                add_err = error.get_error(len(temp_TOF))
                temp_TOF = temp_TOF + add_err
            temp_TDOF = temp_TOF - temp_TOF[0]
            TDOA_measurements.append(temp_TDOF)
            TOF0_jd = temp_TOF[0] * self.one_sec_jDay
            TOA_times.append((temp_time.jd1, temp_time.jd2+TOF0_jd))
        # give the information to all algorithms
        true_positive = np.zeros(len(auth_algs))
        false_positive = np.zeros(len(auth_algs))
        true_negative = np.zeros(len(auth_algs))
        false_negative = np.zeros(len(auth_algs))
        correct_src = np.zeros(len(auth_algs))
        distances = np.zeros(len(auth_algs))
        for i_alg in range(len(auth_algs)):
            if len(receiver_list) < auth_algs[i_alg].minimal_receivers():
                # too few receivers
                pass
            else:
                is_auth, index_est, est_pos = auth_algs[i_alg].authenticate_satellite(receiver_list, TOA_times, TDOA_measurements, self.tle_calculator)
                if is_auth:
                    if sender_is_valid:
                        true_positive[i_alg] += 1
                    else:
                        false_positive[i_alg] += 1
                else:
                    if sender_is_valid:
                        false_negative[i_alg] += 1
                    else:
                        true_negative[i_alg] += 1
                if est_pos is not None:
                    distances[i_alg] = np.sqrt(np.sum(np.power(est_pos - sender_TEME_positions[-1], 2)))
                if index_est is None:
                    # an attack was found
                    if 'attack' in sender_sat_name:
                        correct_src[i_alg] += 1
                else:
                    est_name = self.tle_calculator.sat_names[index_est]
                    if sender_sat_name == est_name:
                        correct_src[i_alg] += 1
        return true_positive, false_negative, false_positive, true_negative, correct_src, distances




