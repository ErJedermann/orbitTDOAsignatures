from Authentication_Algorithms.Authentication_Algorithm_Interface import Authentication_Algorithm
from TLEcalculator import TLEcalculator
import ITRSconverter

import astropy.units as u
from astropy import coordinates as coord
import numpy as np


# This implementation of the TDOA_filtered_comparator computes the TDOA fingerprint of all measurements, for all
# satellites at once.
class TDOA_filtered_comparator(Authentication_Algorithm):
    def __init__(self, valid_sat_names: [str], used_receivers: int = None, used_measurements: int = None):
        self.speed_of_light = 299792458 # meter/sec
        self.used_receivers = used_receivers
        self.valid_sat_names = valid_sat_names
        self.used_measurements = used_measurements

    def __filter_satellite_invisibility(self, receivers: np.ndarray, satellites: np.ndarray) -> [bool]:
        # receivers.shape = (x,3) receivers[rec_x, pos_{xyz}]
        # satellites.shape = (y,3) satellites[sat_y, pos_{xyz}]
        # returns true: satellite is invisible (should be punished), false: satellite is visible (not punished)
        # idea: if the scalar-product of rec * (sat-rec) is < zero, the satellite is behind the horizon (not visible)
        # get the vector (sat-rec)
        vec_x = receivers[:,0][:,None] - satellites[:,0][None,:]  # shape (x,y)
        vec_y = receivers[:,1][:,None] - satellites[:,1][None,:]  # shape (x,y)
        vec_z = receivers[:,2][:,None] - satellites[:,2][None,:]  # shape (x,y)
        # do the dot-product (scalar-product) of dist * receivers
        scalar_pod = vec_x * receivers[:,0][:,None] + vec_y * receivers[:,1][:,None] + vec_z * receivers[:,2][:,None]
        scalar_bool = np.where(scalar_pod < 0, False, True)  # shape (x,y)
        sat_bool = np.bitwise_and.reduce(scalar_bool, 0)  # shape (y)
        return sat_bool

    def __get_rec_pos(self, receivers:[coord.ITRS]) -> [(float, float, float)]:
        pos_list = []
        for i in range(len(receivers)):
            if receivers[i].cartesian.xyz.unit is u.m:
                pos_list.append(receivers[i].cartesian.xyz.value / 1000)  # [m] in [km]
            else:
                pos_list.append(receivers[i].cartesian.xyz.value)
        return pos_list

    def __limit_receivers(self, receivers: [coord.ITRS], TDOA_measurements: [[float]]):
        if self.used_receivers is not None:
            if len(receivers) > self.used_receivers:
                TDOA_measurements = np.array(TDOA_measurements)  # shape (meas, rec)
                temp_rec_itrs = []
                temp_TDOFs = []
                for i in np.linspace(0, len(receivers)-1, self.used_receivers).astype(int):
                    temp_rec_itrs.append(receivers[i])
                    temp_TDOFs.append(TDOA_measurements[:,i])
                receivers = temp_rec_itrs
                TDOA_measurements = temp_TDOFs
                TDOA_measurements = np.array(TDOA_measurements).T  # shape (rec, meas) -> (meas, rec)
        receivers = np.array(receivers)  # shape (rec)
        TDOA_measurements = np.array(TDOA_measurements)  # shape (meas, rec)
        return receivers, TDOA_measurements

    def __limit_measurements(self, receiving_jDay: [(float, float)], TDOA_measurements: [[float]]):
        if self.used_measurements is not None:
            if len(receiving_jDay) > self.used_measurements:
                TDOA_measurements = np.array(TDOA_measurements)  # shape (meas, rec)
                temp_rec_jDay = []
                temp_TDOAs = [[] for i in range(self.used_measurements)]
                counter = 0
                for i in np.linspace(0, len(receiving_jDay)-1, self.used_measurements).astype(int):
                    temp_rec_jDay.append(receiving_jDay[i])
                    temp_TDOAs[counter] = TDOA_measurements[i,:]
                    counter += 1
                receiving_jDay = temp_rec_jDay
                TDOA_measurements = temp_TDOAs  # shape (meas, rec)
                TDOA_measurements = np.array(TDOA_measurements)
        return receiving_jDay, TDOA_measurements

    def minimal_receivers(self):
        if self.used_receivers is None:
            return 4
        else:
            return self.used_receivers

    def authenticate_satellite(self, receivers: [coord.ITRS], receiving_jDay: [(float, float)],
                               TDOA_measurements: [[float]], tle_calculator: TLEcalculator) -> (bool, int, (float, float, float)):
        # limit the number of receivers and measurements
        receivers, TDOA_measurements = self.__limit_receivers(receivers, TDOA_measurements)
        receiving_jDay, TDOA_measurements = self.__limit_measurements(receiving_jDay, TDOA_measurements)
        # calculate the receiver TEME-positions at all measurement-times
        rec_pos_itrs = self.__get_rec_pos(receivers)
        rec_pos_itrs = np.array(rec_pos_itrs).T
        rec_pos_list_TEME = []
        for i in range(len(receiving_jDay)):
            temp_time = receiving_jDay[i]
            pos_teme, vel_teme = ITRSconverter.ITRS_2_TEME(rec_pos_itrs, temp_time[0], temp_time[1])
            rec_pos_list_TEME.append(pos_teme.T)
        rec_pos_list_TEME = np.array(rec_pos_list_TEME)  # shape (meas, rec, 3)
        # get the TEME-sat-positions at the desired times,  sat_pos_list_TEME.shape(sat, meas, 3)
        receiving_jDay = np.array(receiving_jDay)
        sat_pos_list_TEME, sat_vel = tle_calculator.calculate_multiple_positions_all(receiving_jDay[:, 0], receiving_jDay[:, 1])  # pos in TEME [km]
        # calculate the TOF for each sat to each rec at each measurement
        rec_x = rec_pos_list_TEME[:, :, 0][None, :, :]  # shape: (?, meas, rec)
        rec_y = rec_pos_list_TEME[:, :, 1][None, :, :]
        rec_z = rec_pos_list_TEME[:, :, 2][None, :, :]
        sat_x = sat_pos_list_TEME[:, :, 0][:, :, None]  # shape: (sat, meas, ?)
        sat_y = sat_pos_list_TEME[:, :, 1][:, :, None]
        sat_z = sat_pos_list_TEME[:, :, 2][:, :, None]
        dists = np.sqrt(np.power(rec_x - sat_x, 2) + np.power(rec_y - sat_y, 2) + np.power(rec_z - sat_z, 2))  # shape: (sat, meas, rec)
        signal_speed = self.speed_of_light / 1000  # from [m/s] in [km/s]
        TOFs = dists / signal_speed  # in [s]
        TDOFs = TOFs - TOFs[:, :, 0][:, :, None]  # shape: (sat, meas, rec)
        # calculate the disparity (the "root mean squared error")
        diffs = TDOFs - TDOA_measurements[None, :, :]  # shape: (sat, meas, rec)
        disparity = np.sqrt( 1.0 / (len(receivers) * len(receiving_jDay)) * np.sum(np.sum(np.power(diffs, 2), axis=2), axis=1))  # shape:(sat)
        # calculate the visibilities of the satellites at the begin and at the end
        sat_invis_start = self.__filter_satellite_invisibility(rec_pos_list_TEME[0,:,:], sat_pos_list_TEME[:,0,:])
        sat_invis_end = self.__filter_satellite_invisibility(rec_pos_list_TEME[-1,:,:], sat_pos_list_TEME[:,-1,:])
        sat_invis_total = np.bitwise_or(sat_invis_start, sat_invis_end)
        dis_max = np.max(disparity)
        disparity = disparity + (sat_invis_total * dis_max)
        # at the end get the smallest with the smallest over all disparity
        coords = np.where(disparity == np.amin(disparity))
        index_list = list(coords[0])
        best_index = index_list[0]
        # check if the sat at best_index has a valid name
        est_name = tle_calculator.sat_names[best_index]
        if any(est_name in s for s in self.valid_sat_names):
            return True, best_index, sat_pos_list_TEME[best_index,0,:]
        else:
            return False, best_index, sat_pos_list_TEME[best_index,0,:]



