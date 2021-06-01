from Authentication_Algorithms.Authentication_Algorithm_Interface import Authentication_Algorithm
from Authentication_Algorithms.TDOA_filtered_comparator import TDOA_filtered_comparator
from TLEcalculator import TLEcalculator
import ITRSconverter

import astropy.units as u
from astropy import coordinates as coord
import numpy as np
from scipy.spatial.transform import Rotation

# This uses the TDOA_filtered_comparator to get a first estimate. Then is validates this estimate by checking the
# area between the centroid and the estimated satellite.
class TDOA_filter_grid(Authentication_Algorithm):
    def __init__(self, valid_sat_names: [str], used_receivers: int = None, used_measurements: int = None):
        self.speed_of_light = 299792458 # meter/sec
        self.used_receivers = used_receivers
        self.valid_sat_names = valid_sat_names
        self.used_measurements = used_measurements
        self.height_threshold = 100  # km
        self.chg_threshold = 0.4  # disparity-change threshold

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



    def __make_between_grid(self, start_pos: (float, float, float), end_pos: (float, float, float), lenth: float = 1.0) -> [(float, float, float)]:
        # create a normalized mesh with x=[0,1], y=[-0.05,0.05], z=[-0.05,0.05]
        x_mesh = np.linspace(0, 1, num=20)  # the axis that will be rotated towards the end_pos
        y_mesh = np.linspace(-0.05, 0.05, num=3)
        z_mesh = np.linspace(-0.05, 0.05, num=3)
        mesh = np.array(np.meshgrid(x_mesh, y_mesh, z_mesh)).reshape(3, len(x_mesh) * len(y_mesh) * len(z_mesh)).T
        vector_st_en = end_pos - start_pos
        magnitude_st_en = np.linalg.norm(vector_st_en)
        vector_st_en = vector_st_en/magnitude_st_en
        # first make a rotation around the y axis (in the xz-plane), until the z component reaches zero
        magnitude_xz = np.linalg.norm([vector_st_en[0], vector_st_en[2]])
        rate_y = vector_st_en[0] / magnitude_xz
        angle_y = np.arccos(rate_y)
        if vector_st_en[2] < 0:
            angle_y = angle_y * -1
        # rotate the target vector this amount
        rot_y = Rotation.from_rotvec(angle_y * np.array([0, 1, 0]))
        vector_st_en_2 = rot_y.apply(vector_st_en)
        # now make a rotation around the z axis (in the xy-plane), until the rotated target vector matches x=1
        magnitude_xy_2 = np.linalg.norm([vector_st_en_2[0], vector_st_en_2[1]])
        rate_z = vector_st_en_2[0] / magnitude_xy_2
        angle_z = np.arccos(rate_z)
        if vector_st_en_2[1] > 0:
            angle_z = angle_z * -1
        # make the rotations of the mesh in the opposite direction and opposite order (first z, than y)
        rot_z = Rotation.from_rotvec(-angle_z * np.array([0, 0, 1]))
        mesh = rot_z.apply(mesh)
        rot_y = Rotation.from_rotvec(-angle_y * np.array([0, 1, 0]))
        mesh = rot_y.apply(mesh)
        # change the size of the mesh and add the start point
        mesh = mesh * (magnitude_st_en * lenth)
        mesh = mesh + start_pos
        # append the end_pos (so the algorithm has a reference 'best point' to compare against)
        mesh = np.append(mesh, [end_pos], axis=0)
        return mesh

    def minimal_receivers(self):
        if self.used_receivers is None:
            return 4
        else:
            return self.used_receivers

    def authenticate_satellite(self, receivers: [coord.ITRS], receiving_jDay: [(float, float)],
                               TDOA_measurements: [[float]], tle_calculator: TLEcalculator) -> (bool, int, (float, float, float)):
        # limit the number of used receivers
        receivers, TDOA_measurements = self.__limit_receivers(receivers, TDOA_measurements)
        receiving_jDay, TDOA_measurements = self.__limit_measurements(receiving_jDay, TDOA_measurements)
        # give the information to the TDOA_filtered_comparator
        fil_com = TDOA_filtered_comparator(self.valid_sat_names, self.used_receivers)
        est_auth, est_index, est_pos = fil_com.authenticate_satellite(receivers, receiving_jDay, TDOA_measurements, tle_calculator)
        # make this for every measurement
        rec_pos_itrs = self.__get_rec_pos(receivers)
        rec_pos_itrs = np.array(rec_pos_itrs)
        centroid_pos_itrs = np.sum(rec_pos_itrs, axis=0)/len(rec_pos_itrs)
        rec_all_pos = []
        transformed_all_mesh = []
        for i_meas in range(len(receiving_jDay)):
            # Make a grid of 200 x 10 x 10 (H x W x L). It starts at the centroid and goes to the satellite.
            temp_time = receiving_jDay[i_meas]
            err, end_pos_TEME, vel = tle_calculator.calculate_position_single(tle_calculator.satList[est_index], temp_time[0], temp_time[1])
            start_pos_TEME, vel = ITRSconverter.ITRS_2_TEME(centroid_pos_itrs, temp_time[0], temp_time[1])
            # make the grid so long, that it's end will be at 150 km height
            end_pos_itrs = ITRSconverter.TEME_2_ITRS(end_pos_TEME, [0,0,0], temp_time[0], temp_time[1])
            end_lon, end_lat, end_height = ITRSconverter.ITRS_2_LonLatHeight(end_pos_itrs)  # height will be much above 100 km
            length_factor = self.height_threshold/end_height
            # make the grid
            transformed_list = self.__make_between_grid(start_pos_TEME, end_pos_TEME, length_factor)  #shape(7220,3)
            transformed_all_mesh.append(transformed_list)
            rec_list = []
            for i_rec in range(len(rec_pos_itrs)):
                rec_pos_TEME, vel = ITRSconverter.ITRS_2_TEME(rec_pos_itrs[i_rec], temp_time[0], temp_time[1])  #shape(3)
                rec_list.append(rec_pos_TEME)  #shape(rec,3)
            rec_all_pos.append(rec_list)
        transformed_all_mesh = np.array(transformed_all_mesh)  #shape(meas, 7220, 3)
        rec_all_pos = np.array(rec_all_pos)  #shape(meas, rec, 3)
        dists = np.sqrt(np.power(transformed_all_mesh[:,:,0][:,:,None] - rec_all_pos[:,:,0][:,None,:], 2) +
                        np.power(transformed_all_mesh[:,:,1][:,:,None] - rec_all_pos[:,:,1][:,None,:], 2) +
                        np.power(transformed_all_mesh[:,:,2][:,:,None] - rec_all_pos[:,:,2][:,None,:], 2))  #shape(meas,7220,rec)
        signal_speed = self.speed_of_light / 1000  # m/s in km/s
        TOFs = dists / signal_speed  #shape(meas, 7221, rec)
        TOF_zero = TOFs[:,:,0]
        TDOAs = TOFs - TOF_zero[:,:,None]
        diffs = TDOAs - TDOA_measurements[:,None,:]  #shape(meas, 7220, rec)
        disparity = np.sqrt(1.0 / (len(receivers) * len(receiving_jDay)) * np.sum(np.sum(np.power(diffs, 2), axis=2), axis=0))  #shape(7220)
        # get most similar mesh-point
        coords = np.where(disparity == np.amin(disparity))
        index_list = list(coords[0])
        best_index = index_list[0]
        similar_pos_TEME = transformed_all_mesh[0,best_index,:]  #shape(3)
        start_time = receiving_jDay[0]
        # If the similar_position is below 100 km altitude and the signature is much better, it is an attacker.
        similar_pos_ITRS = ITRSconverter.TEME_2_ITRS(similar_pos_TEME, np.array([0,0,0]), start_time[0], start_time[1])
        sim_lon, sim_lat, sim_height = ITRSconverter.ITRS_2_LonLatHeight(similar_pos_ITRS)
        chg_ratio = 1-(disparity[best_index]/disparity[-1])
        if sim_height < self.height_threshold and chg_ratio > self.chg_threshold:
            return False, None, similar_pos_TEME
        else:
            return est_auth, est_index, est_pos





