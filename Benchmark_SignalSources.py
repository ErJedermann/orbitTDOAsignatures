from Benchmark_TLE_TDOAs import Benchmark_TLE_TDOAs
from Authentication_Algorithms.Authentication_Algorithm_Interface import Authentication_Algorithm
from TimingErrors.NormalDistributed_DynamicError import NormalDistributed_DynamicError
from Signal_Source.SignalSourceInterface import SignalSource

import numpy as np


class Benchmark_SingnalSource:
    def __init__(self, repetitions: int, total_TLE_file: str):
        self.tle_file = total_TLE_file
        self.repetitions = repetitions
        self.ben_core = Benchmark_TLE_TDOAs(self.tle_file)
        self.ssList = []
        self.ssUtilization = []
        self.eval_alg_list = []
        self.eval_alg_names = []

    def __prepare_utilization(self, special_utilization: [float]):
        special_utilization = np.array(special_utilization)
        if np.sum(special_utilization) <= 1.0:
            # the utilization is given in percent
            special_repetitions = np.array(special_utilization * self.repetitions).astype(int)
        elif np.sum(special_utilization) == self.repetitions:
            # the utilization is given in repetitions
            special_repetitions = np.array(special_utilization).astype(int)
        else:
            # the utilization is given in a relative rate
            special_repetitions = np.array(special_utilization / np.sum(special_utilization) * self.repetitions).astype(int)
        return special_repetitions

    def make_log_time_scale(self, start: float, stop: float, amount: int):
        return np.logspace(start, stop, num=amount, endpoint=True, base=10)

    def add_SignalSource(self, sSource: SignalSource, utilization: float):
        self.ssList.append(sSource)
        self.ssUtilization.append(utilization)

    def set_Algorithms(self, eval_algs: [Authentication_Algorithm], eval_algs_names: [str]):
        self.eval_alg_list = eval_algs
        self.eval_alg_names = eval_algs_names

    # This method is used to compare dynamic ranges of measurement amounts, receiver amounts and receiver diameter.
    def measurement_dynamics(self, errors_std: float, measurements_amounts: [int], measurements_time_interval: float,
                             receiver_amounts: [int], receiver_dias_min: [float], receiver_dias_max: [float],
                             save_folder_path: str):
        # create the error model
        error_model = NormalDistributed_DynamicError(0, errors_std)
        # convert the special_utilizations to special_repetitions
        special_repetitions = self.__prepare_utilization(self.ssUtilization)
        # prepare the data-collection
        true_positive_list = [[] for i in range(len(measurements_amounts))]  # true_positive_list[measurements][receivers][diameter][alg] = value
        false_positive_list = [[] for i in range(len(measurements_amounts))]
        false_negative_list = [[] for i in range(len(measurements_amounts))]
        true_negative_list = [[] for i in range(len(measurements_amounts))]
        correct_list = [[] for i in range(len(measurements_amounts))]  # correct_list[measurements][receivers][diameter][alg] = value
        dist_total_list = [[] for i in range(len(measurements_amounts))]  # dist_total_list[measurements][receivers][diameter][alg] = avg distance in km
        # go through each measurement_amount
        for i_meas in range(len(measurements_amounts)):
            #print(f"VERBOSE:Ben_sig.md.measurements: {i_meas + 1}/{len(measurements_amounts)}.")
            for i_rec in range(len(receiver_amounts)):
                #print(f"VERBOSE:Ben_sig.md.rec_am: {i_rec + 1}/{len(receiver_amounts)}.")
                temp_tp = []  # temp_true_positive_list[diameter][alg] = value
                temp_fn = []
                temp_fp = []
                temp_tn = []
                temp_co = []
                temp_dist = []
                for i_dia in range(len(receiver_dias_min)):
                    print(f"VERBOSE:Ben_sig.md.measurements - rec_am - rec_dia: {i_meas + 1}/{len(measurements_amounts)} - {i_rec + 1}/{len(receiver_amounts)} - {i_dia + 1}/{len(receiver_dias_min)}.")
                    tp_values = np.zeros(len(self.eval_alg_list))
                    fn_values = np.zeros(len(self.eval_alg_list))
                    fp_values = np.zeros(len(self.eval_alg_list))
                    tn_values = np.zeros(len(self.eval_alg_list))
                    co_values = np.zeros(len(self.eval_alg_list))
                    dist_values = np.zeros(len(self.eval_alg_list))
                    for i_spec in range(len(special_repetitions)):
                        curr_repetition = special_repetitions[i_spec]
                        curr_ss = self.ssList[i_spec]
                        curr_val = curr_ss.is_authentic()
                        for i_rep in range(curr_repetition):
                            receiver_list = self.ben_core.get_receivers(receiver_amounts[i_rec], receiver_dias_min[i_dia],
                                                                        receiver_dias_max[i_dia])
                            measurement_times = self.ben_core.get_measurement_times(measurements_amounts[i_meas],
                                                                                    measurements_time_interval)
                            sender_found, sender_index, sender_pos = curr_ss.get_positions_TEME(receiver_list, measurement_times)
                            if sender_found:
                                sender_valid = self.ben_core.verify_sender_visibility(receiver_list, measurement_times, sender_pos)
                            else:
                                sender_valid = False
                                measurement_times = self.ben_core.get_measurement_times(measurements_amounts[i_meas], measurements_time_interval)
                            while not sender_valid:
                                sender_found, sender_index, sender_pos = curr_ss.get_positions_TEME(receiver_list, measurement_times)
                                if sender_found:
                                    sender_valid = self.ben_core.verify_sender_visibility(receiver_list, measurement_times, sender_pos)
                                else:
                                    sender_valid = False
                                    measurement_times = self.ben_core.get_measurement_times( measurements_amounts[i_meas], measurements_time_interval)
                            sender_name = curr_ss.get_names()[sender_index]
                            # execute the simulation
                            tp, fn, fp, tn, corr, dist = self.ben_core.perform_mesaurement(receiver_list,
                                                                                           self.eval_alg_list, sender_pos,
                                                                                           measurement_times, curr_val,
                                                                                           sender_name,
                                                                                           error_model)
                            tp_values = tp_values + tp
                            fn_values = fn_values + fn
                            fp_values = fp_values + fp
                            tn_values = tn_values + tn
                            co_values = co_values + corr
                            dist_values = dist_values + dist
                    temp_tp.append(tp_values)
                    temp_fn.append(fn_values)
                    temp_fp.append(fp_values)
                    temp_tn.append(tn_values)
                    temp_co.append(co_values)
                    temp_dist.append(dist_values)
                true_positive_list[i_meas].append(temp_tp)
                false_positive_list[i_meas].append(temp_fp)
                false_negative_list[i_meas].append(temp_fn)
                true_negative_list[i_meas].append(temp_tn)
                correct_list[i_meas].append(temp_co)
                dist_total_list[i_meas].append(temp_dist)
        # prepare the results
        correct_list = np.array(correct_list)  # correct_list[measurements][receivers][diameter][alg] = value
        correct_list = correct_list / self.repetitions * 100  # for values in [%]
        true_positive_list = np.array(true_positive_list)
        false_positive_list = np.array(false_positive_list)
        false_negative_list = np.array(false_negative_list)
        true_negative_list = np.array(true_negative_list)
        # save results in files
        entries = len(receiver_dias_max)*len(self.eval_alg_list)*len(measurements_amounts)*len(receiver_amounts)
        true_positive_list = true_positive_list.reshape(entries)
        false_positive_list = false_positive_list.reshape(entries)
        false_negative_list = false_negative_list.reshape(entries)
        true_negative_list = true_negative_list.reshape(entries)
        correct_list = correct_list.reshape(entries)
        np.savetxt(save_folder_path+"/data_tp.csv", true_positive_list, delimiter=",")
        np.savetxt(save_folder_path+"/data_fp.csv", false_positive_list, delimiter=",")
        np.savetxt(save_folder_path+"/data_fn.csv", false_negative_list, delimiter=",")
        np.savetxt(save_folder_path+"/data_tn.csv", true_negative_list, delimiter=",")
        np.savetxt(save_folder_path+"/data_acc.csv", correct_list, delimiter=",")




