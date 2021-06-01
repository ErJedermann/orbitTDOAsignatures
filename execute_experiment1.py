from TLEcalculator import TLEcalculator
from Benchmark_SignalSources import Benchmark_SingnalSource
from Signal_Source.Satellite import SatelliteSource
from Signal_Source.Drone import DroneSource
from Authentication_Algorithms.TDOA_filter_grid import TDOA_filter_grid
import time
import numpy as np

if __name__ == '__main__':
    # Make a run of experiment 1
    all_tle_file = "./data/TLEs/active_2021_01_16.tle"
    repetitions = 900
    ben_sig = Benchmark_SingnalSource(repetitions, all_tle_file)
    ir_sats = True
    if ir_sats:
        folder_path = "./data/experiment1/run_iridium1/"
        iridium_tle_file = './data/TLEs/Iridium.tle'
        iridium_source = SatelliteSource(True, iridium_tle_file)
        iridium_source.set_deviation(30)
        ben_sig.add_SignalSource(iridium_source, 300)
        leo_no_ir_tle_file = "./data/TLEs/LEO_no_iridium.tle"
        leo_no_ir_source = SatelliteSource(False, leo_no_ir_tle_file)
        leo_no_ir_source.set_deviation(30)
        ben_sig.add_SignalSource(leo_no_ir_source, 300)
        valid_sats_tle = TLEcalculator(iridium_tle_file)
    else:
        folder_path = "./data/experiment1/run_starlink1/"
        starlink_tle_file = './data/TLEs/Starlink.tle'
        starlink_source = SatelliteSource(True, starlink_tle_file)
        starlink_source.set_deviation(30)
        ben_sig.add_SignalSource(starlink_source, 300)
        leo_no_sl_tle_file = "./data/TLEs/LEO_no_starlink.tle"
        leo_no_sl_source = SatelliteSource(False, leo_no_sl_tle_file)
        leo_no_sl_source.set_deviation(30)
        ben_sig.add_SignalSource(leo_no_sl_source, 300)
        valid_sats_tle = TLEcalculator(starlink_tle_file)
    drone_source = DroneSource(False, 30, 18000, 20)
    ben_sig.add_SignalSource(drone_source, 300)
    comp_grid = TDOA_filter_grid(valid_sats_tle.sat_names)
    algs = [comp_grid]
    algs_names = ['Comp-grid'] #'Comp-Fil'
    ben_sig.set_Algorithms(algs, algs_names)
    rec_amounts = [4,6,8,10,12,14,16,18,20]  # for exp1
    measurements_per_satellite = [1, 3, 5, 7, 9, 11, 15, 20, 25, 30]  # for exp1
    measurements_time_interval = 0.1
    rec_dia_avg = np.array([1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3])  # used fo both experiments
    rec_dia_min = rec_dia_avg * 0.9
    rec_dia_max = rec_dia_avg * 1.1
    errors_std = 100e-9  # 100ns error
    t_start = time.time()
    ben_sig.measurement_dynamics(errors_std, measurements_per_satellite, measurements_time_interval, rec_amounts, rec_dia_min, rec_dia_max, folder_path)
    t_end = time.time()
    simulations = repetitions * len(rec_amounts) * len(measurements_per_satellite) * len(rec_dia_avg)
    print(f"FINISHED: {simulations} simulations in {t_end-t_start} seconds.")