import numpy as np
import Combine_files as cf
import Contour_plots as cp

if __name__ == '__main__':
    rec_amounts = [4,6,8,10,12,14,16,18,20]
    measurements_per_satellite = [1, 3, 5, 7, 9, 11, 15, 20, 25, 30]
    rec_dia_avg = np.array([1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3])
    algs_names = ['Comp-Grid']
    data_shape = (len(measurements_per_satellite), len(rec_amounts), len(rec_dia_avg), len(algs_names))
    entries = len(algs_names) * len(rec_dia_avg) * len(rec_amounts) * len(measurements_per_satellite)
    input_range = [1,2,3,4,5,6]
    beta_value = 1.0 / 5.0
    save_svg = False

    # make the iridium-plot
    open_path_acc = "./data/experiment1/run_iridium"  # not the whole folder-name, no end-number!
    safe_path_acc = "./data/experiment1/"
    cf.accumulate_data_sets(open_path_acc, safe_path_acc, input_range, data_shape)
    cp.make_receiver_colorplot(rec_amounts, measurements_per_satellite, data_shape, safe_path_acc, beta_value, save_svg)
    cp.make_diameter_colorplot(rec_dia_avg, measurements_per_satellite, data_shape, safe_path_acc, beta_value, save_svg)

    # make the starlink-plot
    open_path_acc = "./data/experiment1/run_starlink"  # not the whole folder-name, no end-number!
    safe_path_acc = "./data/experiment1/"
    cf.accumulate_data_sets(open_path_acc, safe_path_acc, input_range, data_shape)
    cp.make_receiver_colorplot(rec_amounts, measurements_per_satellite, data_shape, safe_path_acc, beta_value, save_svg)
    cp.make_diameter_colorplot(rec_dia_avg, measurements_per_satellite, data_shape, safe_path_acc, beta_value, save_svg)