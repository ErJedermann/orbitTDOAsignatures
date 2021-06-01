import numpy as np
import Combine_files as cf
import Classification_plots as cp

if __name__ == '__main__':
    rec_amounts = [6]
    measurements_per_satellite = [5]
    rec_dia_avg = np.array([1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3])
    algs_names = ['Comp-Grid']
    data_shape = (len(measurements_per_satellite), len(rec_amounts), len(rec_dia_avg), len(algs_names))
    entries = len(algs_names) * len(rec_dia_avg) * len(rec_amounts) * len(measurements_per_satellite)
    input_range = [0,1]
    save_svg = False

    # make the iridium-plot
    open_iridium = "./data/experiment2/iridium/"
    cp.make_classification_plot(rec_dia_avg, algs_names, open_iridium, save_svg)

    # make the starlink-plot
    open_starlink = "./data/experiment2/starlink/"
    cp.make_classification_plot(rec_dia_avg, algs_names, open_starlink, save_svg)