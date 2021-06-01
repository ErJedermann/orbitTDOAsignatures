


import plotly.graph_objects as go
import plotly.colors


import numpy as np

colors = plotly.colors.DEFAULT_PLOTLY_COLORS

def plot_classification_line_graph(x_elements: [float], x_label: str, alg_names: [str],figure_name: str,
                              true_positive: [float, float], false_positive: [float, float],
                              false_negative: [float, float], true_negative: [float, float], save_svg: bool = False):
    # shape(true_positive) = (alg, x_el) = fp = fn = tn
    x_elements = np.array(x_elements)

    if true_positive is not None:
        true_positive = np.array(true_positive)
    if false_positive is not None:
        false_positive = np.array(false_positive)
    if false_negative is not None:
        false_negative = np.array(false_negative)
    if true_negative is not None:
        true_negative = np.array(true_negative)
    for alg_index in range(len(alg_names)):
        figure = go.Figure()
        if true_positive is not None:
            figure.add_trace(go.Scatter(name='correctly authenticated', x=x_elements, y=true_positive[alg_index,:], mode='lines+markers'))
        if false_negative is not None:
            figure.add_trace(go.Scatter(name='falsely rejected', x=x_elements, y=false_negative[alg_index,:], mode='lines+markers'))
        if false_positive is not None:
            figure.add_trace(go.Scatter(name='falsely authenticated', x=x_elements, y=false_positive[alg_index,:], mode='lines+markers'))
        if true_negative is not None:
            figure.add_trace(go.Scatter(name='correctly rejected', x=x_elements, y=true_negative[alg_index,:], mode='lines+markers'))
        figure.update_layout(showlegend=False,
                             xaxis_title=x_label,
                             #title=f'{figure_name} ({alg_names[alg_index]})',
                             yaxis_title='rate [%]',
                             width=500,
                             height=400,
                             font=dict(size=20),
                             margin=dict(l=10, r=10, b=10, t=10, pad=4)
                             )
        figure.show()
        if save_svg:
            plotly.offline.plot(figure, image_filename=f'{figure_name}({alg_names[alg_index]})', image='svg')
    return 0

def make_classification_plot(rec_dia_avg: [float], algs_names: [str], open_folder: str, save_svg: bool = False):
    rec_dia_avg = np.array(rec_dia_avg) / 1000  # m in km (for x-axis)
    my_shape = (len(algs_names), len(rec_dia_avg))
    data_fn = np.loadtxt(open_folder + 'data_fn.csv', delimiter=",")
    data_fp = np.loadtxt(open_folder + 'data_fp.csv', delimiter=",")
    data_tp = np.loadtxt(open_folder + 'data_tp.csv', delimiter=",")
    data_tn = np.loadtxt(open_folder + 'data_tn.csv', delimiter=",")
    data_fn = data_fn.reshape(my_shape)
    data_fp = data_fp.reshape(my_shape)
    data_tp = data_tp.reshape(my_shape)
    data_tn = data_tn.reshape(my_shape)
    # false rejection rate
    frr = data_fn / (data_tn + data_fn) * 100
    # false authentication rate
    far = data_fp / (data_tp + data_fp) * 100
    x_axis = rec_dia_avg
    x_label = 'distribution diameter [km]'
    name = f'frr_and_far'
    plot_classification_line_graph(x_axis, x_label, algs_names, name, frr, far, None, None, save_svg)


