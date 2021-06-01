import plotly.graph_objects as go
import plotly.colors
import numpy as np

colors = plotly.colors.DEFAULT_PLOTLY_COLORS

def plot_2D_colorimage(x_numbers: [float], x_label: str, y_numbers: [], y_label: str, data: [[float]], min_value: float, max_value: float, figure_name: str, do_print: bool):
    # data.shape: (x,y)
    # make the plot
    data = np.array(data)
    length_x, length_y = np.shape(data)
    fig = go.Figure(go.Contour(x=x_numbers, y=y_numbers, z=data,colorscale='Rainbow',
                               contours=dict(start=min_value,end=max_value,
                                             #size=0.01,
                                             ),
                               ))
    # Layout
    fig.update_layout(
        #title=figure_name,
        width=600,
        height=400,
        font=dict(size=20),
        margin=dict(l=10, r=10, b=10, t=10, pad=4)
    )
    fig.update_layout(xaxis_title=x_label)
    fig.update_yaxes(title_text=y_label)
    fig.show()
    if do_print:
        plotly.offline.plot(fig, image_filename=f'{figure_name}', image='svg')

def make_receiver_colorplot(rec_amounts: [int], meas_per_sat: [int], data_shape, open_path: str, beta_value: float, save_svg: bool):
    data_fn = np.loadtxt(open_path + 'data_accumulated_fn.csv', delimiter=",")
    data_fp = np.loadtxt(open_path + 'data_accumulated_fp.csv', delimiter=",")
    data_tp = np.loadtxt(open_path + 'data_accumulated_tp.csv', delimiter=",")
    data_tn = np.loadtxt(open_path + 'data_accumulated_tn.csv', delimiter=",")
    data_fn = data_fn.reshape(data_shape)
    data_fp = data_fp.reshape(data_shape)
    data_tp = data_tp.reshape(data_shape)
    data_tn = data_tn.reshape(data_shape)
    data_sum = ((1+beta_value**2)*data_tp) / ((1+beta_value**2)*data_tp + beta_value**2 * data_fn + data_fp)
    x_axis = rec_amounts
    x_label = 'number of receivers'
    y_axis = meas_per_sat
    y_label = 'number of messages'
    #min_value = np.min(data_sum)  # original
    #max_value = np.max(data_sum)  # original
    min_value = 0.771156462585034  # to have the same colorscale on all plots
    max_value = 0.9998713826366559  # to have the same colorscale on all plots
    data = data_sum[:, :, 0, 0]
    print(f"{np.shape(data)}")
    name = f'rec_response_variable'
    plot_2D_colorimage(x_axis, x_label, y_axis, y_label, data, min_value, max_value, name, save_svg)

def make_diameter_colorplot(rec_dia: [int], meas_per_sat: [int], data_shape, open_path: str, beta_value: float, save_svg: bool):
    rec_dia = np.array(rec_dia) / 1000  # m in km
    data_fn = np.loadtxt(open_path + 'data_accumulated_fn.csv', delimiter=",")
    data_fp = np.loadtxt(open_path + 'data_accumulated_fp.csv', delimiter=",")
    data_tp = np.loadtxt(open_path + 'data_accumulated_tp.csv', delimiter=",")
    data_tn = np.loadtxt(open_path + 'data_accumulated_tn.csv', delimiter=",")
    data_fn = data_fn.reshape(data_shape)
    data_fp = data_fp.reshape(data_shape)
    data_tp = data_tp.reshape(data_shape)
    data_tn = data_tn.reshape(data_shape)
    data_sum = ((1+beta_value**2)*data_tp) / ((1+beta_value**2)*data_tp + beta_value**2 * data_fn + data_fp)
    x_axis = rec_dia
    x_label = 'distribution diameter [km]'
    y_axis = meas_per_sat
    y_label = 'number of messages'
    #min_value = np.min(data_sum)  # original
    #max_value = np.max(data_sum)  # original
    min_value = 0.771156462585034  # to have the same colorscale on all plots
    max_value = 0.9998713826366559  # to have the same colorscale on all plots
    data = data_sum[:, 0, :, 0]
    print(f"{np.shape(data)}")
    name = f'dia_response_variable'
    plot_2D_colorimage(x_axis, x_label, y_axis, y_label, data, min_value, max_value, name, save_svg)
