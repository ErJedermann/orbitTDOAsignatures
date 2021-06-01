
import numpy as np

def load_array(shape: [int], file:str):
    input_data = np.loadtxt(file, delimiter=",")
    input_data = input_data.reshape(shape)
    return input_data

# sum up the values of multiple arrays into one array
def accumulate_arrays(file_path1: str, numbers: [int], file_path2:str, shape: [int]):
    data = np.zeros(shape)
    for i in range(len(numbers)):
        temp_path = file_path1+str(numbers[i])+file_path2
        temp_data = load_array(shape, temp_path)
        data = data + temp_data
    return data

# concatenate multiple arrays of same size behind each other
def concatenate_arrays(file_path1: str, numbers: [int], file_path2:str, shape: [int]):
    temp_path = file_path1 + str(numbers[0]) + file_path2
    data = load_array(shape, temp_path)
    for i in range(1, len(numbers)):
        temp_path = file_path1+str(numbers[i])+file_path2
        temp_data = load_array(shape, temp_path)
        data = np.concatenate((data, temp_data), axis=0)
    return data

def accumulate_data_sets(open_path: str, safe_path: str, input_range: [int], shape: [int]):
    line_shape = np.prod(shape)
    acc_fp = accumulate_arrays(open_path, input_range, '/data_fp.csv', shape)
    acc_fp = acc_fp.reshape(line_shape)
    np.savetxt(str(safe_path + "data_accumulated_fp.csv"), acc_fp, delimiter=",")
    acc_fn = accumulate_arrays(open_path, input_range, '/data_fn.csv', shape)
    acc_fn = acc_fn.reshape(line_shape)
    np.savetxt(str(safe_path + "data_accumulated_fn.csv"), acc_fn, delimiter=",")
    acc_tn = accumulate_arrays(open_path, input_range, '/data_tn.csv', shape)
    acc_tn = acc_tn.reshape(line_shape)
    np.savetxt(str(safe_path + "data_accumulated_tn.csv"), acc_tn, delimiter=",")
    acc_tp = accumulate_arrays(open_path, input_range, '/data_tp.csv', shape)
    acc_tp = acc_tp.reshape(line_shape)
    np.savetxt(str(safe_path + "data_accumulated_tp.csv"), acc_tp, delimiter=",")


def concatenate_data_sets(open_path: str, safe_path: str, input_range: [int], shape: [int]):
    line_shape = np.prod(shape) * len(input_range)
    con_tp = concatenate_arrays(open_path, input_range, '/data_tp.csv', shape)
    con_tp = con_tp.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_tp.csv", con_tp, delimiter=",")
    con_tn = concatenate_arrays(open_path, input_range, '/data_tn.csv', shape)
    con_tn = con_tn.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_tn.csv", con_tn, delimiter=",")
    con_fp = concatenate_arrays(open_path, input_range, '/data_fp.csv', shape)
    con_fp = con_fp.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_fp.csv", con_fp, delimiter=",")
    con_fn = concatenate_arrays(open_path, input_range, '/data_fn.csv', shape)
    con_fn = con_fn.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_fn.csv", con_fn, delimiter=",")


def concatenate_both_data_sets(open_path1: str, open_path2: str, safe_path: str, input_range: [int], shape: [int]):
    line_shape = np.prod(shape) * len(input_range) * 2
    con_tp1 = concatenate_arrays(open_path1, input_range, '/data_tp.csv', shape)
    con_tp2 = concatenate_arrays(open_path2, input_range, '/data_tp.csv', shape)
    con_tp = np.concatenate((con_tp1, con_tp2), axis=0)
    con_tp = con_tp.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_both_tp.csv", con_tp, delimiter=",")
    con_tn1 = concatenate_arrays(open_path1, input_range, '/data_tn.csv', shape)
    con_tn2 = concatenate_arrays(open_path2, input_range, '/data_tn.csv', shape)
    con_tn = np.concatenate((con_tn1, con_tn2), axis=0)
    con_tn = con_tn.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_both_tn.csv", con_tn, delimiter=",")
    con_fp1 = concatenate_arrays(open_path1, input_range, '/data_fp.csv', shape)
    con_fp2 = concatenate_arrays(open_path2, input_range, '/data_fp.csv', shape)
    con_fp = np.concatenate((con_fp1, con_fp2), axis=0)
    con_fp = con_fp.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_both_fp.csv", con_fp, delimiter=",")
    con_fn1 = concatenate_arrays(open_path1, input_range, '/data_fn.csv', shape)
    con_fn2 = concatenate_arrays(open_path2, input_range, '/data_fn.csv', shape)
    con_fn = np.concatenate((con_fn1, con_fn2), axis=0)
    con_fn = con_fn.reshape(line_shape)
    np.savetxt(safe_path + "data_concatenated_both_fn.csv", con_fn, delimiter=",")
