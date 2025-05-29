import os
import numpy as np

def read_np4x4_from_file(filename):
  floats = []
  with open(filename, 'r') as file:
    for line in file:
        values = line.split()
        floats.extend(values)
        if len(floats) >= 16:
            break
  return ' '.join(floats[:16])

def string_to_np4x4(string_of_floats):
    if string_of_floats is None:
        return None
    if os.path.isfile(string_of_floats):
        float_strings = read_np4x4_from_file(string_of_floats).split()
    else:
        float_strings = string_of_floats.split()
    float_values = [float(x) for x in float_strings]
    if len(float_values) != 16:
        raise ValueError("Input string must contain exactly 16 float values")
    np_4x4 = np.array(float_values).reshape(4, 4)
    return np_4x4

def np4x4_to_string(array):
    # convert to string (the dumb way)
    string = ''
    for r in xrange(4):
        for c in xrange(4):
            string += "%0.9f " % array[r][c]
    return string
