__author__ = 'german'

import math

def mean(array):
    return sum(map(float, array)) / len(array)

def var(array):
    sample_mean_squared = mean(array) ** 2
    return sum(list([float(x) ** 2 - sample_mean_squared for x in array])) / (len(array) - 1)

def ssd(array):
    return math.sqrt(var(array))

def median(array):
    sorted_list = sorted(array)
    list_len = len(array)
    index = (list_len - 1) // 2
    if (list_len % 2):
        return sorted_list[index]
    else:
        return (sorted_list[index] + sorted_list[index + 1]) / 2.0


def medianW(input_array):
    array = []
    for i in xrange(len(input_array)):
        for j in xrange(0, i + 1):
            array.append((1.0 * (input_array[i] + input_array[j])) / 2)

    sorted_list = sorted(array)
    list_len = len(array)
    index = (list_len - 1) // 2
    if (list_len % 2):
        return sorted_list[index]
    else:
        return (sorted_list[index] + sorted_list[index + 1]) / 2.0

def himed(array):
    sorted_list = sorted(array)
    list_len = len(array)
    index = int(math.ceil(list_len / 2) + 1)
    return sorted_list[index - 1]

def lomed(array):
    sorted_list = sorted(array)
    list_len = len(array)
    index = int(math.ceil((list_len + 1) / 2))
    return sorted_list[index - 1]

def sn_estimator(array):
    result_array = []
    multiplier = 1.0
    for elem in array:
        tmp_list = [abs(elem - x) for x in array]
        result_array.append(himed(tmp_list))
    if len(array) < 10:
        multipliers = [0.743, 1.851, 0.954, 1.351, 0.993, 1.198, 1.005, 1.131]
        multiplier = multipliers[len(array) - 1]
    elif (len(array) % 2 == 1):
        multiplier = 1.0 * len(array) / (1.0 * len(array) - 0.9)
    return 1.1926 * multiplier * lomed(result_array)


