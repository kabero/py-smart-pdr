import csv
import numpy as np

units = {"milli": 10**3, "micro": 10**6, "nano": 10**9, "pico": 10**12}


def read_data(path, unit="micro", starts_with_zero=True):
    with open(path) as f:
        reader = csv.reader(f, delimiter=" ")
        data = [list(map(float, row))
                for i, row in enumerate(reader) if i >= 1]
    data = np.array(data)
    if starts_with_zero:
        data[:, 0] -= data[0, 0]
    data[:, 0] /= units[unit]
    return data


if __name__ == "__main__":
    path = "./data/huayi_handheld1/gyro.txt"
    data = read_data(path)
    print(data.shape)
