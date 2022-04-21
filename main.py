# a = [ax, ay, az]: 3軸加速度(LCS)
# b = [bx, by, bz]: 3軸地磁気(LCS)
# a_gcs = [ax_gcs, ay_gcs, az_gcs]: 3軸加速度(GCS)

from tkinter import W
import numpy as np
import matplotlib.pyplot as plt
from collections import deque
import pandas as pd

from utils import read_data
import show_graph
from params import *
from detect_step import extract_t_peak


def calc_roll(ax, az):
    return np.arctan(-ax / az)


def calc_pitch(ax, ay, az):
    return np.arctan(ay / np.sqrt(ax ** 2 + az ** 2))


def calc_yaw(bx, by, bz, roll, pitch):
    molecule = -bx * np.sin(pitch) * np.sin(roll) - by * \
        np.cos(roll) + bz * np.sin(roll) * np.cos(pitch)
    denominator = bx * np.cos(roll) + bz * np.sin(pitch)
    return np.arctan(molecule / denominator)


def calc_euler(a, b):
    roll = calc_roll(a[0], a[2])
    pitch = calc_pitch(*a)
    yaw = calc_yaw(*b, roll, pitch)
    return roll, pitch, yaw


def calc_rotations_matrix(roll, pitch, yaw, is_devided=False):
    R_roll = np.array([[np.cos(roll), 0, np.sin(roll)],
                       [0, 1, 0],
                       [-np.sin(roll), 0, np.cos(roll)]])
    R_pitch = np.array([[1, 0, 0],
                        [0, -np.cos(pitch), np.sin(pitch)],
                        [0, np.sin(pitch), np.cos(pitch)]])
    R_yaw = np.array([[np.cos(yaw), np.sin(yaw), 0],
                      [-np.sin(yaw), np.cos(yaw), 0],
                      [0, 0, 1]])
    if is_devided:
        return R_roll, R_pitch, R_yaw
    else:
        return R_yaw @ R_pitch @ R_roll


def calc_a_gcs(a_lcs, mag_lcs):
    eulers = calc_euler(a_lcs, mag_lcs)
    R = calc_rotations_matrix(*eulers)
    a_gcs = (R @ np.array(a_lcs).T)
    return a_gcs


def calc_g(lastg, a_gcs):
    return ALPHA * lastg + (1 - ALPHA) * a_gcs[2]


def calc_a_hpf(a_lcs, mag_lcs, lastg):
    a_gcs = calc_a_gcs(a_lcs, mag_lcs)
    g = calc_g(lastg, a_gcs)
    return a_gcs[2] - g, lastg


def calc_m_gcs(a_lcs, mag_lcs):
    eulers = calc_euler(a_lcs, mag_lcs)
    R_roll, R_pitch, _ = calc_rotations_matrix(*eulers, is_devided=True)
    return R_roll @ R_pitch @ np.array(mag_lcs).T


def calc_h_mag(a_lcs, mag_lcs):
    # print(a_lcs)
    # print(mag_lcs)
    mag_gcs = calc_m_gcs(a_lcs, mag_lcs)
    molecule = -mag_gcs[1]
    denominator = np.sqrt(mag_gcs[0]**2 + mag_gcs[1]**2) + mag_gcs[0]
    return 2 * np.arctan(molecule / denominator) - H_DECLINE


def calc_w_hat(w_lcs, bias_gyro=[]):
    if len(bias_gyro) > 0:
        return w_lcs - bias_gyro
    else:
        return w_lcs


def calc_w_gcs(w_hat, gt):
    # print(gt, w_hat)
    # print(np.array(np.dot(gt.T, w_hat)))
    # exit()
    return np.dot(w_hat, gt) / np.linalg.norm(gt)


def calc_h_gyro(w_gcs_history):
    w_gcs_history = np.array(w_gcs_history)
    return np.sum(-w_gcs_history)


def calc_h_t(h_gyro, h_mag, h_mag_prev, h_prev):
    h_delta_cor = np.abs(h_mag - h_gyro)
    h_delta_mag = np.abs(h_mag - h_mag_prev)

    if h_delta_cor <= H_COR/180*np.pi and h_delta_mag <= H_MAG/180*np.pi:
        w_pmg = 1 / (W_PREV + W_MAG + W_GYRO)
        res = w_pmg * (W_PREV * h_prev + W_MAG * h_mag + W_GYRO * h_gyro)
    elif h_delta_cor <= H_COR/180*np.pi and h_delta_mag > H_MAG/180*np.pi:
        w_mg = 1 / (W_MAG + W_GYRO)
        res = w_mg * (W_MAG * h_mag + W_GYRO * h_gyro)
    elif h_delta_cor > H_COR/180*np.pi and h_delta_mag <= H_MAG/180*np.pi:
        res = h_prev
    else:
        w_pg = 1 / (W_PREV + W_GYRO)
        res = w_pg * (W_PREV * h_prev + W_GYRO * h_gyro)

    return res


def print_start_msg(sampling_rate):
    print("#########################")
    print("sampling rate[s]:", sampling_rate)
    print("sampling rate[Hz]:", 1 / sampling_rate)
    print("#########################")
    print()


def main():
    # データの読み込み
    accpath = "data/huayi_handheld1/acce.txt"
    magpath = "data/huayi_handheld1/magnet.txt"
    gyropath = "data/huayi_handheld1/gyro.txt"
    # accpath = "data/dan_handheld1/acce.txt"
    # magpath = "data/dan_handheld1/magnet.txt"
    # gyropath = "data/dan_handheld1/gyro.txt"
    # accpath = "data/tang_handheld1/acce.txt"
    # magpath = "data/tang_handheld1/magnet.txt"
    # gyropath = "data/tang_handheld1/gyro.txt"

    # めっちゃいい
    # accpath = "data/hao_handheld1/acce.txt"
    # magpath = "data/hao_handheld1/magnet.txt"
    # gyropath = "data/hao_handheld1/gyro.txt"

    # まったくよくない
    # accpath = "data/hao_handheld2/acce.txt"
    # magpath = "data/hao_handheld2/magnet.txt"
    # gyropath = "data/hao_handheld2/gyro.txt"

    accdata = read_data(accpath, "nano")
    magdata = read_data(magpath, "nano")
    gyrodata = read_data(gyropath, "nano")

    #　ダウンサンプリングもどき
    magdata = magdata[::SAMPLING_STEP]

    # 各種情報の表示
    sampling_rate = np.mean(np.diff(magdata[:, 0]))
    print_start_msg(sampling_rate)

    # データの整形(補間)
    t_mag = magdata[:, 0]
    downsampled_accdata = []
    downsampled_gyrodata = []
    for tm in t_mag:
        idx = np.abs(accdata[:, 0] - tm).argmin()
        downsampled_accdata.append(accdata[idx, 1:])
        downsampled_gyrodata.append(gyrodata[idx, 1:])
    downsampled_accdata = np.array(downsampled_accdata)
    downsampled_gyrodata = np.array(downsampled_gyrodata)
    magdata = magdata[:, 1:]

    # ローパス & ハイパスフィルター
    a_step = []
    t_step = []
    history = deque([])
    g = GRAVITY_ACC
    for i, (acc, mag) in enumerate(zip(downsampled_accdata, magdata)):
        a_hpf, g = calc_a_hpf(acc, mag, g)
        history.append(a_hpf)
        if len(history) <= WINDOW_SIZE:
            continue
        else:
            history.popleft()
            moving_average = np.sum(history) / WINDOW_SIZE
            t_step.append(t_mag[i - WINDOW_SIZE // 2])
            a_step.append(moving_average)
    a_step = np.array(a_step)
    t_step = np.array(t_step)

    # ピーク検出
    t_peak = extract_t_peak(t_step, a_step)
    show_graph.show_t_peaks(t_step, a_step)

    # ジャイロセンサのドリフトオフセットを計算
    # 静止状態の平均(最初と最後)
    sum_gyro = np.cumsum(downsampled_gyrodata, axis=0)
    sum_gyro_first = sum_gyro[0:STABLE_LEN, :]
    sum_gyro_end = sum_gyro[-STABLE_LEN:, :]
    delta_t = t_mag[-STABLE_LEN] - t_mag[STABLE_LEN]
    b_gyro = (sum_gyro_first[-1] - sum_gyro_end[0]) / delta_t

    # 方向推定
    tmp = (downsampled_accdata, magdata, downsampled_gyrodata)
    g = GRAVITY_ACC
    h_mag_res = []
    w_gcs_history = []
    for i, (acc, mag, gyro) in enumerate(zip(*tmp)):
        a_gcs = calc_a_gcs(acc, mag)
        g = calc_g(g, a_gcs)
        eulers = calc_euler(acc, mag)
        R = calc_rotations_matrix(*eulers, is_devided=True)
        gt = (R[0] @ R[1]).T @ np.array([[0, 0, g]]).T
        w_hat = calc_w_hat(gyro, b_gyro)  # 補正後の角速度(ドリフトバイアス補正有り)
        # w_hat = calc_w_hat(gyro)  # 補正後の角速度
        w_gcs = calc_w_gcs(w_hat, gt)  # 補正後の角速度(GCS)
        w_gcs_history.append(w_gcs)
        h_mag = calc_h_mag(acc, mag)
        h_mag_res.append(h_mag)
    w_gcs_history = np.array(w_gcs_history)

    # TODO: これは超恣意的に割っているので必ず修正！
    h_gyro_res = np.cumsum(-w_gcs_history) * 0.06
    # h_gyro_res = np.array(h_gyro_res)
    h_mag_res = np.array(h_mag_res)

    # ジャイロセンサに基づく方位推定結果と磁力センサに基づく方位推定結果を
    # 統合して導き出した方向推定結果
    h_t_res = []
    h_t = 0
    for i in range(len(h_gyro_res)):
        if i == 0:
            h_t = calc_h_t(h_gyro_res[i], h_mag_res[i], 0, h_t)
        else:
            h_t = calc_h_t(h_gyro_res[i], h_mag_res[i], h_mag_res[i-1], h_t)
        h_t_res.append(h_t)
    h_t_res = np.array(h_t_res)

    # plt.plot(h_gyro_res % (np.pi * 2) - np.pi)
    def show_heading_direction():
        plt.figure(figsize=(13, 6))
        plt.plot(h_gyro_res, label="gyro", marker="o", markersize=3)
        plt.plot(h_t_res, label="gyro + mag",
                 marker="o", markersize=3, color="g")
        plt.plot(np.unwrap(h_mag_res, discont=np.pi * 1.0),
                 color="orange", label="mag", marker="o", markersize=3)
        plt.legend()
        plt.show()

    show_heading_direction()

    # ステップ長推定
    sorted_t_peak = sorted(list(t_peak))
    a_peak = a_step[sorted_t_peak]
    t_head_list = []
    for i in range(1, len(sorted_t_peak)):
        start = sorted_t_peak[i-1] + 1
        end = sorted_t_peak[i]
        # print(start, end)
        # print(a_peak[start:end])
        # print(i)
        t_head = np.argmin(a_step[start:end]) + start
        t_head_list.append(t_head)

    # show_graph.show_t_peaks_and_valleys(t_step, a_step, t_head_list)
    a_valley = np.hstack([[0], np.array(a_step[t_head_list])])
    a_pp_step = a_peak - a_valley  # ステップイベント間隔

    def estimate_step_length(a_pp_step):
        res = []
        for e in a_pp_step:
            if e < A_STEP:
                res.append(GAMMA_FOURTH_ROOT *
                           np.power(e, 1/4) + DELTA_FOURTH_ROOT)
            else:
                res.append(GAMMA_LOG * np.log(e) + DELTA_LOG)
        return res

    step_length = estimate_step_length(a_pp_step)

    # 位置の推定
    position = INITIAL_POSITION
    pos_history = []
    h_k = h_t_res[sorted_t_peak]
    for i, sl in enumerate(step_length):
        position = [position[0] + sl * np.sin(h_k[i]), position[1] + sl * np.cos(h_k[i])]
        pos_history.append(position)

    pos_history = np.array(pos_history)


    def show_pos_history():
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')

        pos_x = pos_history[:, 0]
        pos_y = pos_history[:, 1]

        plt.scatter(pos_x, pos_y, marker="o")
        plt.show()


    show_pos_history()


if __name__ == "__main__":
    main()
