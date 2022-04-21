from params import *


def extract_t1(t_step, a_step):
    t1 = set()
    for t in range(len(t_step)):
        if (t - SAMPLE_N//2 < 0) or (t + SAMPLE_N//2 >= len(t_step)):
            continue
        if a_step[t] <= A_PEAK:
            continue

        is_ok = True
        for i in range(-SAMPLE_N//2, SAMPLE_N//2 + 1):
            if i == 0:
                continue
            if a_step[t] <= a_step[t+i]:
                is_ok = False
        if is_ok:
            t1.add(t)
    return t1


def extract_t2(t_step, a_step):
    t2 = set()
    for t in range(len(t_step)):
        if t + SAMPLE_N//2 >= len(t_step):
            continue
        max_diff_left = -10000000
        max_diff_right = -10000000
        for i in range(1, SAMPLE_N//2 + 1):
            max_diff_left = max(max_diff_left, a_step[t] - a_step[t-i])
            max_diff_right = max(max_diff_right, a_step[t] - a_step[t+i])

        if max_diff_left > A_PP and max_diff_right > A_PP:
            t2.add(t)
    return t2


def extract_t3(t_step, a_step):
    t3 = set()
    for t in range(len(t_step)):
        if (t - SAMPLE_N//2 < 0) or (t + SAMPLE_N//2 >= len(t_step)):
            continue
        sum_diff_left = 0
        sum_diff_right = 0
        for i in range(t - SAMPLE_N//2, t + SAMPLE_N//2 + 1):
            if i <= t - 1:
                sum_diff_left += (a_step[i+1] - a_step[i])
            elif i >= t + 1:
                sum_diff_right += (a_step[i] - a_step[i-1])

        if sum_diff_left > 0 and sum_diff_right < 0:
            t3.add(t)
    return t3


def extract_t_peak(t_step, a_step):
    t1 = extract_t1(t_step, a_step)
    t2 = extract_t2(t_step, a_step)
    t3 = extract_t3(t_step, a_step)
    t_peak = t1 & t2 & t3
    return t_peak
