import matplotlib.pyplot as plt
from detect_step import extract_t_peak, extract_t1, extract_t2, extract_t3


def show_magdata(magdata):
    plt.plot(magdata[:, 1])
    plt.plot(magdata[:, 2])
    plt.plot(magdata[:, 3])
    plt.show()


def show_t_peaks(t_step, a_step):
    t_peak = extract_t_peak(t_step, a_step)
    t1 = extract_t1(t_step, a_step)
    t2 = extract_t2(t_step, a_step)
    t3 = extract_t3(t_step, a_step)

    fig = plt.figure(figsize=[20, 10])
    plt.title(
        f"t_peak: {len(t_peak)}, t1: {len(t1)}, t2: {len(t2)}, t3: {len(t3)}"
    )

    ax1 = fig.add_subplot(4, 1, 1)
    ax1.plot(t_step, a_step, marker="o")
    ax1.plot(t_step[list(t_peak)], a_step[list(t_peak)],
             linestyle="", marker="o")

    ax2 = fig.add_subplot(4, 1, 2)
    ax2.plot(t_step, a_step, marker="o")
    ax2.plot(t_step[list(t1)], a_step[list(t1)], linestyle="", marker="o")

    ax3 = fig.add_subplot(4, 1, 3)
    ax3.plot(t_step, a_step, marker="o")
    ax3.plot(t_step[list(t2)], a_step[list(t2)], linestyle="", marker="o")

    ax4 = fig.add_subplot(4, 1, 4)
    ax4.plot(t_step, a_step, marker="o")
    ax4.plot(t_step[list(t3)], a_step[list(t3)], linestyle="", marker="o")

    plt.show()


def show_t_peaks_and_valleys(t_step, a_step, t_head):
    t_peak = extract_t_peak(t_step, a_step)
    plt.title(
        f"t_peak: {len(t_peak)}, t_head: {len(t_head)}"
    )

    plt.plot(t_step, a_step, marker="o")
    plt.plot(t_step[list(t_peak)], a_step[list(t_peak)],
             linestyle="", marker="o")
    plt.plot(t_step[list(t_head)], a_step[list(t_head)],
             linestyle="", marker="o")
    plt.show()
