import matplotlib.pyplot as plt
from math import exp, sqrt, pi, log10

MAXVELL_DISTRIBUTION_TRASH_HOLD = 0.01
INTERPOLATION_SCALE = 4

#arr is an input data (pixel cords), range_ is min and max values on axis, and max_ is picture size
def to_linear_form(arr : list[int], range_ : list[float], max_ : int):
    return [10 ** (range_[0] + float(i) / float(max_) * (range_[1] - range_[0])) for i in arr]

#need to supress log10(0)
def temp_log(x: float):
    if x <= 0: return None
    return log10(x)

#arr is set of points in linear cords system
def to_logarithmic_form(arr: list[float | int] | range):
    return [temp_log(i) for i in arr]

#linear interpolation for array of points, just plugin more points between each
def linear_interpolate(arr: list[float], multiplier: int):
    ret = [arr[0]]
    for i in range(len(arr) - 1):
        ret += [arr[i]]
        step = (arr[i + 1] - arr[i]) / (multiplier + 2)
        temp = [arr[i] + step * (j + 1) for j in range(multiplier + 1)]
        ret += temp
    ret += [arr[-1]]
    return ret

# maxvell distribution function
def maksvell_dist(t: float, e: float):
    if t == 0:
        return 0
    f_p = (2 * pi) / (sqrt((pi * t) ** 3))
    s_p = sqrt(e) * exp(-e / t)
    return f_p * s_p


# finding closer energy in input data
def find_closest_energy(e: float, initial_distribution: list[list[float]]):
    if e < initial_distribution[0][0]:
        return 0
    if e > initial_distribution[0][-1]:
        return initial_distribution[1][-1]
    for i in range(len(initial_distribution[0]) - 1):
        if initial_distribution[0][i] < e < initial_distribution[0][i+1]:
            t = (e - initial_distribution[0][i])/(initial_distribution[0][i+1] - initial_distribution[0][i])
            return initial_distribution[1][i]+(initial_distribution[1][i+1] - initial_distribution[1][i])*t

#running through energy range and integrating maxvell distribution multiplied by initial distribution value for de
def integrate_distribution(energy_range_: range, initial_distribution: list[list[float]]):
    integrated_distributions = []

    for i in energy_range_:
        print(str(i/energy_range_[-1]*100)[:5]+"%")
        integral_sum = 0
        check_sum = 0
        j = 0
        while 1 - check_sum > MAXVELL_DISTRIBUTION_TRASH_HOLD:
            current_dist = (maksvell_dist(i, j))
            check_sum += current_dist
            integral_sum += current_dist * find_closest_energy(j, initial_distribution)
            j += 1

        integrated_distributions.append(integral_sum)
    return integrated_distributions

if __name__ == "__main__":
    # input data (pixel coordinates)
    file = [list(map(int, i.split(';'))) for i in open("input_data.txt").read().split("\n")]

    picture_sizeX = 447
    picture_sizeY = 363

    # [min, max] 10 powers on axis, eg x belongs to [10^2, 10^5]
    x_range = [2, 5]
    y_range = [-4, -1]

    # separated pixel coordinates
    x_cords = [i[0] for i in file]  # e-v
    y_cords = [picture_sizeY - i[1] for i in file]  # Y

    # converting pixel coordinates to linear coordinate system set of points
    x_cords = to_linear_form(x_cords, x_range, picture_sizeX)
    y_cords = to_linear_form(y_cords, y_range, picture_sizeY)


    energy_range = range(10, 10000, 50)
    int_dist = integrate_distribution(energy_range, [x_cords, y_cords])
    plt.plot(to_logarithmic_form(x_cords), to_logarithmic_form(y_cords))
    plt.plot(to_logarithmic_form(energy_range), to_logarithmic_form(int_dist))

    plt.xlabel('E, эВ')
    plt.ylabel('Y, атомов/ион')

    plt.show()
