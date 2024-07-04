import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import math
from scipy.signal import find_peaks


def system_dynamics(t, y):
    x, x_dot, alpha, alpha_dot = y
    dxdt = x_dot
    dx_dotdt = (x * (alpha_dot**2) + G * np.sin(PHI_0_RAD) * np.cos(alpha) + (((C * M_GONDEL) * (R_KARUSELL - x))))
    dalphadt = alpha_dot
    dalpha_dotdt = (-(2 * alpha_dot * x_dot * x + G * np.sin(PHI_0_RAD) * np.sin(alpha) * x) /
                    (x**2 + (5/3) * (B_GONDEL**2 + H_GONDEL**2) + 20 * R_KARUSELL**2))
    return [dxdt, dx_dotdt, dalphadt, dalpha_dotdt]


def run_simulation():
    t_span = [0, SIMULATION_TIME]
    solution = solve_ivp(system_dynamics, t_span, initial_conditions, method='RK45', rtol=1e-3, atol=1e-6) 

    """ Solver keeps local error estimates less than `atol + rtol * abs(y) `
        - rtol controls relative accuracy (number of correct digits)
        - atol controls absolute accuracy (number of correct decimal places)

        To achieve desired `rtol`, set `atol` smaller than the smallest value of `rtol * abs(y)`
        Conversely, to achieve desired `atol`, set `rtol` smaller than smallest value of `rtol * abs(y)`

        Default values used for rtol and atol.
    """
    time_array = solution.t
    x, x_dot, alpha, alpha_dot = solution.y
    return time_array, x, x_dot, alpha, alpha_dot



def plot_results() -> None:
    fig, axs = plt.subplots(4, 1, figsize=(10, 12))

    axs[0].plot(time_array, x, label=r'$x(t)$')
    axs[0].set_ylabel(r'$x$ [m]')
    axs[0].legend(loc='upper right')
    axs[0].grid()

    axs[1].plot(time_array, x_dot, label=r'$\dot{x}(t)$')
    axs[1].set_ylabel(r'$\dot{x}$ [m/s]')
    axs[1].legend(loc='upper right')
    axs[1].grid()

    axs[2].plot(time_array, alpha, label=r'$\alpha(t)$')
    axs[2].set_ylabel(r'$\alpha$ [rad]')
    axs[2].legend(loc='upper right')
    axs[2].grid()

    axs[3].plot(time_array, alpha_dot, label=r'$\dot{\alpha}(t)$')
    axs[3].set_ylabel(r'$\dot{\alpha}$ [rad/s]')
    axs[3].legend(loc='upper right')
    axs[3].grid()

    fig.text(0.5, 0.04, 'Time (t)', ha='center', fontsize=12)
    fig.suptitle(f"{C=} N/m", fontsize=16)
    plt.show()


def get_local_extremes(data_series):
    local_maxima = []
    local_minima = []
    for i in range(1, len(data_series) - 1):
        if data_series[i] > data_series[i - 1] and data_series[i] > data_series[i + 1]:
            local_maxima.append(data_series[i])
        if data_series[i] < data_series[i - 1] and data_series[i] < data_series[i + 1]:
            local_minima.append(data_series[i])
    return local_minima, local_maxima


def get_global_maxima(local_minima, local_maxima):
    try:
        global_minima_below_zero = round(min(local_minima), 5)
        global_maxima_above_zero = round(max(local_maxima), 5)
        difference = round(global_maxima_above_zero - global_minima_below_zero, 5)

        print(f"[C: {C:<4} N/m] Global maxima: {global_maxima_above_zero:<8} | Global minima: {global_minima_below_zero:<8} | Difference is {difference:<8}")
    except ValueError:
        print(f"[C: {C:<3} N/m] No Global minima or maxima found!")


def find_optimal_C(start, end, step_size):
    for C in np.arange(start, end, step_size):  # Iterate over a range of C values
        time_array, x, x_dot, alpha, alpha_dot = run_simulation(C)
        local_minima, local_maxima = get_local_extremes(x)
        difference_above, difference_below = get_global_maxima(local_minima, local_maxima)
        if difference_above <= MAX_RADIAL_DISPLACEMENT and difference_below <= MAX_RADIAL_DISPLACEMENT:
            return C, difference_above, difference_below
    return None, None, None


if __name__ == "__main__":
    G = 9.81                            # Gravity acceleration [m/s^2]
    PHI_0_DEG = 30                      # Angle of tilted carousel [°]
    PHI_0_RAD = math.radians(PHI_0_DEG) # Angle of tilted carousel [radians]
    R_KARUSELL = 6                      # Radius of carousel [m]
    C = 150                          # Spring stiffness [N/m]
    M_GONDEL = 300                      # Mass of gondola [kg]
    N_GONDEL = 0.2333                   # Rotational speed of gondola [radians/second]
    B_GONDEL = 1.5                      # Width of gondola [m]
    H_GONDEL = 1.5                      # Height of gondola [m]
    LS = R_KARUSELL * np.tan(PHI_0_RAD)

    SIMULATION_TIME = 20  # [s]
    MAX_RADIAL_DISPLACEMENT = 0.005      # Maximum radial displacement [m]

    initial_conditions = [R_KARUSELL, 0, 0, 2 * np.pi * N_GONDEL]

    for C in np.arange(31.1, 31.2, 0.01):
        time_array, x, x_dot, alpha, alpha_dot = run_simulation()
        #plot_results()
        local_minima, local_maxima = get_local_extremes(x)
        get_global_maxima(local_minima, local_maxima)