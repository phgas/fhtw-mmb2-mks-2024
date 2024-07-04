import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import math
from termcolor import colored
import colorama

colorama.init()

def system_dynamics(t, y):
    x, x_dot, alpha, alpha_dot = y
    dxdt = x_dot
    dx_dotdt = (x * (alpha_dot**2) + G * np.sin(PHI_0_RAD) *
                np.cos(alpha) + ((C / M_GONDEL) * (R_KARUSELL - x)))
    dalphadt = alpha_dot
    dalpha_dotdt = (-(2 * alpha_dot * x_dot * x + G * np.sin(PHI_0_RAD) * np.sin(alpha) * x) /
                    (x**2 + (5/3) * (B_GONDEL**2 + H_GONDEL**2) + 20 * R_KARUSELL**2))
    return [dxdt, dx_dotdt, dalphadt, dalpha_dotdt]


def run_simulation():
    t_span = [0, SIMULATION_TIME]
    solution = solve_ivp(system_dynamics, t_span,
                         initial_conditions, method='RK45', rtol=1e-3, atol=1e-6)

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


def plot_results(time_array, x, x_dot, alpha, alpha_dot) -> None:
    fig, axs = plt.subplots(4, 1, figsize=(10, 12))

    axs[0].plot(time_array, x, label=r'$x(t)$')
    axs[0].set_ylabel(r'$x$ [m]')
    axs[0].legend(loc='upper right')
    axs[0].grid()
    axs[0].set_xlim(min(time_array), max(time_array))

    axs[1].plot(time_array, x_dot, label=r'$\dot{x}(t)$')
    axs[1].set_ylabel(r'$\dot{x}$ [m/s]')
    axs[1].legend(loc='upper right')
    axs[1].grid()
    axs[1].set_xlim(min(time_array), max(time_array))

    degrees_alpha = np.degrees(alpha)
    axs[2].plot(time_array, degrees_alpha, label=r'$\alpha(t)$')
    axs[2].set_ylabel(r'$\alpha$ [deg]')
    axs[2].legend(loc='upper right')
    axs[2].grid()
    axs[2].set_xlim(min(time_array), max(time_array))

    degrees_alpha_dot = np.degrees(alpha_dot)
    axs[3].plot(time_array, degrees_alpha_dot, label=r'$\dot{\alpha}(t)$')
    axs[3].set_ylabel(r'$\dot{\alpha}$ [deg/s]')
    axs[3].legend(loc='upper right')
    axs[3].grid()
    axs[3].set_xlim(min(time_array), max(time_array))

    for i in range(1, int(max(degrees_alpha) / 360) + 1):
        idx = np.where(np.isclose(degrees_alpha, 360 * i, atol=1))[0]
        if len(idx) > 0:
            axs[2].plot(time_array[idx[0]], degrees_alpha[idx[0]], 'ro')
            for ax in axs:
                ax.axvline(x=time_array[idx[0]],
                           color='r', linestyle='--', linewidth=1)

    y_ticks = np.arange(0, max(degrees_alpha) + 360, 360)
    axs[2].set_yticks(y_ticks)
    axs[2].set_yticklabels([f'{int(y)}' for y in y_ticks])

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.text(0.5, 0.01, 'Time (t)', ha='center', fontsize=12)
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


def get_global_extremes(local_minima, local_maxima):
    try:
        global_minima = round(min(local_minima), 5)
        global_maxima = round(max(local_maxima), 5)
        difference = round(global_maxima - global_minima, 5)
        return global_minima, global_maxima, difference

    except ValueError:

        return None, None, None


def find_optimal_C(start, step_size, MAX_RADIAL_DISPLACEMENT, amount_of_results):
    global C
    results_counter = 0

    for C in np.arange(start, 999_999_999_999, step_size):
        try:
            time_array, x, x_dot, alpha, alpha_dot = run_simulation()
            local_minima, local_maxima = get_local_extremes(x)
            global_minima, global_maxima, difference = get_global_extremes(
                local_minima, local_maxima)
            if difference is not None and difference <= MAX_RADIAL_DISPLACEMENT:
                print(colored(
                    f"[C: {C:<8} N/m] Global maxima: {global_maxima:<8} | Global minima: {global_minima:<8} | Difference is {difference:<8}", 'green'))
                results_counter += 1
                if results_counter >= amount_of_results:
                    break
            elif difference is None:
                print(colored(f"[C: {C:<8} N/m] no Global minima or maxima found!", 'red'))
            else:
                print(colored(
                    f"[C: {C:<8} N/m] Global maxima: {global_maxima:<8} | Global minima: {global_minima:<8} | Difference is {difference:<8}", 'red'))
                results_counter = 0

        except Exception as e:
            print(colored(f"[C: {C:<8} N/m] Error: {str(e)}", 'red'))
            results_counter = 0


if __name__ == "__main__":
    G = 9.81                            # Gravity acceleration [m/s^2]
    PHI_0_DEG = 30                      # Angle of tilted carousel [°]
    PHI_0_RAD = math.radians(PHI_0_DEG) # Angle of tilted carousel [radians]
    R_KARUSELL = 6                      # Radius of carousel [m]
    # C = 150                            # Spring stiffness [N/m]
    M_GONDEL = 300                      # Mass of gondola [kg]
    # Rotational speed of gondola [radians/second]
    N_GONDEL = 0.2333
    B_GONDEL = 1.5                      # Width of gondola [m]
    H_GONDEL = 1.5                      # Height of gondola [m]
    LS = R_KARUSELL * np.tan(PHI_0_RAD)

    SIMULATION_TIME = 20                 # [s]
    MAX_RADIAL_DISPLACEMENT = 0.005     # Maximum radial displacement [m]

    initial_conditions = [R_KARUSELL, 0, 0, 2 * np.pi * N_GONDEL]


    find_optimal_C(start=0, step_size=100_000,
                   MAX_RADIAL_DISPLACEMENT=MAX_RADIAL_DISPLACEMENT, amount_of_results=2)
