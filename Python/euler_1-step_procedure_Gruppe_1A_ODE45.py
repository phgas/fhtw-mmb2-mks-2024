import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import math


def system_dynamics(t, y):
    x, x_dot, alpha, alpha_dot = y
    dxdt = x_dot
    dx_dotdt = (x * (alpha_dot**2) + G * np.sin(PHI_0_RAD) * np.cos(alpha) + ((C * x / M_GONDEL) * 
                (((LS**2 + R_KARUSELL**2) / np.sqrt(LS**4 + (LS**2 * R_KARUSELL) + (x**2 * LS**2) + x**2 * R_KARUSELL**2)) 
                - 1)))
    dalphadt = alpha_dot
    dalpha_dotdt = (-(2 * alpha_dot * x_dot + G * np.sin(PHI_0_RAD) * np.sin(alpha) * x) /
                    (x**2 + (5/3) * (B_GONDEL**2 + H_GONDEL**2) + 20 * R_KARUSELL**2))
    return [dxdt, dx_dotdt, dalphadt, dalpha_dotdt]


def run_simulation():
    t_span = [0, SIMULATION_TIME]
    t_eval = np.arange(0, SIMULATION_TIME, STEP_TIME)
    solution = solve_ivp(system_dynamics, t_span, initial_conditions, t_eval=t_eval, method='RK45')
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


def get_global_maxima(local_minima, local_maxima) -> None:
    global_minima_above_zero = min(local_maxima)
    global_minima_below_zero = min(local_minima)
    global_maxima_above_zero = max(local_maxima)
    global_maxima_below_zero = max(local_minima)

    print(f"Global maxima: {global_minima_above_zero} to {global_maxima_above_zero} | Difference is {global_maxima_above_zero-global_minima_above_zero}")
    print(f"Global minima: {global_minima_below_zero} to {global_maxima_below_zero} | Difference is {global_maxima_below_zero-global_minima_below_zero}")






if __name__ == "__main__":
    G = 9.81                            # Gravity acceleration [m/s^2]
    PHI_0_DEG = 30                      # Angle of tilted carousel [°]
    PHI_0_RAD = math.radians(PHI_0_DEG) # Angle of tilted carousel [radians]
    R_KARUSELL = 6                      # Radius of carousel [m]
    C = 100000                           # Spring stiffness [N/m]
    M_GONDEL = 300                      # Mass of gondola [kg]
    N_GONDEL = 0.2333                   # Rotational speed of gondola [radians/second]
    B_GONDEL = 1.5                      # Width of gondola [m]
    H_GONDEL = 1.5                      # Height of gondola [m]
    LS = R_KARUSELL * np.tan(PHI_0_RAD)

    SIMULATION_TIME = 10  # [s]
    STEP_TIME = 0.001     # [s]

    initial_conditions = [R_KARUSELL, 0, 0, 2 * np.pi * N_GONDEL]
    time_array, x, x_dot, alpha, alpha_dot = run_simulation()
    local_minima, local_maxima = get_local_extremes(x)

    plot_results()
    get_global_maxima(local_minima, local_maxima)
    