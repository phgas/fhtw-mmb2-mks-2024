import numpy as np
import matplotlib.pyplot as plt
import math


def setup_simulation(simulation_time, step_time):
    number_of_steps = int(simulation_time / step_time)
    time_array = np.zeros(number_of_steps)
    x = np.zeros(number_of_steps)          # Displacement in meter
    x_dot = np.zeros(number_of_steps)      # Velocity in meter/second
    alpha = np.zeros(number_of_steps)      # Angular displacement in radians
    alpha_dot = np.zeros(number_of_steps)  # Angular velocity in radians/second
    return number_of_steps, time_array, x, x_dot, alpha, alpha_dot


def setup_inital_conditions(initial_x, initial_x_dot, initial_alpha, initial_alpha_dot):
    x[0] = initial_x
    x_dot[0] = initial_x_dot
    alpha[0] = initial_alpha
    alpha_dot[0] = initial_alpha_dot


def run_simulation():
    for i in range(0, number_of_steps-1):

        x[i + 1] = x_dot[i] * STEP_TIME + x[i]
        # is it even correct to only use alpha_dot in second place??
        x_dot[i + 1] = (x[i] * (alpha_dot[i]**2) + G * np.sin(PHI_0_RAD) * np.cos(alpha[i]) + ((C * x[i] / M_GONDEL) * 
                       (((LS**2 + R_KARUSELL**2) / np.sqrt(LS**4 + (LS**2 * R_KARUSELL) + (x[i]**2 * LS**2) + x[i]**2 * R_KARUSELL**2)) 
                       - 1))) * STEP_TIME + x_dot[i]

        # is it even correct to only use alpha_dot and x_dot in second place??
        alpha[i + 1] = alpha_dot[i] * STEP_TIME + alpha[i]

        alpha_dot[i + 1] = (-(2 * alpha_dot[i] * x_dot[i] + G * np.sin(PHI_0_RAD) * np.sin(alpha[i]) * x[i]) /
                             (x[i]**2 + (5/3) * (B_GONDEL**2 + H_GONDEL**2) + 20 * R_KARUSELL**2)) * STEP_TIME + alpha_dot[i]

        time_array[i + 1] = (i + 1) * STEP_TIME


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
    C = 60000                           # Spring stiffness [N/m]
    M_GONDEL = 300                      # Mass of gondola [kg]
    N_GONDEL = 0.2333                   # Rotational speed of gondola [radians/second]
    B_GONDEL = 1.5                      # Width of gondola [m]
    H_GONDEL = 1.5                      # Height of gondola [m]
    LS = R_KARUSELL * np.tan(PHI_0_RAD)

    SIMULATION_TIME = 10  # [s]
    STEP_TIME = 0.0001     # [s]
    number_of_steps, time_array, x, x_dot, alpha, alpha_dot = setup_simulation(
        SIMULATION_TIME, STEP_TIME)
    setup_inital_conditions(initial_x=R_KARUSELL,               # Initial displacement (6 meter)
                            initial_x_dot=0,                    # Initial velocity ()
                            initial_alpha=0,                    # Initial angular displacement ()
                            initial_alpha_dot=2*np.pi*N_GONDEL) # Initial angular velocity ()
    run_simulation()
    plot_results()

    local_minima, local_maxima = get_local_extremes(x)
    get_global_maxima(local_minima, local_maxima)