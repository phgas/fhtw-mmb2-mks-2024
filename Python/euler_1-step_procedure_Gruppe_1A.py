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
        x_dot[i + 1] = (x[i] * (alpha_dot[i]**2) + G * np.sin(PHI_0_RAD) * 
                        np.cos(alpha[i]) + ((C / M_GONDEL) * 
                            (R_KARUSELL - x[i]))) * STEP_TIME + x_dot[i]
        alpha[i + 1] = alpha_dot[i] * STEP_TIME + alpha[i]
        alpha_dot[i + 1] = (-(2 * alpha_dot[i] * x_dot[i] * x[i] + G * 
                              np.sin(PHI_0_RAD) * np.sin(alpha[i]) * x[i]) /
                                (x[i]**2 + (5/3) * (B_GONDEL**2 + H_GONDEL**2) + 
                                    20 * R_KARUSELL**2)) * STEP_TIME + alpha_dot[i]
        time_array[i + 1] = (i + 1) * STEP_TIME

def plot_results() -> None:
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

if __name__ == "__main__":
    G = 9.81                            # Gravity acceleration [m/s^2]
    PHI_0_DEG = 30                      # Angle of tilted carousel [°]
    PHI_0_RAD = math.radians(PHI_0_DEG) # Angle of tilted carousel [radians]
    R_KARUSELL = 6                      # Radius of carousel [m]
    C = 3000000                         # Spring stiffness [N/m]
    M_GONDEL = 300                      # Mass of gondola [kg]
    N_GONDEL = 0.2333                   # Rotational speed [radians/second]
    B_GONDEL = 1.5                      # Width of gondola [m]
    H_GONDEL = 1.5                      # Height of gondola [m]

    SIMULATION_TIME = 10  # [s]
    STEP_TIME = 0.000001     # [s]
    number_of_steps, time_array, x, x_dot, alpha, alpha_dot = setup_simulation(
        SIMULATION_TIME, STEP_TIME)
    setup_inital_conditions(initial_x=R_KARUSELL, # Initial displacement (meter)
                            initial_x_dot=0,      # Initial velocity (meter/second)
                            initial_alpha=0,      # Initial angular displacement (radians)
                            initial_alpha_dot=2*np.pi*N_GONDEL) # Initial angular velocity (radians/second)
    run_simulation()
    plot_results()