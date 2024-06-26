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
        x_dot[i + 1] = (x[i] * (alpha_dot[i]**2) - ((C / M_GONDEL) * (x[i] + DELTA_L)) - G * np.sin(PHI_0_RAD) * np.sin(alpha[i]))* STEP_TIME + x_dot[i]
        
        alpha[i + 1] = alpha_dot[i] * STEP_TIME + alpha[i]
        alpha_dot[i + 1] = (-(M_GONDEL * G * np.sin(PHI_0_RAD) * np.cos(alpha[i]) * x[i]) / 
                             (M_GONDEL * x[i]**2 + I_STANGE + 20*I_GONDEL)) * STEP_TIME + alpha_dot[i]

        time_array[i + 1] = (i + 1) * STEP_TIME


def plot_results() -> None:
    plt.plot(time_array, x, label=r'$x(t)$')
    plt.plot(time_array, x_dot, label=r'$\dot{x}(t)$')
    plt.plot(time_array, alpha, label=r'$\alpha(t)$')
    plt.plot(time_array, alpha_dot, label=r'$\dot{\alpha}(t)$')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


if __name__ == "__main__":
    G = 9.81                               # Gravity acceleration [m/s^2]
    PHI_0_DEG = 30                         # Angle of tilted carousel [°]
    PHI_0_RAD = math.radians(PHI_0_DEG)    # Angle of tilted carousel [radians]
    C = 107.22                             # Spring stiffness [N/m]
    M_GONDEL = 300                         # Mass of gondola [kg]
    N_GONDEL = 0.2333                      # Mass of gondola [radians/second]
    DELTA_L = 0.005                        #  [m]
    I_STANGE = 500                         # Moment of inertia of rod [kg*m^2]
    I_GONDEL = 218250                      # Moment of inertia of gondola [kg*m^2]

    SIMULATION_TIME = 20  # [s]
    STEP_TIME = 0.001     # [s]
    number_of_steps, time_array, x, x_dot, alpha, alpha_dot = setup_simulation(
        SIMULATION_TIME, STEP_TIME)
    setup_inital_conditions(initial_x=6,                 # Initial displacement (6 meter)
                            initial_x_dot=0,              # Initial velocity (0 degrees, horizontal position)
                            initial_alpha=0,             # Initial angular displacement (at rest)
                            initial_alpha_dot=2*np.pi*N_GONDEL)   # Initial angular velocity (90 degrees, horizontal position)
    run_simulation()
    plot_results()
