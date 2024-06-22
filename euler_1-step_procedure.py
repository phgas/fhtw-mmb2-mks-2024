import numpy as np
import matplotlib.pyplot as plt


def setup_simulation(simulation_time, step_time):
    number_of_steps = int(simulation_time / step_time)
    time_array = np.zeros(number_of_steps)
    phi = np.zeros(number_of_steps)   # angular displacement in radians
    omega = np.zeros(number_of_steps)  # angular velocity in radians/second
    return number_of_steps, time_array, phi, omega


def setup_inital_conditions(initial_phi, initial_omega):
    phi[0] = initial_phi
    omega[0] = initial_omega
    return phi, omega


def run_simulation():
    for i in range(0, number_of_steps-1):
        phi[i + 1] = omega[i] * STEP_TIME + phi[i]
        omega[i + 1] = (-G / L * np.sin(phi[i])) * STEP_TIME + omega[i]
        time_array[i + 1] = (i + 1) * STEP_TIME
    return time_array, phi, omega


def plot_results() -> None:
    plt.plot(time_array, phi, label='$\phi(t)$')
    plt.plot(time_array, omega, label='$\dot{\phi}(t)$')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


if __name__ == "__main__":
    L = 1                 # Length of pendulum [m]
    G = 9.81              # Gravity Acceleration  [m/s^2]
    SIMULATION_TIME = 20  # [s]
    STEP_TIME = 0.001     # [s]
    number_of_steps, time_array, phi, omega = setup_simulation(
        SIMULATION_TIME, STEP_TIME)

    setup_inital_conditions(initial_phi=np.pi/2, # Initial angular velocity (at rest)
                            initial_omega=0)     # Initial angular displacement (90 degrees, horizontal position)

    run_simulation()
    plot_results()
