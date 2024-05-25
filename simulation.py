import numpy as np
import matplotlib.pyplot as plt

def simulate(angle_deg, initial_velocity, plot_trajectory=False):
    # Constants
    C_d = 0.47  # Drag coefficient for a sphere
    rho = 1.225  # Air density in kg/m^3
    r = 0.02  # Radius of ping pong ball in meters (approx. 2 cm)
    A = np.pi * r**2  # Cross-sectional area
    m = 0.0027  # Mass of ping pong ball in kg (approx. 2.7 grams)
    g = 9.81  # Acceleration due to gravity in m/s^2

    # Convert angle to radians
    theta = np.radians(angle_deg)

    # Initial conditions
    v_x = initial_velocity * np.cos(theta)
    v_y = initial_velocity * np.sin(theta)
    x, y = 0, 0  # Initial position

    # Time step
    dt = 0.01  # Time step in seconds

    # Lists to store the trajectory
    x_vals = [x]
    y_vals = [y]

    # Simulation loop
    while y >= 0:
        v = np.sqrt(v_x**2 + v_y**2)
        F_d = 0.5 * C_d * rho * A * v**2
        F_d_x = F_d * (v_x / v)
        F_d_y = F_d * (v_y / v)

        # Update velocities
        v_x -= (F_d_x / m) * dt
        v_y -= (g + F_d_y / m) * dt

        # Update positions
        x += v_x * dt
        y += v_y * dt

        # Store values
        x_vals.append(x)
        y_vals.append(y)

    if plot_trajectory:
        plt.figure(figsize=(10, 5))
        plt.plot(x_vals, y_vals, label=f'Angle: {angle_deg}°, Velocity: {initial_velocity} m/s')
        plt.ylim(-0.1, 1)
        plt.xlabel('Distance (m)')
        plt.ylabel('Height (m)')
        plt.title('Ping Pong Ball Trajectory')
        plt.legend()
        plt.grid(True)
        plt.show()
    return x

def guessLaunchAngle(distance, initial_velocity):
    # Starting values
    angle_deg = 0
    max_angle = 90
    angle_step = 0.01
    distance_travelled = simulate(angle_deg, initial_velocity)
    i = 0

    list_solution = []

    # Loop until the distance is close enough
    while abs(distance_travelled - distance) > 0.01:

        angle_deg += angle_step
        
        distance_travelled = simulate(angle_deg, initial_velocity)

        if (abs(distance_travelled - distance) < 0.01):
            list_solution.append((angle_deg, abs(distance_travelled - distance)))

        #print(angle_deg)
        #print(distance_travelled)

        if angle_deg > max_angle or angle_deg < 0 or i > 100000:
            # print("Error: Could not find a suitable launch angle.")
            break
        i += 1

    if (len(list_solution) == 0):
        print("Error: Could not find a suitable launch angle.")
        return -1
    else:
        sorted_list = sorted(list_solution, key=lambda x: x[1])
        angle_deg = sorted_list[0][0]
    print(f"Launch angle: {angle_deg:.2f}")
    return angle_deg

# Guess the launch angle required to hit a target distance
initial_velocity = 8  # Set a constant initial velocity
target_distance = 2  # Target distance to hit
launch_angle = guessLaunchAngle(target_distance, initial_velocity)
if launch_angle == -1:
    exit()
distance = simulate(angle_deg=launch_angle, initial_velocity=initial_velocity, plot_trajectory=True)
print(f"Distance traveled: {distance:.2f} meters")
