# importing libraries
import numpy as np
import csv
import os
from config import *
from solver import initialize_fields, step

# function to calculate the wall-averaged Nu & maximum mid-plane velocities
def calculate_validation_metrics(T, u, v):
    # computing the wall-averaged Nusselt Number
    Nu_sum = 0.0
    T_wall = 1.0
    for j in range(1, ny + 1):
        Nu_local = (T_wall - T[1, j]) / (dx / 2.0)
        Nu_sum += Nu_local
    Nu_avg = Nu_sum / ny
    
    # computing maximum mid-plane velocities
    u_max = 0.0
    for j in range(1, ny + 1):
        val = abs(u[nx // 2, j])
        if val > u_max:
            u_max = val
            
    v_max = 0.0
    for i in range(1, nx + 1):
        val = abs(v[i, ny // 2])
        if val > v_max:
            v_max = val
            
    return Nu_avg, u_max, v_max

# function to store the calculated wall-averaged Nu & maximum mid-plane velocities for simulations
def export_validation_metrics(Ra_val, Nu_avg, u_max, v_max):
    # appending calculated values to csv file
    file_exists = os.path.isfile('Results/validation_metrics.csv')
    with open('Results/validation_metrics.csv', mode='a', newline='') as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(['Rayleigh Number', 'Nu_avg', 'Max U (mid-plane)', 'Max V (mid-plane)'])
        writer.writerow([Ra_val, round(Nu_avg, 4), round(u_max, 4), round(v_max, 4)])

# function to store the final grid variables
def export_field_data(Ra_val, T, u, v):
    # creating a unique csv file for each Ra run for contour plotting
    filename = f'Results/results_Ra_{int(Ra_val)}.csv'
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['x', 'y', 'T', 'u', 'v'])
        
        # interpolating staggered grid variables to cell centers for field extraction
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                x_c = (i - 0.5) * dx
                y_c = (j - 0.5) * dy
                T_c = T[i, j]
                
                u_c = 0.5 * (u[i-1, j] + u[i, j])
                v_c = 0.5 * (v[i, j-1] + v[i, j])
                
                writer.writerow([round(x_c, 5), round(y_c, 5), round(T_c, 5), round(u_c, 5), round(v_c, 5)])

# function to run the simulation
def run_simulation():
    print(f"Starting FVM SIMPLE Solver Run | Ra = {Ra}, Grid = {nx}x{ny}")
    u, v, p, T = initialize_fields()
    
    for iteration in range(1, max_iters + 1):
        err = step(u, v, p, T)
        
        if iteration % 100 == 0:
            # calculating the wall-averaged Nu & maximum mid-plane velocities
            Nu_avg, u_max, v_max = calculate_validation_metrics(T, u, v)

            print(f"Iter: {iteration:<6} | dT Err: {err:.2e} | Nu_avg: {Nu_avg:.4f} | U_max: {u_max:.4f} | V_max: {v_max:.4f}")
            
        if err < tol:
            print(f"steady-state convergence reached at iteration {iteration}.")
            break
            
    # calculating the wall-averaged Nu & maximum mid-plane velocities
    Nu_avg, u_max, v_max = calculate_validation_metrics(T, u, v)
    
    print("\n--- Simulation Results ---")
    print(f"Average Wall Nusselt Number (Nu_avg): {Nu_avg:.4f}")
    print(f"Max U-Velocity (mid-plane): {u_max:.4f}")
    print(f"Max V-Velocity (mid-plane): {v_max:.4f}")

    print("\n--- Exporting Data ---")
    export_validation_metrics(Ra, Nu_avg, u_max, v_max)
    export_field_data(Ra, T, u, v)
    print(f"Parameters appended to 'validation_metrics.csv'.")
    print(f"Field data saved to 'results_Ra_{int(Ra)}.csv'.")

# executes the entire script
if __name__ == "__main__":
    run_simulation()