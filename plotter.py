# importing libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# function to plot the results
def plot_results(ra_val):
    # specifying result file name
    filename = f'Results/results_Ra_{int(ra_val)}.csv'

    #loading data file
    print(f"Loading data from {filename}...")
    df = pd.read_csv(filename)

    # using pivot to reshape the 1D CSV columns into 2D spatial grids (rows = y coordinates, columns = x coordinates)
    T_grid = df.pivot(index='y', columns='x', values='T')
    U_grid = df.pivot(index='y', columns='x', values='u')
    V_grid = df.pivot(index='y', columns='x', values='v')

    # extracting 1D axes and 2D value arrays
    X = T_grid.columns.values
    Y = T_grid.index.values
    T = T_grid.values
    U = U_grid.values
    V = V_grid.values

    # ---------------------------------------------------------
    # 1. plotting the temperature contours (Isotherms)
    # ---------------------------------------------------------
    plt.figure(figsize=(7, 6))

    contour = plt.contourf(X, Y, T, levels=50, cmap='inferno') 
    plt.colorbar(contour, label='Dimensionless Temperature (T)')
    
    # adding contour lines
    plt.contour(X, Y, T, levels=15, colors='black', linewidths=0.5, alpha=0.5)

    plt.title(f'Isotherms (Temperature) | Ra = 10^{int(np.log10(ra_val))}')
    plt.xlabel('X (Width)')
    plt.ylabel('Y (Height)')
    plt.tight_layout()
    
    t_filename = f'Results/Plot_Temperature_Ra_{int(ra_val)}.png'
    plt.savefig(t_filename, dpi=300)
    print(f"Saved: {t_filename}")
    plt.show()

    # ---------------------------------------------------------
    # 2. plotting the mid-plane velocity Profiles
    # ---------------------------------------------------------
    # finding the column/row indices closest to the geometric center (0.5)
    mid_x_idx = np.argmin(np.abs(X - 0.5))
    mid_y_idx = np.argmin(np.abs(Y - 0.5))

    # extracting U-Velocity along the vertical mid-plane (x = 0.5)
    U_profile = U[:, mid_x_idx]
    # extracting V-Velocity along the horizontal mid-plane (y = 0.5)
    V_profile = V[mid_y_idx, :]

    plt.figure(figsize=(12, 5))

    # creating subplot 1: Vertical Mid-plane (U vs Y)
    plt.subplot(1, 2, 1)
    plt.plot(U_profile, Y, color='blue', linewidth=2, label='U-velocity at x=0.5', marker='o', markersize=4)
    plt.title(f'Horizontal Velocity Profile (Ra = 10^{int(np.log10(ra_val))})')
    plt.xlabel('Dimensionless U-Velocity')
    plt.ylabel('Y (Height)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.axvline(0, color='black', linewidth=1) # Zero-line reference
    plt.legend()

    # creating subplot 2: Horizontal Mid-plane (V vs X)
    plt.subplot(1, 2, 2)
    plt.plot(X, V_profile, color='red', linewidth=2, marker='o', markersize=4, label='V-velocity at y=0.5')
    plt.title(f'Vertical Velocity Profile (Ra = 10^{int(np.log10(ra_val))})')
    plt.xlabel('X (Width)')
    plt.ylabel('Dimensionless V-Velocity')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.axhline(0, color='black', linewidth=1) # Zero-line reference
    plt.legend()


    plt.tight_layout()
    v_filename = f'Results/Plot_Velocities_Ra_{int(ra_val)}.png'
    plt.savefig(v_filename, dpi=300)
    print(f"Saved: {v_filename}")
    plt.show()

# executing the plotting script
if __name__ == "__main__":

    # specifying Rayleigh number values for plotting results
    target_Ra = 1000 # switch to 10000 or 100000
    
    plot_results(target_Ra)