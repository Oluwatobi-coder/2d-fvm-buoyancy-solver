# importing libraries
import numpy as np
from config import *

# funtion to initialize the staggered grid architecture
def initialize_fields():
    # u on vertical faces, v on horizontal, p and T at cell centers
    u = np.zeros((nx + 1, ny + 2)) 
    v = np.zeros((nx + 2, ny + 1))
    p = np.zeros((nx + 2, ny + 2))
    T = np.zeros((nx + 2, ny + 2))
    
    for i in range(nx + 2):
        for j in range(ny + 2):
            T[i, j] = T_ref
            
    return u, v, p, T

# function for specifying boundary conditions
def apply_bcs(T, u, v):
    # Temperature BCs: hot left wall, cold right wall, adiabatic top/bottom
    for j in range(ny + 2):
        T[0, j] = 2.0 * 1.0 - T[1, j]        # Hot left wall
        T[nx + 1, j] = 2.0 * 0.0 - T[nx, j]  # Cold right wall
        
    for i in range(nx + 2):
        T[i, 0] = T[i, 1]                    # Bottom adiabatic
        T[i, ny + 1] = T[i, ny]              # Top adiabatic
        
    # Velocity BCs: no-slip walls
    for j in range(ny + 2):
        u[0, j] = 0.0
        u[nx, j] = 0.0
        
    for i in range(nx + 2):
        v[i, 0] = 0.0
        v[i, ny] = 0.0

# function for solving the pressure-correction Poisson equation
def solve_poisson(p, div):
    p_corr = np.zeros((nx + 2, ny + 2))
    p_old = np.zeros((nx + 2, ny + 2))
    
    # setting the iterative solver for pressure-correction equation
    for _ in range(50): 
        for i in range(nx + 2):
            for j in range(ny + 2):
                p_old[i, j] = p_corr[i, j]
                
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                p_corr[i, j] = 0.5 * (dy**2 * (p_old[i+1, j] + p_old[i-1, j]) +
                                      dx**2 * (p_old[i, j+1] + p_old[i, j-1]) -
                                      dx**2 * dy**2 * div[i, j]) / (dx**2 + dy**2)
                                      
        # applying Neumann BCs for pressure correction
        for j in range(1, ny + 1):
            p_corr[0, j] = p_corr[1, j]
            p_corr[nx + 1, j] = p_corr[nx, j]
        for i in range(1, nx + 1):
            p_corr[i, 0] = p_corr[i, 1]
            p_corr[i, ny + 1] = p_corr[i, ny]
            
    return p_corr

# function to execute solver loop
def step(u, v, p, T):
    apply_bcs(T, u, v)

    # specifying the intermediate velocities fields
    u_star = np.zeros((nx + 1, ny + 2))
    v_star = np.zeros((nx + 2, ny + 1))
    
    for i in range(nx + 1):
        for j in range(ny + 2):
            u_star[i, j] = u[i, j]
    for i in range(nx + 2):
        for j in range(ny + 1):
            v_star[i, j] = v[i, j]

    # setting up the 2D SIMPLE Loop (momentum prediction with Central Differencing scheme)
    
    # solve X-Momentum for intermediate velocity u_star
    for i in range(1, nx):
        for j in range(1, ny + 1):
            # interpolating the diffusion terms 
            u_xx = (u[i+1, j] - 2*u[i, j] + u[i-1, j]) / dx**2
            u_yy = (u[i, j+1] - 2*u[i, j] + u[i, j-1]) / dy**2
            
            # interpolating neighbor velocities at the cell faces
            u_east = 0.5 * (u[i, j] + u[i+1, j])
            u_west = 0.5 * (u[i-1, j] + u[i, j])
            v_north = 0.5 * (v[i, j] + v[i+1, j])
            v_south = 0.5 * (v[i, j-1] + v[i+1, j-1])
            
            # interpolating velocities to cell faces via Central Differencing
            u_e_cd = 0.5 * (u[i, j] + u[i+1, j])
            u_w_cd = 0.5 * (u[i-1, j] + u[i, j])
            u_n_cd = 0.5 * (u[i, j] + u[i, j+1])
            u_s_cd = 0.5 * (u[i, j-1] + u[i, j])
            
            # calculating the convective flux terms: d(u^2)/dx + d(uv)/dy
            conv_x = (u_east * u_e_cd - u_west * u_w_cd) / dx
            conv_y = (v_north * u_n_cd - v_south * u_s_cd) / dy
            
            # advancing u_star 
            u_star[i, j] = u[i, j] + dt * (nu * (u_xx + u_yy) - (conv_x + conv_y))

    # solving Y-Momentum for intermediate velocity v_star
    for i in range(1, nx + 1):
        for j in range(1, ny):
            # interpolating the diffusion terms 
            v_xx = (v[i+1, j] - 2*v[i, j] + v[i-1, j]) / dx**2
            v_yy = (v[i, j+1] - 2*v[i, j] + v[i, j-1]) / dy**2
            
            # interpolating neighbor velocities at the cell faces
            u_east = 0.5 * (u[i, j] + u[i, j+1])
            u_west = 0.5 * (u[i-1, j] + u[i-1, j+1])
            v_north = 0.5 * (v[i, j] + v[i, j+1])
            v_south = 0.5 * (v[i, j-1] + v[i, j])
            
            # interpolating velocities to cell faces via Central Differencing
            v_e_cd = 0.5 * (v[i, j] + v[i+1, j])
            v_w_cd = 0.5 * (v[i-1, j] + v[i, j])
            v_n_cd = 0.5 * (v[i, j] + v[i, j+1])
            v_s_cd = 0.5 * (v[i, j-1] + v[i, j])
            
            # calculating the convective flux terms: d(uv)/dx + d(v^2)/dy[cite: 1]
            conv_x = (u_east * v_e_cd - u_west * v_w_cd) / dx
            conv_y = (v_north * v_n_cd - v_south * v_s_cd) / dy
            
            # calculating the explicit buoyancy force using Boussinesq approximation
            T_face = 0.5 * (T[i, j] + T[i, j+1])
            buoyancy_term = Ra * Pr * (T_face - T_ref)
            
            # advancing v_star
            v_star[i, j] = v[i, j] + dt * (nu * (v_xx + v_yy) - (conv_x + conv_y) + buoyancy_term)
            
    # calculating divergence for Continuity equation (div u = 0)
    div = np.zeros((nx + 2, ny + 2))
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            div[i, j] = (u_star[i, j] - u_star[i-1, j]) / dx + \
                        (v_star[i, j] - v_star[i, j-1]) / dy
                        
    # calculating corrected pressure
    p_corr = solve_poisson(p, div / dt)
    
    # updating the velocities and pressure
    for i in range(1, nx):
        for j in range(1, ny + 1):
            u[i, j] = u_star[i, j] - dt * (p_corr[i+1, j] - p_corr[i, j]) / dx
            
    for i in range(1, nx + 1):
        for j in range(1, ny):
            v[i, j] = v_star[i, j] - dt * (p_corr[i, j+1] - p_corr[i, j]) / dy
            
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            p[i, j] += p_corr[i, j]
            
    # interpolating for Temperature
    T_old = np.zeros((nx + 2, ny + 2))
    for i in range(nx + 2):
        for j in range(ny + 2):
            T_old[i, j] = T[i, j]
            
    max_dt_error = 0.0
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            # calculating thermal energy transport terms (Central Differencing)
            T_xx = (T_old[i+1, j] - 2*T_old[i, j] + T_old[i-1, j]) / dx**2
            T_yy = (T_old[i, j+1] - 2*T_old[i, j] + T_old[i, j-1]) / dy**2
            
            u_c = 0.5 * (u[i-1, j] + u[i, j])
            v_c = 0.5 * (v[i, j-1] + v[i, j])
            
            T_x = (T_old[i+1, j] - T_old[i-1, j]) / (2*dx)
            T_y = (T_old[i, j+1] - T_old[i, j-1]) / (2*dy)
            
            T[i, j] = T_old[i, j] + dt * (alpha * (T_xx + T_yy) - (u_c * T_x + v_c * T_y))
            
            diff = abs(T[i, j] - T_old[i, j])
            if diff > max_dt_error:
                max_dt_error = diff
                
    return max_dt_error