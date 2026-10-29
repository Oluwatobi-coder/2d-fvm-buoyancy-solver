
# specifying fluid & flow parameters
Ra = 1e3      # Rayleigh number, change to 1e4 and 1e5 for different test cases
Pr = 0.71      # Prandtl number (air)

# specifying spatial discretization & domain geometry
nx = 64        # no. of control volumes in x-direction
ny = 64        # no. of control volumes in y-direction
Lx = 1.0       # domain width [m]
Ly = 1.0       # domain height [m]
dx = Lx / nx   # grid spacing in x
dy = Ly / ny   # grid spacing in y

# specifying dimensionless physical properties
nu = Pr        # kinematic viscosity (scaled for dimensionless formulation)
alpha = 1.0    # thermal diffusivity
T_ref = 0.5    # reference temperature

# setting temporal Stability & time-step control
# specifying explicit 2D diffusion von Neumann stability constraint: dt <= 0.5 / [ D * (1/dx^2 + 1/dy^2) ]
safety_factor = 0.8                                           
diffusivity_limit = max(alpha, nu)                            # maximum diffusion coefficient
max_stable_dt = 0.5 / (diffusivity_limit * (1/dx**2 + 1/dy**2)) # theoretical stability threshold
dt = safety_factor * max_stable_dt                            # simulation time step size

# setting max no. of iterations & convergence criteria
max_iters = 100000  # Maximum allowable outer loop iterations
tol = 1e-5          # Convergence residual tolerance