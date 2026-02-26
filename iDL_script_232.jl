# This script is used to execute the inflated Dynamic Laplacian (iDL) method on the 
# 2-3-2 (coherent-mixing-coherent) Double Gyre system, as defined in §6.3 of AFK24

# Load in the general functions necessary to execute the method
# (e.g. generating the trajectory data, constructing the diffusion-map matrices etc.)
println("Loading Julia packages and iDL functions...")
include("./iDL_functions.jl")

# Load in a function file containing code to plot the eigenvectors produced for the 2-3-2 Double Gyre system
include("./dg_232_plot.jl")

# Define the 2-3-2 Double Gyre velocity system from §6.3 of AFK24
function dg_velocity_232(du, u, p, t)

    # Define the time-dependent double gyre vector field F(t,x)
    A, ω = 0.25, π

    # Define ε in a continuous fashion, so that ε(5) = ε(10) = 0.25, while also ensuring
    # that the ODE solver used later will not produce incorrect trajectory data (which can happen
    # when a function that is not smooth and differentiable across the entire time interval
    # is passed through to it)
    ε(t) = 0.25 * ((1 / 2)*(tanh(1000 * (t - 4.9)) - tanh(1000 * (t - 10.1))))

    f(t, u) = ε(t) * sin(ω * t) * u[1]^2 + (1 - 2 * ε(t) * sin(ω * t)) * u[1]
    df(t, u) = (2 * ε(t) * sin(ω * t) * u[1]) + (1 - 2 * ε(t) * sin(ω * t))

    du[1] = -π * A * sin(π * f(t, u)) * cos(π * u[2])
    du[2] = π * A * cos(π * u[1]) * sin(π * u[2]) * df(t, u)
    
end
# Add the path to the folder to which we wish to save the results data and Figures
pathname = "./DG_232/"

# Define a grid of initial conditions for each trajectory in space
# along with a range of time steps
println("Defining initial conditions and time range...")
# init_conds is a vector of length N, where N = nx*ny is the number of trajectories we compute
# for the 2-3-2 DG flow system. Each element of init_conds is a two-element vector
# [x, y], containing the x and y coordinates of the initial conditions for each
# trajectory.

nx, ny = 40, 20 # Numbers of points in x and y, for a total of N = nx*ny
x_min, x_max = 0, 2 # The minimum and maximum values in x
y_min, y_max = 0, 1 # The minimum and maximum values in y

# Define the lengths of our spatial domain M in each dimension, and the mesh sizes in each direction
Lx, Ly = (x_max-x_min), (y_max-y_min)
Δx, Δy = Lx/nx, Ly/ny

# Define the initial condition coordinate ranges in x and y...
x_range = (x_min+(Δx/2)):Δx:(x_max-(Δx/2))
y_range = (y_min+(Δy/2)):Δy:(y_max-(Δy/2))

# ...and use these to build the vector of initial conditions.
init_conds = [[x, y] for x ∈ x_range for y ∈ y_range]

# time_steps is a vector of length T, where each element is the t-th time step between
# t0 and τ (inclusive)
t₀, Δt, τ = 0, 0.2, 15 # Define the initial time, time step and final time
time_steps = t₀:Δt:τ

# Define the boundary conditions in space for this system using this dirichlet Boolean variable 
# Set dirichlet to "true" for Dirichlet BCs, or "false" for Neumann BCs
dirichlet = true

# Solve the 2-3-2 Double Gyre system to obtain the trajectory data
println("Generating trajectory data...")
@time traj_solutions = solve_ode(dg_velocity_232, init_conds, time_steps)

# Convert the trajectory data into a more user-friendly format
# i.e. a T-length vector of N-length vectors of length 2, each one containing
# the [x,y] coordinates of trajectory n = 1,...,N at time t = 1,...,T

N = length(traj_solutions)
T = length(time_steps)

traj_data_full = [[traj_solutions[n](time_steps[t]) for n ∈ 1:N] for t ∈ 1:T]

# Define values for the parameters ϵ and a 
# (the diffusion bandwidth parameter and temporal diffusion strength parameter respectively)
# For now, use the parameters listed in AFK24 for this system
ϵ = 0.032/4
# Even though ϵ is quoted as 0.032 in AFK24 and defined as 0.032 in the MATLAB code,
# the reason we divide this parameter by 4 is that in the MATLAB code, the kernel
# used to construct the diffusion-map matrices is defined as:
# k(x,y) = exp(||x-y||^2) / ϵ)
# rather than
# k(x,y) = exp(||x-y||^2) / 4ϵ)
# so 0.032 is still the denominator in the exponential argument. Later on, however,
# the exponential time Laplacian matrix is defined as
# exp((ϵ/4) * ((a^2) / 2) * (Lᵗ / (Δt)^2))
# and when the eigenvalues Λ of the exponentiated, discrete iDL are back transformed 
# later to obtain the actual iDL eigenvalues 0 = Υ_1 > Υ_2 > Υ_3 > ..., 
# this formula is used:
# Υ = (1/(ϵ/4))*log(Λ)
# Hence, to avoid altering any code in iDL_functions.jl, ϵ = 0.032/4 = 0.008 here,
# which should still give the same results as prepared by the MATLAB code 
# and as shown in AFK24.
#a = 30.0601 # In the MATLAB code, this is the a value used for Neumann boundary conditions...
a = 67.2165 # ...and this is the value used for Dirichlet boundary conditions

# Construct the diffusion-map matrices Pˣ(ϵ) for each time step
# Pvec is a T-length vector of N × N matrices, each matrix representing Pˣ(ϵ) for tₖ, k = 1,...,T
println("Constructing the diffusion-map matrices...")
@time Pvec = make_operators(traj_data_full, ϵ, T)

# Construct the exponential time Laplacian matrix for the inflated Dynamic Laplacian

println("Constructing the time Laplacian matrix...")
L_exp = make_time_Laplacian(T, Δt, ϵ, a)

# Obtain the leading num_of_Λ eigenvalues and eigenvectors of the inflated Dynamic Laplacian
# The eigenvalues of the inflated Dynamic Laplacian Υ = (1/ϵ)*log(Λ)

println("Eigensolving the inflated dynamic Laplacian...")
num_of_Λ = 10
tol = 1e-08 # Define a tolerance for Arnoldi method convergence

# If we are using Dirichlet BCs, identify the boundary points on each time step first and use these to build the boundary indicator matrices Bmat[t] on each time step to be used in the multiplication scheme for the iDL
if dirichlet == true
    Bmat = find_boundary_points(traj_data_full, N, T)
    @time Λ, V, Υ = eigensolve_iDL(Pvec, L_exp, ϵ, num_of_Λ, tol, Bmat)
elseif dirichlet == false
    @time Λ, V, Υ = eigensolve_iDL(Pvec, L_exp, ϵ, num_of_Λ, tol)
end

# Classify the eigenvalues Υ as spatial or temporal and plot the spectrum
# We shouldn't encounter any complex eigenvalues/eigenvectors in this example, but just in case pass through the real parts of Υ and V only
println("Plotting the spectrum...")
figname = pathname * "iDL_Spectrum_DG_232.png"
spat_inds, temp_inds = classify_eigs_and_plot_spectrum(real.(Υ), real.(V), N, T, figname)

# Plot some leading spatial eigenvectors of the iDL
# Pass through the real part of V to the function in case some eigenvectors are complex (though they shouldn't be in this example)
inds_to_plot = spat_inds[2:4]
println("Plotting some eigenvectors...")
figname = pathname * "iDL_DG_232_First3Evecs"
plot_eigvecs(real.(V), inds_to_plot, traj_data_full, time_steps, figname)

# Plot SEBA vectors from some of the leading spatial eigenvectors and their augmented counterparts.
# This produces Figures similar to what we see in Figure 23 (left and right) of AFK24.
# Again, pass through the real part of V to the function in case some eigenvectors are complex (though they shouldn't be in this example)
inds_to_use = spat_inds[2:4]
figname = pathname * "iDL_DG_232_SEBA_6Vecs"
Σ = plot_SEBA(real.(V), inds_to_use, traj_data_full, time_steps, figname)

# Save the trajectory data, the eigenvalue/eigenvector data, the time steps taken and key parameters relevant to the iDL method to a JLD2 file

println("Saving the results...")
filenamesave = pathname * "iDL_Results_DG_232.jld2"
jldsave(filenamesave; Λ, V, Υ, Σ, traj_data_full, time_steps, ϵ, a, spat_inds, temp_inds)

println("The inflated dynamic Laplacian calculations are complete!")