# This script is used to execute the inflated Dynamic Laplacian (iDL) method on the 
# switching Double Gyre system, as defined in §6.1 of AFK24

# Load in the general functions necessary to execute the method
# (e.g. generating the trajectory data, constructing the diffusion-map matrices etc.)
println("Loading Julia packages and iDL functions...")
include("./iDL_functions_dg.jl")

# Load in a function file containing code to plot the eigenvectors produced for the switching Double Gyre system
include("./dg_switching_plot.jl")

# Define the switching Double Gyre velocity system from §6.1 of AFK24
function dg_velocity_switching(du, u, p, t)

    # Define the time-dependent switching double gyre vector field F(t,x) see [Atnip/Froyland/Koltai, 2024]
    r(t) = (1 / 2) * (1 + tanh(40 * (t - (1 / 2))))
    α(t) = (1 - 2 * r(t)) / (3 * (r(t) - 2) * (r(t) + 1))
    β(t) = (2 - 9 * α(t)) / 3

    du[1] = 20 * (-π / 2) * sin(π * (α(t) * u[1]^2 + β(t) * u[1])) * cos(π * u[2] / 2)
    du[2] = 20 * π * (2 * α(t) * u[1] + β(t)) * cos(π * (α(t) * u[1]^2 + β(t) * u[1])) * sin(π * u[2] / 2)

end

# Add the path to the folder to which we wish to save the results data and Figures
pathname = "./"

# Define a grid of initial conditions for each trajectory in space
# along with a range of time steps
println("Defining initial conditions and time range...")
# init_conds is a vector of length N, where N = nx*ny is the number of trajectories we compute
# for the switching DG flow system. Each element of init_conds is a two-element vector
# [x, y], containing the x and y coordinates of the initial conditions for each
# trajectory.

nx, ny = 45, 30 # Numbers of points in x and y, for a total of N = nx*ny
x_min, x_max = 0, 3 # The minimum and maximum values in x
y_min, y_max = 0, 2 # The minimum and maximum values in y

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
t₀, Δt, τ = 0, 0.01, 1 # Define the initial time, time step and final time
time_steps = t₀:Δt:τ

# Define the boundary conditions in space for this system using this dirichlet Boolean variable 
# Set dirichlet to "true" for Dirichlet BCs, or "false" for Neumann BCs
dirichlet = false

# Solve the switching Double Gyre system to obtain the trajectory data
println("Generating trajectory data...")
@time traj_solutions = solve_ode(dg_velocity_switching, init_conds, time_steps)

# Convert the trajectory data into a more user-friendly format
# i.e. a T-length vector of N-length vectors of length 2, each one containing
# the [x,y] coordinates of trajectory n = 1,...,N at time t = 1,...,T

N = length(traj_solutions)
T = length(time_steps)

traj_data_full = [[traj_solutions[n](time_steps[t]) for n ∈ 1:N] for t ∈ 1:T]

# Define values for the parameters ϵ and a 
# (the diffusion bandwidth parameter and temporal diffusion strength parameter respectively)
# For now, use the parameters listed in AFK24 for this system
ϵ = 0.0115
a = 1.0

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
figname = pathname * "iDL_Spectrum_DG_switching.png"
spat_inds, temp_inds = classify_eigs_and_plot_spectrum(real.(Υ), real.(V), N, T, figname)

# Plot some leading spatial eigenvectors of the iDL
# Pass through the real part of V to the function in case some eigenvectors are complex (though they shouldn't be in this example)
inds_to_plot = spat_inds[2:4]
println("Plotting some eigenvectors...")
figname = pathname * "iDL_DG_switching_First3Evecs"
plot_eigvecs(real.(V), inds_to_plot, traj_data_full, time_steps, figname)

# Plot some leading spatial eigenvectors, this time only for the trajectories corresponding to the top 5% of eigenvector variance values over time
# This produces Figures similar to what we see in Figure 11 (left and right) of AFK24.
# Again, pass through the real part of V to the function in case some eigenvectors are complex (though they shouldn't be in this example)
inds_to_plot_topvar = spat_inds[2:4]
figname = pathname * "iDL_DG_switching_First3Evecs_TopVar"
plot_eigvecs_topvar(real.(V), inds_to_plot_topvar, traj_data_full, time_steps, figname)

# Save the trajectory data, the eigenvalue/eigenvector data, the time steps taken and key parameters relevant to the iDL method to a JLD2 file

println("Saving the results...")
filenamesave = pathname * "iDL_Results_DG_switching.jld2"
jldsave(filenamesave; Λ, V, Υ, traj_data_full, time_steps, ϵ, a, spat_inds, temp_inds)

println("The inflated dynamic Laplacian calculations are complete!")