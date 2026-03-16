# This script is used to execute the inflated Dynamic Laplacian (iDL) method on the 
# twice coherent, twice mixing (cmcm) Childress-Soward system, as defined in §7.3 of 
# FK23

# Load in the general functions necessary to execute the method
# (e.g. generating the trajectory data, constructing the diffusion-map matrices etc.)
println("Loading Julia packages and iDL functions...")
include("./iDL_functions_cmcm.jl")

# Load in a function file containing code to plot the eigenvectors produced for the CMCM Childress-Soward system
include("./cmcm_plot.jl")

# Define the CMCM Childress-Soward velocity system from §7.3 of FK23
function cs_velocity_cmcm(du, u, p, t)

    # Define a smooth, continuous indicator function equalling 1 during a mixing regime and 0 during a coherent regime
    tanh_fun(t) = (((1 / 2)*(tanh(1000 * (t - 0.6)) - tanh(1000 * (t - 1.4)))) + ((1 / 2)*(tanh(1000 * (t - 3.21)) - tanh(1000 * (t - 4.01)))))
    # We use this continuous definition to avoid incorrect trajectory data being produced
    # by the ODE solver, which can happen when a function that is not smooth and 
    # differentiable across the entire time interval is passed through to it.

    # Define functions for r (the modulation parameter) and A (the velocity amplitude)
    r(t) = sign(cos(5*π*t))*tanh_fun(t)
    A(t) = 60 - 20*tanh_fun(t)

    # Return the CS velocity
    du[1] = A(t) * (sin(u[1])*cos(u[2]) - r(t)*cos(u[1])*sin(u[2]))
    du[2] = A(t) * (-cos(u[1])*sin(u[2]) + r(t)*sin(u[1])*cos(u[2]))
    
end
# Add the path to the folder to which we wish to save the results data and Figures
pathname = "./CMCM/"

# Define a grid of initial conditions for each trajectory in space
# along with a range of time steps
println("Defining initial conditions and time range...")
# init_conds is a vector of length N, where N = nx*ny is the number of trajectories we compute
# for the CMCM CS flow system. Each element of init_conds is a two-element vector
# [x, y], containing the x and y coordinates of the initial conditions for each
# trajectory.

nx, ny = 30, 30 # Numbers of points in x and y, for a total of N = nx*ny
x_min, x_max = 0, 2π # The minimum and maximum values in x
y_min, y_max = 0, 2π # The minimum and maximum values in y

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
t₀, Δt, τ = 0, 2/75, 4 # Define the initial time, time step and final time
time_steps = t₀:Δt:τ

# Solve the CMCM Childress-Soward system to obtain the trajectory data
println("Generating trajectory data...")
@time traj_solutions = solve_ode(cs_velocity_cmcm, init_conds, time_steps)

# Convert the trajectory data into a more user-friendly format
# i.e. a T-length vector of N-length vectors of length 2, each one containing
# the [x,y] coordinates of trajectory n = 1,...,N at time t = 1,...,T

N = length(traj_solutions)
T = length(time_steps)

# This Childress-Soward flow is defined on a 2π-torus (periodic in both spatial 
# dimensions), so apply mod(⋅,Lx=2π) to the trajectory data at each time step
traj_data_full = [[mod.(traj_solutions[n](time_steps[t]), Lx) for n ∈ 1:N] for t ∈ 1:T]

# Define values for the parameters ϵ and a 
# (the diffusion bandwidth parameter and temporal diffusion strength parameter respectively)
# Define a as listed in FK23, and use the findϵ function on the initial trajectory data to determine a value for ϵ
ϵ = findϵ(traj_data_full[1])
a = 4/π

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
num_of_Λ = 20
tol = 1e-08 # Define a tolerance for Arnoldi method convergence

@time Λ, V, Υ = eigensolve_iDL(Pvec, L_exp, ϵ, num_of_Λ, tol)

# Classify the eigenvalues Υ as spatial or temporal and plot the spectrum
# We shouldn't encounter any complex eigenvalues/eigenvectors in this example, but just in case pass through the real parts of Υ and V only
println("Plotting the spectrum...")
figname = pathname * "iDL_Spectrum_CMCM.png"
spat_inds, temp_inds = classify_eigs_and_plot_spectrum(real.(Υ), real.(V), N, T, figname)

# Plot some leading spatial eigenvectors of the iDL
# Pass through the real part of V to the function in case some eigenvectors are complex (though they shouldn't be in this example)
println("Plotting some eigenvectors...")
inds_to_plot = spat_inds[2:5]
figname = pathname * "iDL_CMCM_First4Evecs"
plot_eigvecs(real.(V), inds_to_plot, traj_data_full, time_steps, figname)

inds_to_plot = spat_inds[6:9]
figname = pathname * "iDL_CMCM_Next4Evecs"
plot_eigvecs(real.(V), inds_to_plot, traj_data_full, time_steps, figname)

# Plot the maxima of SEBA vectors from some of the leading spatial eigenvectors.
# This produces a Figure similar to what we see in Figure 12 of FK23.
# Again, pass through the real part of V to the function in case some eigenvectors are complex (though they shouldn't be in this example)
inds_to_use = spat_inds[2:9]
augment_vecs = false # Decide if you want to calculate SEBA vectors using the eigenvector augmentation technique (true) or using the iDL eigenvectors as they are (false)
figname = pathname * "iDL_CMCM_SEBAMax_8Vecs_NoAug"
Σ_noAug = plot_max_SEBA(real.(V), inds_to_use, traj_data_full, time_steps, figname, augment_vecs)

# Try making another SEBA Max plot with eigenvector augmentation this time
inds_to_use = spat_inds[2:7]
augment_vecs = true
figname = pathname * "iDL_CMCM_SEBAMax_12Vecs_Aug"
Σ_Aug = plot_max_SEBA(real.(V), inds_to_use, traj_data_full, time_steps, figname, augment_vecs)

# Save the trajectory data, the eigenvalue/eigenvector data, the time steps taken and key parameters relevant to the iDL method to a JLD2 file

println("Saving the results...")
filenamesave = pathname * "iDL_Results_CMCM.jld2"
jldsave(filenamesave; Λ, V, Υ, Σ_noAug, Σ_Aug, traj_data_full, time_steps, ϵ, a, spat_inds, temp_inds)

println("The inflated dynamic Laplacian calculations are complete!")