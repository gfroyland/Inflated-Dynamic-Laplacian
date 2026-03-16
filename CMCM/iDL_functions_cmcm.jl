# Load in the required packages
using GLMakie, LinearAlgebra, Distances, Statistics, SparseArrays, OrdinaryDiffEq, LinearMaps, ArnoldiMethod, JLD2, ConcaveHull, NearestNeighbors

# Solve the system of interest from our initial points to produce trajectory data
# velocity is a function representing the velocity of the flow system, defined as velocity(du, u, p, t) to be compatible with Julia's ODE solvers.
# init_conds is a vector of vectors, with each element being a vector of length 2 [x, y] containing the coordinates of each initial condition.
# time_steps is a vector containing the equispaced time points between 0 and τ (inclusive), on each of which a diffusion matrix is computed from the trajectory data. 
function solve_ode(velocity, init_conds, time_steps)

    # Define a singular "sample" ODEProblem for the velocity system
    t_span = (time_steps[1], time_steps[end])
    init_cond_sample = init_conds[1]
    ode_sample = ODEProblem(velocity, init_cond_sample, t_span)

    # Define prob_func, used to create an ensemble of ODEProblems (one for each
    # initial point). prob_func tells the EnsembleProblem operator below what to
    # change with each remade ODEProblem.
    function prob_func(prob, i, repeat)
        remake(prob, u0=init_conds[i])
    end

    # Build the full ensemble of ODEProblems (one for each initial point), and 
    # solve the ensemble of problems
    ode_ensemble = EnsembleProblem(ode_sample, prob_func=prob_func)
    traj_solutions = solve(ode_ensemble, Tsit5(), EnsembleThreads(), trajectories=length(init_conds))

    # Return the ensemble of ODE solutions
    return traj_solutions

end

# Find the indices of the boundary points at each time step when we use Dirichlet BCs
# traj_data_full is a vector of length T (the number of time steps); each element of this vector is a vector of length N (the number of trajectories), and each element of this vector is a vector of length 2 [x, y] representing the coordinates of trajectory n (n = 1,...,N) at time time_steps(t) (t = 1,...,T).
# N is the number of trajectories
# T is the number of time steps
function find_boundary_points(traj_data_full, N, T)

    # Define the number of neighbouring points used to determine whether or not a trajectory point is a boundary point at time t
    # Taking roughly 1% of the total number of trajectories should work
    num_of_neighbours = trunc(Int64, round(N/100))
    
    # Initialise Bmat, a vector of T sparse diagonal matrices of size N × N
    # The n-th (n = 1,...,N) diagonal element of Bmat[t] (Bmat[t][n,n]) will take a value of 0 if the n-th trajectory is a boundary point at time t, or 1 if it is not. 
    Bmat = Vector{SparseMatrixCSC{Float64,Int64}}(undef, T)

    # Determine the boundary points on each time step
    Threads.@threads for t = 1:T
        # First, compute the concave hull for the set of trajectory points at time t
        hull = concave_hull(traj_data_full[t],num_of_neighbours)

        # Then, use indexin to find the indices (ranging from 1 to N) of the trajectories whose points correspond to vertices of the concave hull boundary
        boundary_indices = indexin(hull.vertices, traj_data_full[t])

        # Initialise an N-length vector of ones, and set the values of the vector at the indices corresponding to trajectories on the boundary to 0
        Bvec = ones(N)
        Bvec[boundary_indices] .= 0
        
        # Turn this vector into a sparse diagonal matrix, and remove the zero elements from it
        Bmat_now = spdiagm(Bvec)
        dropzeros!(Bmat_now)
        Bmat[t] = Bmat_now
    end
    
    # Return the boundary indicator matrices
    return Bmat
    
end

# Calculate mean nearest neighbour distances over all trajectory points in space at a particular moment in time
# Code prepared by Gary Froyland
function nndist(X)

    # Input: array X is an n vector of d-vectors;  n data points and d dimensions
    # Output: average distance to the (2d)^th nearest neighbour (averaged over all points)

    d = size(X[1])[1]

    # calculate full pairwise distance array (nxn)
    D = pairwise(Euclidean(), X)

    # find distance from each point to its (2d)^th nearest neighbour
    nn_dist = [sort(D[:, i])[2d+1] for i ∈ eachindex(X)]

    # take the mean of these distances
    average_nn_dist = mean(nn_dist)

    return average_nn_dist

end

# Calculate a suitable value of the spatial diffusion bandwidth parameter ϵ
# To be used later, for now use the values of ϵ quoted in AFK24.
# Code prepared by Gary Froyland
function findϵ(X)

    #X is an N-vector of d-vectors
    N = length(X)

    #cheap way to compute distances, you can use a faster way
    D = pairwise(Euclidean(), X)

    baseϵ = nndist(X)
    #adjust the numbers 10 below to something larger if required (if you don't see flat sections at the start and end)
    loggridϵ = -10*log(2)+log(baseϵ):log(2):10*log(2)+log(baseϵ)

    logsumvec = []
    for logϵ = loggridϵ
        #S is your matrix P using exp(logϵ) as the ϵ
        S = exp.(-D .^ 2 / exp(logϵ))
        push!(logsumvec, log(sum(S) / N^2))
    end

    #plot(loggridϵ, logsumvec, label=false)
    #display(scatter!(loggridϵ, logsumvec, label=false, xlabel="log(ϵ)", ylabel="log(ΣP / N²)"))

    centdiff = (logsumvec[3:end] - logsumvec[1:end-2]) / (2 * log(2))

    #plot(loggridϵ[2:end-1], centdiff, label=false)
    #display(scatter!(loggridϵ[2:end-1], centdiff, xlabel="log(ϵ)", label=false, ylabel="central difference of sigmoid function"))

    println("Maximum slope of sigmoid derivative occurs at approximately ϵ = ", exp(loggridϵ[argmax(centdiff)+1]))

    return exp(loggridϵ[argmax(centdiff)+1])

end

# Build the diffusion-map matrices on each time step
# Code prepared by Gary Froyland
function diffusion_matrix(X, ϵ)

    #X is a vector of vectors, e.g. 1000 random points in ℝ²:  X=[rand(2).*[2, 1] for i=1:1000] 
    #ϵ is the bandwidth parameter
    #Example call with X as above:  P, λ, λscaled, v, vscaled = diffusion_matrix(X,0.5)

    n = length(X)

    k = exp.(-pairwise(SqEuclidean(), X) / 4ϵ)
    k[k.<0.007] .= 0

    # normalise for density of points
    Σk = sum(k, dims=2)
    kbar = [k[i, j] / (Σk[i] * Σk[j]) for i = 1:n, j = 1:n]

    # normalise to make a Markov matrix
    P = sparse(stack(normalize!.(eachcol(kbar), 1)))

    return P

end

"""
    sparse_gaussian(X; eps, tau)

X: dim×numpts matrix with points as columns (dim=2 or 3 typical).
Returns sparse S with S[i,j] = exp(-||xi-xj||^2 / 4eps) for entries > tau.
"""
function sparse_gaussian(X, ϵ, τ=0.007)
    ϵ > 0 || throw(ArgumentError("ϵ must be > 0"))
    (0 < τ < 1) || throw(ArgumentError("tau must satisfy 0 < τ < 1"))

    dim, numpts = size(X)
    # compute distance threshold r corresponding to kernel value threshold τ
    r = sqrt(4ϵ * log(inv(τ)))      

    #construct Tree
    tree = KDTree(X)
    pertree = PeriodicTree(tree, [0.0, 0.0], [2π, 2π])

    I = Int[]
    J = Int[]
    V = eltype(X)[]
    sizehint!(I, 20numpts)
    sizehint!(J, 20numpts)
    sizehint!(V, 20numpts)  # cheap heuristic

    for j in 1:numpts
        for i in inrange(pertree, X[:,j], r)
            # build upper triangle only
            if i < j
                continue
            end
            d² = peuclidean(X[:, i], X[:, j], [2π, 2π])^2
            # distance below for cylinder [0,20] x [-3,3] with x-coord periodic
            # d² = peuclidean(X[:, i], X[:, j], [20, Inf])^2

            #compute kernel entry
            s = exp(-d² / 4ϵ)
            # if the value of s is too small, don't include (effectively truncate to zero)
            if s ≤ τ
                continue
            end
            push!(I, i)
            push!(J, j)
            push!(V, s)
            if i ≠ j
                (push!(I, j); push!(J, i); push!(V, s)) # mirror to make symmetric
            end
        end
    end

    K = sparse(I, J, V, numpts, numpts)
    
    # normalise for density of points
    ΣK = sum(K, dims=2)
    Vbar = zeros(length(V))

    for k ∈ eachindex(V)
        Vbar[k] = V[k] / (ΣK[I[k]] * ΣK[J[k]])
    end

    Kbar = sparse(I, J, Vbar, numpts, numpts)

    # normalise to make a Markov matrix
    Kbar_colsums = sum(Kbar, dims=1)
    for n ∈ eachindex(Vbar)
        if Kbar_colsums[J[n]] != 0
            Vbar[n] /= Kbar_colsums[J[n]]
        end
    end

    P = sparse(I, J, Vbar, numpts, numpts)
    return P

end

# This function is used to build diffusion-map matrices on each time step using the trajectory data
# traj_data_full is a vector of length T (the number of time steps); each element of this vector is a vector of length N (the number of trajectories), and each element of this vector is a vector of length 2 [x, y] representing the coordinates of trajectory n (n = 1,...,N) at time time_steps(t) (t = 1,...,T).
# ϵ is the diffusion bandwidth parameter, taken either from the findϵ function defined earlier or as listed in AFK24.
# T is the number of time steps taken from 0 to τ (inclusive), and hence the number of diffusion-map matrices we need to compute.
function make_operators(traj_data_full, ϵ, T)

    # Compute the diffusion-map matrices for each of the T time steps, which will be stored
    # in a T length vector of matrices Pvec.
    Pvec = Vector{SparseMatrixCSC{Float64,Int64}}(undef, T)
    Threads.@threads for t = 1:T
        Pvec[t] = diffusion_matrix(traj_data_full[t], ϵ)
    end

    # Pvec is then returned to the script
    return Pvec

end

function make_operators_new(traj_data_full, ϵ, T)

    n_dims = length(traj_data_full[1][1])
    n_pts = length(traj_data_full[1])

    # Compute the diffusion-map matrices for each of the T time steps, which will be stored
    # in a T length vector of matrices Pvec.
    Pvec = Vector{SparseMatrixCSC{Float64,Int64}}(undef, T)
    Threads.@threads for t = 1:T
        traj_data_now = zeros(n_dims, n_pts)
        for k = 1:n_pts
            traj_data_now[:, k] = traj_data_full[t][k]
        end
        Pvec[t] = sparse_gaussian(traj_data_now, ϵ)
    end

    # Pvec is then returned to the script
    return Pvec

end

# In this function, we build the exponential time Laplacian matrix used to finish building the inflated dynamic Laplacian.
# T is the number of time steps taken between 0 and τ (inclusive)
# Δt is the step size between each time step, i.e. τ/(T-1)
# ϵ is the diffusion bandwidth parameter used earlier to build the diffusion-map matrices
# a is the temporal diffusion strength parameter
function make_time_Laplacian(T, Δt, ϵ, a)

    # Build the second order, central-difference approximation of the Laplacian for time
    L = Tridiagonal(ones(T - 1), -2 * ones(T), ones(T - 1))
    L[1, 1] = -1
    L[T, T] = -1

    # Build the exponential Laplacian, as defined in AFK24, and return it
    L_exp = exp(ϵ * ((a^2) / 2) * (L / (Δt)^2)) # Based on eq. 5.3 of AFK24
    return L_exp

end

# In this function, we eigensolve the discretised inflated Dynamic Laplacian (iDL) using the Arnoldi method and the matrix-vector multiplication technique utilised in eqs. 5.3 and 5.4 of AFK24.
# Pvec is a vector of length T (the number of time steps), with each element of this vector being the diffusion matrix at time time_steps(t) (t = 1,...,T)
# L_exp is the exponential Laplacian operator acting in time
# ϵ is the diffusion bandwidth parameter
# num_of_Λ is the number of leading eigenvalues/eigenvectors we wish to compute for the inflated Dynamic Laplacian
# tol is the convergence tolerance level for the Arnoldi method
# If we are using Dirichlet BCs, we pass through Bmat to this function, a T-length vector of sparse diagonal matrices with the diagonal elements Bmat[t][n,n] equalling 0 if the n-th trajectory is on the boundary at time t, or 1 otherwise
# If we are using Neumann BCs, Bmat is defined as a nothing object and only the five previously listed variables are passed through to the function
function eigensolve_iDL(Pvec, L_exp, ϵ, num_of_Λ, tol, Bmat=nothing)

    # Recover the number of time steps T and the number of trajectories N
    T = size(L_exp, 1)
    N = size(Pvec[1], 1)

    # Define a function for efficient matrix-vector multiplication for the discrete
    # inflated Dynamic Laplacian as described in AFK24 (eqs. (5.3) and (5.4)).
    # The definition of the multiplication scheme depends on the boundary conditions
    # defined for our spatial domain.
    function A(x)

        # Rearrange the NT-length vector x into a T × N matrix F, as defined in eq. (5.4) of AFK24
        F = permutedims(reshape(x, (N, T)))

        # Perform the Strang splitting matrix multiplication detailed in eqs. (5.3)
        # and (5.4) of AFK24 and return the result
        if Bmat ≡ nothing # If we are using Neumann BCs
            temp = L_exp * F # First multiplication step
            Pₐ = reshape(permutedims(L_exp * stack((Pvec[k] * temp[k, :] for k ∈ 1:T), dims=1)), N * T) # Remaining multiplication steps
        else # If we are using Dirichlet BCs
            # 1. Apply the boundary indicator matrix Bmat[k] for time k to each row of F, then multiply the result by L_exp
            Mat1 = L_exp * stack((Bmat[k] * F[k, :] for k ∈ 1:T), dims=1)
            # 2. Perform the next few multiplication steps (all except the last one)
            Mat2 = L_exp * stack((Bmat[k] * Pvec[k] * Bmat[k] * Mat1[k, :] for k ∈ 1:T), dims=1)
            # 3. Perform one last multiplication with Bmat, reshape the matrix back to a vector and return the result
            Pₐ = reshape(permutedims(stack((Bmat[k] * Mat2[k, :] for k ∈ 1:T), dims=1)), N * T)
        end
        return Pₐ

    end

    # Define a linear map for the above multiplication function
    A_map = LinearMap(A, N * T)

    # Pass the linear map through to partialschur to execute the Arnoldi method
    Decomp, History = partialschur(A_map, which=:LM, nev=num_of_Λ, tol=tol)
    # Extract the eigenvalues and eigenvectors (Λ is of length num_of_Λ, V is a matrix 
    # of size (NT) × num_of_Λ, where each column is a spacetime eigenvector for the 
    # inflated Dynamic Laplacian)
    Λ, V = partialeigen(Decomp)

    # Sort the eigenvalues/eigenvectors by their magnitude so that they are indexed in the correct order 
    # (i.e. decreasing in magnitude with Λ[1] ≈ 1) and return them to the script

    Λ_mag = abs.(Λ)
    Λ_ord = sortperm(Λ_mag, rev=true)

    Λ = Λ[Λ_ord]
    V = V[:, Λ_ord]

    # Alternative code using Arpack eigs (no need to sort Λ and V, though Arnoldi generates the same results and faster so we use that method)
    #Λ, V = eigs(A_map, nev=num_of_Λ, which=:LM, tol=tol)

    # Back transform Λ to obtain the eigenvalues of the inflated dynamic Laplacian 0 = Υ_1 > Υ_2 > Υ_3 > ⋯
    # rather than the eigenvalues of the exponentiated discrete version of this operator
    Υ = (1 / ϵ) * log.(Λ)

    # Return the eigenbasis (use real() later for the case(s) where complex eigenvalues/eigenvectors are computed, this shouldn't happen in the examples detailed in AFK24, but just in case)
    return Λ, V, Υ

end

# Use this function to distinguish spatial inflated Dynamic Laplacian eigenvalues from temporal ones, then plot the spectrum with these distinctions illustrated (using different colours for different eigenvalue types)
# Υ and V are the eigenvalues and eigenvectors (respectively) of the inflated Dynamic Laplacian
# N is the number of trajectories
# T is the number of time steps
# figname is the name of the image file (including path) used to save the spectrum figure to
function classify_eigs_and_plot_spectrum(Υ, V, N, T, figname)

    # Retrieve the number of eigenvalues/eigenvectors computed 
    K = length(Υ)

    # Compute the mean spatial variance for each eigenvector over time
    averagespatialvariance = [mean([var(V[(t-1)*N+1:t*N, k]) for t = 1:T]) for k = 1:K]

    # Use the mean spatial variances to distinguish spatial eigenvalues from temporal ones
    tol = 1e-06
    spat_inds = findall(x -> x > tol, averagespatialvariance)
    temp_inds = findall(x -> x < tol, averagespatialvariance)

    # Make sure trivial Υ_1 is classified as spatial, not temporal
    if (~isempty(temp_inds))
        if (abs(temp_inds[1] - 1) < 1e-10)
            popfirst!(temp_inds)
            prepend!(spat_inds, 1)
        end
    end

    # Plot and save the spectrum

    with_theme(theme_latexfonts()) do
        fig = Figure(size=(800, 400))
        ax = Axis(fig[1, 1], xlabel=L"$k$", ylabel=L"$Λ_k$")

        xs = spat_inds
        ys = Υ[spat_inds]
        points = Point2f.(xs, ys)
        scatter!(ax, points, label=L"$\mathrm{Spatial}$", marker=:circle, color=:blue, markersize=10)

        xs = temp_inds
        ys = Υ[temp_inds]
        points = Point2f.(xs, ys)
        scatter!(ax, points, label=L"$\mathrm{Temporal}$", marker=:circle, color=:red, markersize=10)

        hidespines!(ax)
        fig[1, 2] = Legend(fig, ax, "Eigenvalue Types", framevisible=false)

        display(fig)
        save(figname, fig, px_per_unit=2.0)
    end

    # Return the eigenvalue classification information
    return spat_inds, temp_inds

end