# This file contains functions used to plot eigenvectors of the inflated Dynamic Laplacian for the CMCM Childress-Soward system defined in FK23.

# We use this function to plot some of the leading non-trivial spatial eigenvectors of the inflated Dynamic Laplacian
# V is an NT × num_of_Λ matrix containing the spacetime eigenvector data for the leading num_of_Λ eigenvectors, given N trajectories and T time steps
# inds_to_plot is a vector listing the indices of the eigenvectors to plot
# traj_data_full is a vector of length T (the number of time steps); each element of this vector is a vector of length N (the number of trajectories), and each element of this vector is a vector of length 2 [x, y] representing the coordinates of trajectory n (n = 1,...,N) at time time_steps(t) (t = 1,...,T).
# time_steps is a vector of length T containing the equispaced time points between 0 and τ (inclusive). 
# figname is the filename (including path) for the Figures produced by this function. The Figures will be saved as pngs (with the file extension to be added in this function), with one version of the eigenvector Figures prepared for 𝕄₀ and another for 𝕄₁ (see below for definitions of 𝕄₀ and 𝕄₁)
function plot_eigvecs(V, inds_to_plot, traj_data_full, time_steps, figname)

    # Recover the number of time steps T and the number of trajectories N
    T = length(time_steps)
    N = length(traj_data_full[1])

    # Arrange the trajectory data into vectors of coordinate triples for plotting, (t, x, y) in two forms: 
    # 𝕄₀ (the initial conditions copied T times for each time slice, (t, x(0), y(0))) and   
    # 𝕄₁ (the full spacetime trajectory data for each initial condition, i.e. (t, x(t), y(t)))
    traj_data_M0 = [(time_steps[t], traj_data_full[1][x][1], traj_data_full[1][x][2]) for t ∈ 1:T for x ∈ 1:N]
    traj_data_M1 = [(time_steps[t], traj_data_full[t][x][1], traj_data_full[t][x][2]) for t ∈ 1:T for x ∈ 1:N]

    # Plot the eigenvectors of choice in 𝕄₀
    with_theme(theme_latexfonts()) do
        figwidth = 400*length(inds_to_plot)+50 # Choose an appropriate Figure width for the number of eigenvectors to be plotted
        fig = Figure(size=(figwidth, 300))
        for ind ∈ eachindex(inds_to_plot) # Loop over each of the eigenvectors
            # Retrieve the data for eigenvector inds_to_plot[ind]
            eigvec_data = V[:, inds_to_plot[ind]]
            traj_pts = traj_data_M0

            # Rescale the eigenvector data to have maximum value 1 (or minimum value -1)
            v_max = maximum(abs.(eigvec_data))
            eigvec_data = eigvec_data./v_max
            col_lims = (-1, 1)

            # Apply a threshold to the eigenvector data to isolate key coherent set behaviour identified through more extremal values of the eigenvector
            thresh = 0.5 # The threshold used in AFK24
            plotindices = findall(abs.(eigvec_data) .> thresh)
            eigvec_data = eigvec_data[plotindices]
            traj_pts = traj_pts[plotindices]

            # Set Figure column to the current eigenvector index (1,...,length(inds_to_plot))
            ax = Axis3(fig[1, ind], xlabel=L"$t$", ylabel=L"$x$", zlabel=L"$y$")
            scatterplot = scatter!(ax, traj_pts, color=real.(eigvec_data), colormap=Reverse(:RdBu), colorrange=col_lims, transparency=false, markersize=5) #alpha=0.1
            ax.limits = ((0, 4), (0, 2π), (0, 2π))
            ax.azimuth = 245π/180
            ax.elevation = 20π/180

            # Add a colorbar to the right of the last eigenvector plot
            if (ind == length(inds_to_plot))
                Colorbar(fig[1, ind+1], scatterplot, colorrange=col_lims)
            end
        end

        display(fig)
        save(figname * "_M0.png", fig, px_per_unit=2.0)
    end

    # Plot the eigenvectors of choice in 𝕄₁
    with_theme(theme_latexfonts()) do
        figwidth = 400*length(inds_to_plot)+50 # Choose an appropriate Figure width for the number of eigenvectors to be plotted
        fig = Figure(size=(figwidth, 300))
        for ind ∈ eachindex(inds_to_plot)
            # Retrieve the data for eigenvector inds_to_plot[ind]
            eigvec_data = V[:, inds_to_plot[ind]]
            traj_pts = traj_data_M1

            # Rescale the eigenvector data to have maximum value 1 (or minimum value -1)
            v_max = maximum(abs.(eigvec_data))
            eigvec_data = eigvec_data./v_max
            col_lims = (-1, 1)

            # Apply a threshold to the eigenvector data to isolate key coherent set behaviour identified through more extremal values of the eigenvector
            thresh = 0.5 # The threshold used in AFK24
            plotindices = findall(abs.(eigvec_data) .> thresh)
            eigvec_data = eigvec_data[plotindices]
            traj_pts = traj_pts[plotindices]

            # Set Figure column to the current eigenvector index (1,...,length(inds_to_plot))
            ax = Axis3(fig[1, ind], xlabel=L"$t$", ylabel=L"$x$", zlabel=L"$y$")
            scatterplot = scatter!(ax, traj_pts, color=real.(eigvec_data), colormap=Reverse(:RdBu), colorrange=col_lims, transparency=false, markersize=5) #alpha=0.1
            ax.limits = ((0, 4), (0, 2π), (0, 2π))
            ax.azimuth = 245π/180
            ax.elevation = 20π/180

            # Add a colorbar to the right of the last eigenvector plot
            if (ind == length(inds_to_plot))
                Colorbar(fig[1, ind+1], scatterplot, colorrange=col_lims)
            end
        end

        display(fig)
        save(figname * "_M1.png", fig, px_per_unit=2.0)
    end

end

# Load in the SEBA function for use below
include("../SEBA.jl")

# This function plots the maxima of SEBA vectors for the iDL in similar fashion to Figure 12 of FK23
# The inputs for this function are similar to those passed in to plot_eigvecs above (inds_to_plot is now called inds_to_use, as instead of plotting eigenvectors we use them to prepare SEBA vectors in this function)
# We can prepare SEBA vectors either using spatial iDL eigenvectors on their own (in which case the Boolean variable augment_vecs is set to false), or by using eigenvector augmentation described in AFK24 (in which case augment_vecs is set to true)
function plot_max_SEBA(V, inds_to_use, traj_data_full, time_steps, figname, augment_vecs=false)

    # Recover the number of time steps T and the number of trajectories N
    T = length(time_steps)
    N = length(traj_data_full[1])

    # If we are using the augmented SEBA vectors technique,
    # add the augmented counterparts to the eigenvectors of choice
    if augment_vecs == true
    
        num_SEBA_vecs = 2*length(inds_to_use)
        aug_vecs_for_SEBA = zeros(N*T, num_SEBA_vecs) # There will be 2*length(inds_to_use) vectors prepared to hand in to SEBA in this case
        for η ∈ eachindex(inds_to_use)
            # Every odd numbered vector is an iDL eigenvector
            aug_vecs_for_SEBA[:, (2*η-1)] = V[:, inds_to_use[η]]

            # Every even numbered vector is an augmented companion vector for each eigenvector
            # The vectors are constant on each time slice, with values equal to the L²-norm of the eigenvector of interest on each time slice
            time_slice_norms = [norm(V[(τ-1)*N.+(1:N), inds_to_use[η]]) for τ ∈ 1:T]
            aug_vecs_for_SEBA[:, (2*η)] = repeat(time_slice_norms, inner=N)
        end

        # Compute the SEBA vectors using the eigenvectors and their augmented counterparts
        println("Computing SEBA vectors...")
        @time Σ, ℛ = SEBA(aug_vecs_for_SEBA)
    
    elseif augment_vecs == false
        # Otherwise, compute the SEBA vectors from the length(inds_to_use) eigenvectors passed in
        println("Computing SEBA vectors...")
        @time Σ, ℛ = SEBA(V[:, inds_to_use])
    end

    # Report the SEBA minima back to the user
    println("The respective SEBA vector minima are ", minimum(Σ, dims=1))

    # Take the maxima of the SEBA vectors
    Σ_max = maximum(Σ, dims=2)

    # Plot the maxima of the SEBA vectors (𝕄₁ only this time)
    # Arrange the trajectory data into vectors of coordinate triples for plotting, (t, x, y) in
    # 𝕄₁ (the full spacetime trajectory data for each initial condition, i.e. (t, x(t), y(t)))
    traj_data_M1 = [(time_steps[t], traj_data_full[t][x][1], traj_data_full[t][x][2]) for t ∈ 1:T for x ∈ 1:N]

    # Plot the SEBA vector maxima in 𝕄₁
    with_theme(theme_latexfonts()) do
        figwidth = 450 
        fig = Figure(size=(figwidth, 300))
        
        # Retrieve the SEBA maxima data
        eigvec_data = Σ_max
        traj_pts = traj_data_M1
        col_lims = (0, 1)

        # Apply a threshold to the SEBA data to isolate key coherent set behaviour identified through more extremal values of the vectors
        thresh = 0.5 
        plotindices = findall(eigvec_data .> thresh)
        eigvec_data = eigvec_data[plotindices]
        traj_pts = traj_pts[plotindices]

        ax = Axis3(fig[1, 1], xlabel=L"$t$", ylabel=L"$x$", zlabel=L"$y$")
        scatterplot = scatter!(ax, traj_pts, color=eigvec_data, colormap=:Reds, colorrange=col_lims, transparency=false, markersize=5) #alpha=0.1
        ax.limits = ((0, 4), (0, 2π), (0, 2π))

        ax.azimuth = 245π/180
        ax.elevation = 20π/180

        # Add a colorbar to the right of the plot
        Colorbar(fig[1, 2], scatterplot, colorrange=col_lims)

        display(fig)
        save(figname * ".png", fig, px_per_unit=2.0)
    end

    # Return the SEBA vector data to the script for future use
    return Σ

end