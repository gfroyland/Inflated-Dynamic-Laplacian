# This file contains functions used to plot eigenvectors of the inflated Dynamic Laplacian for the switching Double Gyre system defined in AFK24.

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

            # Apply a threshold to the eigenvector data to isolate key coherent set behaviour identified through more extremal values of the eigenvector (positive or negative)
            thresh = 0.25 # The threshold used in AFK24
            plotindices = findall(abs.(eigvec_data) .> thresh)
            eigvec_data = eigvec_data[plotindices]
            traj_pts = traj_pts[plotindices]

            # Set Figure column to the current eigenvector index (1,...,length(inds_to_plot))
            ax = Axis3(fig[1, ind], xlabel=L"$t$", ylabel=L"$x$", zlabel=L"$y$")
            scatterplot = scatter!(ax, traj_pts, color=real.(eigvec_data), colormap=Reverse(:RdBu), colorrange=col_lims, transparency=false, markersize=5) #alpha=0.1
            ax.limits = ((0, 1), (0, 3), (0, 2))
            ax.azimuth = 320π/180
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
        for ind ∈ eachindex(inds_to_plot) # Loop over each of the eigenvectors
            # Retrieve the data for eigenvector inds_to_plot[ind]
            eigvec_data = V[:, inds_to_plot[ind]]
            traj_pts = traj_data_M1

            # Rescale the eigenvector data to have maximum value 1 (or minimum value -1)
            v_max = maximum(abs.(eigvec_data))
            eigvec_data = eigvec_data./v_max
            col_lims = (-1, 1)

            # Apply a threshold to the eigenvector data to isolate key coherent set behaviour identified through more extremal values of the eigenvector (positive or negative)
            thresh = 0.25 # The threshold used in AFK24
            plotindices = findall(abs.(eigvec_data) .> thresh)
            eigvec_data = eigvec_data[plotindices]
            traj_pts = traj_pts[plotindices]

            # Set Figure column to the current eigenvector index (1,...,length(inds_to_plot))
            ax = Axis3(fig[1, ind], xlabel=L"$t$", ylabel=L"$x$", zlabel=L"$y$")
            scatterplot = scatter!(ax, traj_pts, color=real.(eigvec_data), colormap=Reverse(:RdBu), colorrange=col_lims, transparency=false, markersize=5) #alpha=0.1
            ax.limits = ((0, 1), (0, 3), (0, 2))
            ax.azimuth = 320π/180
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

# In this function, we plot select eigenvectors of the inflated Dynamic Laplacian as above, but this time we only plot these eigenvectors for the trajectories corresponding to the top 5% of variance in their eigenvector values over time
# V is an NT × num_of_Λ matrix containing the spacetime eigenvector data for the leading num_of_Λ eigenvectors, given N trajectories and T time steps
# inds_to_plot is a vector listing the indices of the eigenvectors to plot
# traj_data_full is a vector of length T (the number of time steps); each element of this vector is a vector of length N (the number of trajectories), and each element of this vector is a vector of length 2 [x, y] representing the coordinates of trajectory n (n = 1,...,N) at time time_steps(t) (t = 1,...,T).
# time_steps is a vector of length T containing the equispaced time points between 0 and τ (inclusive). 
# figname is the filename (including path) for the Figures produced by this function. The Figures will be saved as pngs (with the file extension to be added in this function), with one version of the eigenvector Figures prepared for 𝕄₀ and another for 𝕄₁ (see below for definitions of 𝕄₀ and 𝕄₁)
function plot_eigvecs_topvar(V, inds_to_plot, traj_data_full, time_steps, figname)

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

            # Find the indices of the trajectories corresponding to the top 5% of variances of their eigenvector values over time
            eigvec_data_mat = reshape(eigvec_data,(N,T)) # Rearrange the eigenvector data into an N × T matrix to make the variance calculations easier
            eigvec_data_var = var(eigvec_data_mat,dims=2)[:,1] # Take the variances of the eigenvector values for each trajectory over time
            num_trajs_to_plot = round(Int64, 0.05 * N) # Calculate 5% of the total number of trajectories
            traj_indices = sort(partialsortperm(eigvec_data_var, 1:num_trajs_to_plot, rev=true)) # Find the indices corresponding to the num_trajs_to_plot largest variances, and sort these so that the indices are increasing
            
            # Recover the eigenvector and trajectory data corresponding to the traj_indices found above
            eigvec_data = reshape(eigvec_data_mat[traj_indices,:],T*num_trajs_to_plot)
            traj_pts_mat = reshape(traj_pts,(N,T))
            traj_pts = reshape(traj_pts_mat[traj_indices,:],T*num_trajs_to_plot)

            # Set Figure column to the current eigenvector index (1,...,length(inds_to_plot))
            ax = Axis3(fig[1, ind], xlabel=L"$t$", ylabel=L"$x$", zlabel=L"$y$")
            scatterplot = scatter!(ax, traj_pts, color=real.(eigvec_data), colormap=Reverse(:RdBu), colorrange=col_lims, transparency=false, markersize=5) #alpha=0.1
            ax.limits = ((0, 1), (0, 3), (0, 2))
            ax.azimuth = 320π/180
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
        for ind ∈ eachindex(inds_to_plot) # Loop over each of the eigenvectors
            # Retrieve the data for eigenvector inds_to_plot[ind]
            eigvec_data = V[:, inds_to_plot[ind]]
            traj_pts = traj_data_M1

            # Rescale the eigenvector data to have maximum value 1 (or minimum value -1)
            v_max = maximum(abs.(eigvec_data))
            eigvec_data = eigvec_data./v_max
            col_lims = (-1, 1)

            # Find the indices of the trajectories corresponding to the top 5% of variances of their eigenvector values over time
            eigvec_data_mat = reshape(eigvec_data,(N,T)) # Rearrange the eigenvector data into an N × T matrix to make the variance calculations easier
            eigvec_data_var = var(eigvec_data_mat,dims=2)[:,1] # Take the variances of the eigenvector values for each trajectory over time
            num_trajs_to_plot = round(Int64, 0.05 * N) # Calculate 5% of the total number of trajectories
            traj_indices = sort(partialsortperm(eigvec_data_var, 1:num_trajs_to_plot, rev=true)) # Find the indices corresponding to the num_trajs_to_plot largest variances, and sort these so that the indices are increasing

            # Recover the eigenvector and trajectory data corresponding to the traj_indices found above
            eigvec_data = reshape(eigvec_data_mat[traj_indices,:],T*num_trajs_to_plot)
            traj_pts_mat = reshape(traj_pts,(N,T))
            traj_pts = reshape(traj_pts_mat[traj_indices,:],T*num_trajs_to_plot)

            # Set Figure column to the current eigenvector index (1,...,length(inds_to_plot))
            ax = Axis3(fig[1, ind], xlabel=L"$t$", ylabel=L"$x$", zlabel=L"$y$")
            scatterplot = scatter!(ax, traj_pts, color=real.(eigvec_data), colormap=Reverse(:RdBu), colorrange=col_lims, transparency=false, markersize=5) #alpha=0.1
            ax.limits = ((0, 1), (0, 3), (0, 2))
            ax.azimuth = 320π/180
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