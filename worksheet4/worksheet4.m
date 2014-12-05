function worksheet4

    clear all; close all; clc;
    
    % Define the grid and timestep sizes
    Ns = [3 7 15 31];
    dts = [1/64 1/128 1/256 1/512 1/1024 1/2048 1/4096];
    
    % Simulation end time
    t_end = 4/8;
    
    % Specify the times at which we want to take the snapshots and initialize structures to store them
    snapshot_times = [1/8, 2/8, 3/8, 4/8];
    snapshots_implicit = struct([]);
    snapshots_explicit = struct([]);
    
    % Set the folder for storing the images (create if it does not exist)
    [folder, name, ext] = fileparts(which(mfilename));
    images_folder = [folder, '/images'];
    if ~exist(images_folder, 'dir')
        mkdir(images_folder);
    end

    tic;
    
    % DO FOR EACH OF THE GRID SIZES
    for i = 1:length(Ns)

        N = Ns(i);
        Nx = N;
        Ny = N;
        hx = 1/(Nx + 1);
        hy = 1/(Nx + 1);

        % DO FOR EACH OF THE STEP SIZES
        for j = 1:length(dts)
            
            dt = dts(j);
            num_steps = t_end / dt;
            
            disp(['N = ', num2str(N), '    dt = 1/', num2str(1/dt)]);

            % Setup initial and boundary conditions
            T_explicit = zeros(Ny + 2, Nx + 2);
            T_explicit(2:Ny+1, 2:Nx+1) = ones(Ny, Nx);
            T_implicit = T_explicit;
            
            % Define the weights of the neighbors used by Explicit Euler
            mask = [0, dt/hy^2, 0; dt/hx^2, 1-2*dt/hx^2-2*dt/hy^2, dt/hx^2; 0, dt/hy^2, 0];
            
            % DO THE TIMESTEPPING
            for n = 1:num_steps
                
                % Compute T_n+1
                T_explicit(2:end-1, 2:end-1) = filter2(mask, T_explicit, 'valid');
                T_implicit = gauss_seidel(T_implicit, dt, hx, hy);
                
                % Store snapshot at the specified times (snapshot_times)
                timestamp = dt*n;
                if find(timestamp == snapshot_times)
                    snapshots_explicit = add_snapshot(snapshots_explicit, T_explicit, 'Explicit Euler', timestamp, N, dt);
                    snapshots_implicit = add_snapshot(snapshots_implicit, T_implicit, 'Implicit Euler', timestamp, N, dt);
                end
                
            end % Timestepping loop

        end % Stepsize (dt) loop

    end % Grid size (N) loop
    
    toc
    
    % PLOTTING ============================================================
    
    subplot_rows = length(Ns);
    subplot_cols = length(dts);
    
    tic;
    
    for i = 1:length(snapshots_explicit)
        
        disp(['Plotting and saving image ', num2str(i), '/', num2str(length(snapshots_explicit))]);
        
        % Plots for Explicit Euler
        snapshot = snapshots_explicit(i);
        fig_idx = determine_figure_idx(snapshot.t, snapshot_times);
        subplot_idx = determine_subplot_idx(snapshot.N, Ns, snapshot.dt, dts);
        plot_snapshot(snapshot, fig_idx, subplot_rows, subplot_cols, subplot_idx);
        save_snapshot(snapshot, build_image_name(snapshot.method, i, snapshot.N, snapshot.dt, snapshot.t));
        
        % Plots for Implicit Euler
        snapshot = snapshots_implicit(i);
        fig_idx = fig_idx + length(snapshot_times);
        subplot_idx = determine_subplot_idx(snapshot.N, Ns, snapshot.dt, dts);
        plot_snapshot(snapshot, fig_idx, subplot_rows, subplot_cols, subplot_idx);
        save_snapshot(snapshot, build_image_name(snapshot.method, i, snapshot.N, snapshot.dt, snapshot.t));
        
    end
    
    % Set all figures visible again...
    for i = 1:length(snapshot_times)*2
        maximize_figure(figure(i));
    end
    
    toc
    
    % FUNCTIONS ===========================================================
    
    % Function to add a snapshot into a given array of snapshots
    function snapshots_array = add_snapshot(snapshots_array, new_snapshot, method, timestamp, N, dt)
        snapshots_array(end+1).t = timestamp;
        snapshots_array(end).N = N;
        snapshots_array(end).dt = dt;
        snapshots_array(end).data = new_snapshot;
        snapshots_array(end).method = method;
    end
    
    % Function to plot a snapshot in the specified figure and subplot
    function plot_snapshot(snapshot, fig_idx, subplot_rows, subplot_cols, subplot_idx)
        fig = figure(fig_idx);
        set(fig, 'visible', 'off');
        set(fig, 'Name', [snapshot.method, ' (t = ', num2str(snapshot.t), ')'], 'NumberTitle', 'off');
        subplot(subplot_rows, subplot_cols, subplot_idx);
        [x, y] = meshgrid(linspace(0, 1, snapshot.N+2), linspace(0, 1, snapshot.N+2));
        meshc(x, y, snapshot.data);
        title(['N=', num2str(snapshot.N), ' dt=1/', num2str(1/snapshot.dt)], 'FontSize', 8);
    end
    
    % Function to save a snapshot as an image
    function save_snapshot(snapshot, filename)
        temp_fig = figure(100);
        set(temp_fig, 'visible', 'off'); 
        [x, y] = meshgrid(linspace(0, 1, snapshot.N+2), linspace(0, 1, snapshot.N+2));
        surf(x, y, snapshot.data);
        xlabel('X');
        ylabel('Y');
        zlabel('Temperature');
        title([snapshot.method, ' (N=', num2str(snapshot.N), ') (dt=1/', num2str(1/snapshot.dt), ') (t=', num2str(snapshot.t), ')']);
        set(temp_fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 6, 4]);
        saveas(gcf, filename);
        delete(temp_fig);
    end
    
    function figure_idx = determine_figure_idx(t, snapshot_times)
        figure_idx = find(t == snapshot_times);
    end

    function subplot_idx = determine_subplot_idx(N, Ns, dt, dts)
        [row, col] = meshgrid(Ns, dts);
        subplot_idx = find(row == N & col == dt);
    end

    function maximize_figure(fig)
        set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    end

    function image_name = build_image_name(method, idx, N, dt, t)
        image_name = [images_folder, '/', method, ' ', num2str(idx), ' (N=', num2str(N), ' dt=1_', num2str(1/dt), ' t=', num2str(t), ').png'];
    end
    
end