function worksheet4

    clear all; close all; clc;
    
    % Define the grid and timestep sizes
    Ns = [3 7 15 31];
    dts = [1/64 1/128 1/256 1/512 1/1024 1/2048 1/4096];
    
    % Simulation end time
    t_end = 4/8;
    
    % Specify the times at which we want to take the snapshots
    % and initialize structures to store them
    snapshot_times = [1/8, 2/8, 3/8, 4/8];
    snapshots_implicit = struct([]);
    snapshots_explicit = struct([]);

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
                    snapshots_explicit = add_snapshot(snapshots_explicit, T_explicit, timestamp, N, dt);
                    snapshots_implicit = add_snapshot(snapshots_implicit, T_implicit, timestamp, N, dt);
                end
                
            end % Timestepping loop

        end % Stepsize (dt) loop

    end % Grid size (N) loop
    
    toc
    
    % Now do the plotting
    
    disp('Plotting...');
    subplot_rows = length(Ns);
    subplot_cols = length(dts);
    
    tic;
    
    for i = 1:length(snapshots_explicit)
        
        % Plots for Explicit Euler
        snapshot = snapshots_explicit(i);
        fig_idx = determine_figure_idx(snapshot.t, snapshot_times);
        subplot_idx = determine_subplot_idx(snapshot.N, Ns, snapshot.dt, dts);
        fig = figure(fig_idx);
        set(fig, 'visible', 'off');
        set_figure_title(fig, 'Explicit Euler', snapshot.t);
        plot_snapshot(snapshot, subplot_rows, subplot_cols, subplot_idx);
        
        % Plots for Implicit Euler
        snapshot = snapshots_implicit(i);
        fig = figure(fig_idx + length(snapshot_times));
        subplot_idx = determine_subplot_idx(snapshot.N, Ns, snapshot.dt, dts);
        set(fig, 'visible', 'off');
        set_figure_title(fig, 'Implicit Euler', snapshot.t);  
        plot_snapshot(snapshot, subplot_rows, subplot_cols, subplot_idx);
        
    end
    
    % Set all figures visible again...
    for i = 1:length(snapshot_times)*2
        fig = figure(i);
        maximize_figure(fig);
    end
    
    toc
    
    % Function to add a snapshot into a given array of snapshots
    function snapshots_array = add_snapshot(snapshots_array, new_snapshot, timestamp, N, dt)
        snapshots_array(end+1).t = timestamp;
        snapshots_array(end).N = N;
        snapshots_array(end).dt = dt;
        snapshots_array(end).data = new_snapshot;
    end
    
    % Function to plot a snapshot in the active figure
    % Subplot index is determined based on the values of N and dt
    function plot_snapshot(snapshot, subplot_rows, subplot_cols, subplot_idx)
        subplot(subplot_rows, subplot_cols, subplot_idx);
        meshc(snapshot.data);
        title(['N=', num2str(snapshot.N), ' dt=1/', num2str(1/snapshot.dt)], 'FontSize', 8);
    end
    
    % Function to save a snapshot as an image
    function save_snapshot(snapshot)
        temp_fig = figure(100);
        set(temp_fig, 'visible', 'off'); 
        surf(snapshot.data);
        title(['N=', num2str(snapshot.N), ' dt=1/', num2str(1/snapshot.dt)]);
        saveas(gcf, ['implicit_N(', num2str(snapshot.N), ')_dt(', num2str(snapshot.dt), ')_t(', num2str(snapshot.t), ').png']);    
        delete(temp_fig);
    end
    
    function figure_idx = determine_figure_idx(t, snapshot_times)
        figure_idx = find(t == snapshot_times);
    end

    function subplot_idx = determine_subplot_idx(N, Ns, dt, dts)
        subplot_row = find(N == Ns);
        subplot_col = find(dt == dts);
        subplot_idx = subplot_col + (subplot_row-1) * length(dts);
    end

    function set_figure_title(fig, method_name, timestamp)
        set(fig, 'Name', [method_name, ' (t = ', num2str(timestamp), ')'], 'NumberTitle', 'off');
    end

    function maximize_figure(fig)
        set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    end
    
end