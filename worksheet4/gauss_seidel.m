function T = gauss_seidel(T_prev, dt, hx, hy)

    T = T_prev;
    num_rows = size(T, 1);
    num_cols = size(T, 2);
    
    k = (1 + 2*dt/hx^2 + 2*dt/hy^2);
    
    residual_norm = Inf;
    iteration = 0;
    
    while residual_norm > 1e-6
    
        for j = 2:num_rows-1
           for i = 2:num_cols-1
               T(j, i) = (((T(j, i-1)+T(j, i+1))/hx^2 + (T(j-1, i)+T(j+1, i))/hy^2)*dt + T_prev(j, i))/k;
           end
        end

        residual = 0;        
        for j = 2:num_rows-1
           for i = 2:num_cols-1
               residual = residual + ((T(j, i-1) - 2*T(j, i) + T(j, i+1))/hx^2 + (T(j-1, i) - 2*T(j, i) + T(j+1, i))/hy^2 - (T(j, i) - T_prev(j, i))/dt)^2;
           end
        end
        residual_norm = sqrt(residual/numel(T(2:end-1,2:end-1)));
        
        iteration = iteration + 1;
        
    end
%     disp(iteration);

end