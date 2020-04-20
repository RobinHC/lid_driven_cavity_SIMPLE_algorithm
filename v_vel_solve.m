% v velocity solve function

function [v_vel] = v_vel_solve(v_vel, v_coeff, IVMAX, JVMAX, w_v, dx, pressure)
    Ap_v = v_coeff{1};
    Aw_v = v_coeff{2};
    Ae_v = v_coeff{3};
    As_v = v_coeff{4};
    An_v = v_coeff{5};
    
    % x-sweep
    A = zeros(JVMAX-2, 3);
    
    for i = 2:IVMAX-1
        A(:, 1) = -As_v(i-1, :);
        A(:, 2) = Ap_v(i-1, :) / w_v;
        A(:, 3) = -An_v(i-1, :);
        b = Aw_v(i-1, :).* v_vel(i-1, 2:end-1);
        b = b + Ae_v(i-1, :).* v_vel(i+1, 2:end-1);
        b = b + Ap_v(i-1, :).* v_vel(i, 2:end-1) * ((1-w_v)/w_v);
        b = b + dx * (pressure(i, 2:end-2) - pressure(i, 3:end-1));
        b(1) = b(1) + As_v(i-1, 1) * v_vel(i, 1);
        b(end) = b(end) + An_v(i-1, end) * v_vel(i, end);
        
        v_vel(i, 2:end-1) = TDMA_solver(JVMAX-2, A, b);
    end
    
    % y-sweep
    A = zeros(IVMAX-2, 3);
    
    for j = 2:JVMAX-1
        A(:, 1) = -Aw_v(:, j-1);
        A(:, 2) = Ap_v(:, j-1) / w_v;
        A(:, 3) = -Ae_v(:, j-1);
        b = As_v(:, j-1).* v_vel(2:end-1, j-1);
        b = b + An_v(:, j-1).* v_vel(2:end-1, j+1);
        b = b + Ap_v(:, j-1).* v_vel(2:end-1, j) * ((1-w_v)/w_v);
        b = b + dx * (pressure(2:end-1, j) - pressure(2:end-1, j+1));
        b(1) = b(1) + Aw_v(1, j-1) * v_vel(1, j);
        b(end) = b(end) + Ae_v(end, j-1) * v_vel(end, j);
        
        v_vel(2:end-1, j) = TDMA_solver(IVMAX-2, A, b);
    end
end