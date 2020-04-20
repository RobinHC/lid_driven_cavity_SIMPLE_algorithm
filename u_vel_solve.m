% u velocity solve function

function [u_vel] = u_vel_solve(u_vel, u_coeff, IUMAX, JUMAX, w_u, dy, pressure)
    Ap_u = u_coeff{1};
    Aw_u = u_coeff{2};
    Ae_u = u_coeff{3};
    As_u = u_coeff{4};
    An_u = u_coeff{5};
    
    % x-sweep
    A = zeros(JUMAX-2, 3);
    
    for i = 2:IUMAX-1
        A(:, 1) = -As_u(i-1, :);
        A(:, 2) = Ap_u(i-1, :) / w_u;
        A(:, 3) = -An_u(i-1, :);
        b = Aw_u(i-1, :).* u_vel(i-1, 2:end-1);
        b = b + Ae_u(i-1, :).* u_vel(i+1, 2:end-1);
        b = b + Ap_u(i-1, :).* u_vel(i, 2:end-1) * ((1-w_u)/w_u);
        b = b + dy * (pressure(i, 2:end-1) - pressure(i+1, 2:end-1));
        b(1) = b(1) + As_u(i-1, 1) * u_vel(i, 1);
        b(end) = b(end) + An_u(i-1, end) * u_vel(i, end);
        
        u_vel(i, 2:end-1) = TDMA_solver(JUMAX-2, A, b);
    end
    
    % y-sweep
    A = zeros(IUMAX-2, 3);
    
    for j = 2:JUMAX-1
        A(:, 1) = -Aw_u(:, j-1);
        A(:, 2) = Ap_u(:, j-1) / w_u;
        A(:, 3) = -Ae_u(:, j-1);
        b = As_u(:, j-1).* u_vel(2:end-1, j-1);
        b = b + An_u(:, j-1).* u_vel(2:end-1, j+1);
        b = b + Ap_u(:, j-1).* u_vel(2:end-1, j) * ((1-w_u)/w_u);
        b = b + dy * (pressure(2:end-2, j) - pressure(3:end-1, j));
        b(1) = b(1) + Aw_u(1, j-1) * u_vel(1, j);
        b(end) = b(end) + Ae_u(end, j-1) * u_vel(end, j);
        
        u_vel(2:end-1, j) = TDMA_solver(IUMAX-2, A, b);
    end
end