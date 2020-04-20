% pressure solve function

function [pressure, p_corr] = p_correction(pressure, p_corr, p_coeff, IPMAX, JPMAX, dy, dx, w_p, rho, u_vel, v_vel)
    Ap_p = p_coeff{1};
    Aw_p = p_coeff{2};
    Ae_p = p_coeff{3};
    As_p = p_coeff{4};
    An_p = p_coeff{5};
    
    % x-sweep
    A = zeros(JPMAX-2, 3);

    for i = 2:IPMAX-1
        A(:, 1) = -As_p(i-1, :);
        A(:, 2) = Ap_p(i-1, :) / w_p;
        A(:, 3) = -An_p(i-1, :);
        b = Aw_p(i-1, :).* p_corr(i-1, 2:end-1);
        b = b + Ae_p(i-1, :).* p_corr(i+1, 2:end-1);
        b = b + Ap_p(i-1, :).* p_corr(i, 2:end-1) * ((1-w_p)/w_p);
        b = b + rho * (dy * (u_vel(i-1, 2:end-1) - u_vel(i, 2:end-1)) + dx * (v_vel(i, 1:end-1) - v_vel(i, 2:end)));
        b(1) = b(1) + As_p(i-1, 1) * p_corr(i, 1);
        b(end) = b(end) + An_p(i-1, end) * p_corr(i, end);
        
        p_corr(i, 2:end-1) = TDMA_solver(JPMAX-2, A, b);
    end
    
    % y-sweep
    A = zeros(IPMAX-2, 3);
    
    for j = 2:JPMAX-1
        A(:, 1) = -Aw_p(:, j-1);
        A(:, 2) = Ap_p(:, j-1) / w_p;
        A(:, 3) = -Ae_p(:, j-1);
        b = As_p(:, j-1).* p_corr(2:end-1, j-1);
        b = b + An_p(:, j-1).* p_corr(2:end-1, j+1);
        b = b + Ap_p(:, j-1).* p_corr(2:end-1, j) * ((1-w_p)/w_p);
        b = b + rho * (dy * (u_vel(1:end-1, j) - u_vel(2:end, j)) + dx * (v_vel(2:end-1, j-1) - v_vel(2:end-1, j)));
        b(1) = b(1) + Aw_p(1, j-1) * p_corr(1, j);
        b(end) = b(end) + Ae_p(end, j-1) * p_corr(end, j);
        
        p_corr(2:end-1, j) = TDMA_solver(IPMAX-2, A, b);
    end
    
    % calculating pressure values
    pressure = pressure + 0.8 * p_corr;
    
    % setting up boundary values for pressure field
    pressure(1, 2:end-1) = pressure(2, 2:end-1);
    pressure(end, 2:end-1) = pressure(end-1, 2:end-1);
    pressure(2:end-1, 1) = pressure(2:end-1, 2);
    pressure(2:end-1, end) = pressure(2:end-1, end-1);
end