% residuals calculation function

function [u_res, v_res, p_res] = convergence(u_vel, u_coeff, v_vel, v_coeff, pressure, u_ref, dx, dy, rho, L)
    Ap_u = u_coeff{1};
    Aw_u = u_coeff{2};
    Ae_u = u_coeff{3};
    As_u = u_coeff{4};
    An_u = u_coeff{5};

    Ap_v = v_coeff{1};
    Aw_v = v_coeff{2};
    Ae_v = v_coeff{3};
    As_v = v_coeff{4};
    An_v = v_coeff{5};
    
    % u velocity residual calculation
    numerator = Ap_u.*u_vel(2:end-1, 2:end-1);
    numerator = numerator - (Aw_u.*u_vel(1:end-2, 2:end-1) + Ae_u.*u_vel(3:end, 2:end-1) + As_u.*u_vel(2:end-1, 1:end-2) + An_u.*u_vel(2:end-1, 3:end));
    numerator = numerator - dx * (pressure(2:end-2, 2:end-1) - pressure(3:end-1, 2:end-1));
    numerator = sum(sum(abs(numerator)));
    denom = sum(sum(abs(Ap_u.*u_vel(2:end-1, 2:end-1))));
    u_res = numerator / denom * 100;
    
    % v velocity residual calculation
    numerator = Ap_v.*v_vel(2:end-1, 2:end-1);
    numerator = numerator - (Aw_v.*v_vel(1:end-2, 2:end-1) + Ae_v.*v_vel(3:end, 2:end-1) + As_v.*v_vel(2:end-1, 1:end-2) + An_v.*v_vel(2:end-1, 3:end));
    numerator = numerator - dy*(pressure(2:end-1, 2:end-2) - pressure(2:end-1, 3:end-1));
    numerator = sum(sum(abs(numerator)));
    denom = sum(sum(abs(Ap_v.*v_vel(2:end-1, 2:end-1))));
    v_res = numerator / denom * 100;
    
    % pressure residual calculation
    numerator = rho * (u_vel(1:end-1, 2:end-1) - u_vel(2:end, 2:end-1)) * dy;
    numerator = numerator + rho * (v_vel(2:end-1, 1:end-1) - v_vel(2:end-1, 2:end)) * dx;
    numerator = sum(sum(abs(numerator)));
    denom = rho * u_ref * L;
    p_res = numerator / denom * 1000;
end