% u velocity correction function

function [u_vel] = u_vel_correction(u_vel, Ap_u, p_corr, dy)

    u_vel_corr = (dy ./ Ap_u) .* (p_corr(2:end-2, 2:end-1) - p_corr(3:end-1, 2:end-1));
    u_vel(2:end-1, 2:end-1) = u_vel(2:end-1, 2:end-1) + u_vel_corr;
end