% v velocity correction function

function [v_vel] = v_vel_correction(v_vel, Ap_v, p_corr, dx)

    v_vel_corr = (dx ./ Ap_v) .* (p_corr(2:end-1, 2:end-2) - p_corr(2:end-1, 3:end-1));
    v_vel(2:end-1, 2:end-1) = v_vel(2:end-1, 2:end-1) + v_vel_corr;
end