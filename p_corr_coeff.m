% Setup pressure correction coefficients function

function [p_coeff] = p_corr_coeff(dx, dy, rho, Ap_u, Ap_v, IPMAX, JPMAX)
    
    % calculate d_u and d_v
    d_u = dy ./ Ap_u;
    d_v = dx ./ Ap_v;

    % calculate a coefficients
    Aw_p = zeros(IPMAX-2, JPMAX-2);
    Ae_p = zeros(IPMAX-2, JPMAX-2);
    As_p = zeros(IPMAX-2, JPMAX-2);
    An_p = zeros(IPMAX-2, JPMAX-2);
    Aw_p(1:end-1, :) = rho * d_u * dy;
    Ae_p(2:end, :) = rho * d_u * dy;
    As_p(:, 2:end) = rho * d_v * dx;
    An_p(:, 1:end-1) = rho * d_v * dx;
    
    Ap_p = Aw_p + Ae_p + As_p + An_p;
    
    p_coeff = {Ap_p, Aw_p, Ae_p, As_p, An_p};
end