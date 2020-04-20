% Setup v_vel coefficients function

function [v_coeff] = v_vel_coeff(u_vel, v_vel, dx, dy, rho, gamma, IVMAX, JVMAX)

    IVCV = IVMAX - 2;
    JVCV = JVMAX - 2;
    
    % calculate physical distance delta_x and delta_y
    delta_x = dx * ones(IVCV, JVCV);
    delta_y = dy * ones(IVCV, JVCV);
    delta_y(:, 1) = 3*dy/2;
    delta_y(:, end) = 3*dy/2;
    
    % calculate diffusion length del_x and del_y
    del_x_w = dx * ones(IVCV, JVCV);
    del_x_w(1, :) = dx/2;
    del_x_e = dx * ones(IVCV, JVCV);
    del_x_e(end, :) = dx/2;
    del_y_s = dy * ones(IVCV, JVCV);
    del_y_n = dy * ones(IVCV, JVCV);
    
    % calculate diffusion length array
    Dw = (delta_y * gamma ./ del_x_w);
    De = (delta_y * gamma ./ del_x_e);
    Ds = (delta_x * gamma ./ del_y_s);
    Dn = (delta_x * gamma ./ del_y_n);
    
    % calculate flow strength
    Fw = rho * ((u_vel(1:end-1, 2:end-2)+u_vel(1:end-1, 3:end-1))/2) .* delta_y;
    Fw(:, 1) = rho * (((u_vel(1:end-1, 1)+u_vel(1:end-1, 2))/2)*dy/2 + ((u_vel(1:end-1, 2)+u_vel(1:end-1, 3))/2)*dy);
    Fw(:, end) = rho * (((u_vel(1:end-1, end-2)+u_vel(1:end-1, end-1))/2)*dy + ((u_vel(1:end-1, end-1)+u_vel(1:end-1, end))/2)*dy/2);
    Fe = rho * ((u_vel(2:end, 2:end-2)+u_vel(2:end, 3:end-1))/2) .* delta_y;
    Fe(:, 1) = rho * (((u_vel(2:end, 1)+u_vel(2:end, 2))/2)*dy/2 + ((u_vel(2:end, 2)+u_vel(2:end, 3))/2)*dy);
    Fe(:, end) = rho * (((u_vel(2:end, end-2)+u_vel(2:end, end-1))/2)*dy + ((u_vel(2:end, end-1)+u_vel(2:end, end))/2)*dy/2);
    Fs = rho * ((v_vel(2:end-1, 1:end-2)+v_vel(2:end-1, 2:end-1))/2) .* delta_x;
    Fs(:, 1) = rho * v_vel(2:end-1, 1) .* delta_x(:, 1);
    Fn = rho * ((v_vel(2:end-1, 2:end-1)+v_vel(2:end-1, 3:end))/2) .* delta_x;
    Fn(:, end) = rho * v_vel(2:end-1, end) .* delta_x(:, end);
    
    % calculate Peclet number
    Pw = Fw .* Dw.^-1;
    Pe = Fe .* De.^-1;
    Ps = Fs .* Ds.^-1;
    Pn = Fn .* Dn.^-1;
        
    % compute coefficient arrays
    Aw_v = Dw.*peclet_function(abs(Pw), "power_law") + max(Fw, 0);
    Ae_v = De.*peclet_function(abs(Pe), "power_law") + max(-Fe, 0);
    As_v = Ds.*peclet_function(abs(Ps), "power_law") + max(Fs, 0);
    An_v = Dn.*peclet_function(abs(Pn), "power_law") + max(-Fn, 0);
    Ap_v = Aw_v + Ae_v + As_v + An_v;
    
    v_coeff = {Ap_v, Aw_v, Ae_v, As_v, An_v};
end