% Setup u_vel coefficients function

function [u_coeff] = u_vel_coeff(u_vel, v_vel, dx, dy, rho, gamma, IUMAX, JUMAX)

    IUCV = IUMAX - 2;
    JUCV = JUMAX - 2;
    
    % calculate physical distance delta_x and delta_y
    delta_x = dx * ones(IUCV, JUCV);
    delta_x(1, :) = (3*dx/2);
    delta_x(end, :) = (3*dx/2);
    delta_y = dy * ones(IUCV, JUCV);
    
    % calculate diffusion length del_x and del_y
    del_x_w = dx * ones(IUCV, JUCV);
    del_x_e = dx * ones(IUCV, JUCV);
    del_y_s = dy * ones(IUCV, JUCV);
    del_y_s(:, 1) = dy/2;
    del_y_n = dy * ones(IUCV, JUCV);
    del_y_n(:, end) = dy/2;
    
    % calculate diffusion length array
    Dw = (delta_y * gamma ./ del_x_w);
    De = (delta_y * gamma ./ del_x_e);
    Ds = (delta_x * gamma ./ del_y_s);
    Dn = (delta_x * gamma ./ del_y_n);
    
    % calculate flow strength
    Fw = rho * ((u_vel(1:end-2, 2:end-1)+u_vel(2:end-1, 2:end-1))/2) .* delta_y;
    Fw(1, :) = rho * u_vel(1, 2:end-1) .* delta_y(1, :);
    Fe = rho * ((u_vel(2:end-1, 2:end-1)+u_vel(3:end, 2:end-1))/2) .* delta_y;
    Fe(end, :) = rho * u_vel(end, 2:end-1) .* delta_y(1, :);
    Fs = rho * ((v_vel(2:end-2, 1:end-1)+v_vel(3:end-1, 1:end-1))/2) .* delta_x;
    Fs(1, :) = rho * (((v_vel(1, 1:end-1)+v_vel(2, 1:end-1))/2)*dx/2 + ((v_vel(2, 1:end-1)+v_vel(3, 1:end-1))/2)*dx);
    Fs(end, :) = rho * (((v_vel(end-1, 1:end-1)+v_vel(end, 1:end-1))/2)*dx/2 + ((v_vel(end-2, 1:end-1)+v_vel(end-1, 1:end-1))/2)*dx);
    Fn = rho * ((v_vel(2:end-2, 2:end)+v_vel(3:end-1, 2:end))/2) .* delta_x;
    Fn(1, :) = rho * (((v_vel(1, 2:end)+v_vel(2, 2:end))/2)*dx/2 + ((v_vel(2, 2:end)+v_vel(3, 2:end))/2)*dx);
    Fn(end, :) = rho * (((v_vel(end-1, 2:end)+v_vel(end, 2:end))/2)*dx/2 + ((v_vel(end-2, 2:end)+v_vel(end-1, 2:end))/2)*dx);
    
    % calculate Peclet number
    Pw = Fw .* Dw.^-1;
    Pe = Fe .* De.^-1;
    Ps = Fs .* Ds.^-1;
    Pn = Fn .* Dn.^-1;
        
    % compute coefficient arrays
    Aw_u = Dw.*peclet_function(abs(Pw), "power_law") + max(Fw, 0);
    Ae_u = De.*peclet_function(abs(Pe), "power_law") + max(-Fe, 0);
    As_u = Ds.*peclet_function(abs(Ps), "power_law") + max(Fs, 0);
    An_u = Dn.*peclet_function(abs(Pn), "power_law") + max(-Fn, 0);
    Ap_u = Aw_u + Ae_u + As_u + An_u;
    
    u_coeff = {Ap_u, Aw_u, Ae_u, As_u, An_u};
end