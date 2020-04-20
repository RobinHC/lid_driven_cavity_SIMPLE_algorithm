% 2-D SIMPLE algorithm

function [u_vel, v_vel, pressure, u_res, v_res, p_res] = solve_SIMPLE(Nx, Ny, Re_top, Re_bottom)
    
    % call geometry function
    [IPCV, JPCV, dx, dy, L, H] = geometry_function(Nx, Ny);
    IPMAX = IPCV + 2;
    JPMAX = JPCV + 2;
    pressure = zeros(IPMAX, JPMAX);
    p_corr = zeros(IPMAX, JPMAX);
    IUMAX = IPCV + 1;
    JUMAX = JPCV + 2;
    u_vel = zeros(IUMAX, JUMAX);
    IVMAX = IPCV + 2;
    JVMAX = JPCV + 1;
    v_vel = zeros(IVMAX, JVMAX);
    
    % call property funciton
    [rho, k, mu, Cp] = property_function();
    
    % call initial conditions function
    u_top_ref = Re_top * mu / (rho * L);
    u_bottom_ref = Re_bottom * mu / (rho * L);
    [u_vel] = ics_function(u_vel, u_top_ref, u_bottom_ref);
    
    iter = 1;
    u_res = ones(800, 1);
    v_res = ones(800, 1);
    p_res = ones(800, 1);
    while (u_res(iter) > 1e-6 && v_res(iter) > 1e-6 && p_res(iter) > 1e-5)
        % call u_vel coefficient function
        u_coeff = u_vel_coeff(u_vel, v_vel, dx, dy, rho, mu, IUMAX, JUMAX);

        % call u_vel solve function
        w_u = 0.5;
        u_vel = u_vel_solve(u_vel, u_coeff, IUMAX, JUMAX, w_u, dy, pressure);

        % call v_vel coefficient function
        v_coeff = v_vel_coeff(u_vel, v_vel, dx, dy, rho, mu, IVMAX, JVMAX);

        % call v_vel solve function
        w_v = 0.5;
        v_vel = v_vel_solve(v_vel, v_coeff, IVMAX, JVMAX, w_v, dx, pressure);

        % call pressure correction coefficient function
        p_coeff = p_corr_coeff(dx, dy, rho, u_coeff{1}, v_coeff{1}, IPMAX, JPMAX);

        % call pressure correction solve function
        w_p = 1;
        [pressure, p_corr] = p_correction(pressure, p_corr, p_coeff, IPMAX, JPMAX, dy, dx, w_p, rho, u_vel, v_vel);

        % call u_vel correction function
        u_vel = u_vel_correction(u_vel, u_coeff{1}, p_corr, dy);

        % call v_vel correction function
        v_vel = v_vel_correction(v_vel, v_coeff{1}, p_corr, dx);
        
        p_corr = zeros(size(p_corr));

        iter = iter + 1;

        % call convergence check function
        [u_res(iter), v_res(iter), p_res(iter)] = convergence(u_vel, u_coeff, v_vel, v_coeff, pressure, u_top_ref, dx, dy, rho, L);
    end 
end