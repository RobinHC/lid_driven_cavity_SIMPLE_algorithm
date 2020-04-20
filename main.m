% take home exam - 2

%% Part - A

% given data
Nx = 5;
Ny = 5;
Re_top = 1000;
Re_bottom = 1000;

% solution
[u_vel, v_vel, pressure, u_res, v_res, p_res] = solve_SIMPLE(Nx, Ny, Re_top, Re_bottom);

% export to excel
filename = 'results_1.xlsx';
writematrix('u_vel', filename, 'Sheet', 1, 'Range', 'A1')
writematrix(u_vel', filename, 'Sheet', 1, 'Range', 'A2')
writematrix('v_vel', filename, 'Sheet', 1, 'Range', 'A9')
writematrix(v_vel', filename, 'Sheet', 1, 'Range', 'A10')
writematrix('pressure', filename, 'Sheet', 1, 'Range', 'A17')
writematrix(pressure', filename, 'Sheet', 1, 'Range', 'A18')

%% Part - B

% given data
L = 0.5;
H = 0.5;
Nx = [16, 32, 64, 128];
Ny = [16, 32, 64, 128];
Re_top = 1000;
Re_bottom = 0;

for i = 1:4
    [u_vel, v_vel, pressure] = solve_SIMPLE(Nx(i), Ny(i), Re_top, Re_bottom);
    
    % centerline velocities
    u_vel_center = (u_vel(Nx(i)/2 + 1, :) + u_vel(Nx(i)/2 + 2, :))/2;
    v_vel_center = (v_vel(:, Ny(i)/2 + 1) + v_vel(:, Ny(i)/2 + 2))/2;
    
    % define axis for plotting
    x_axis = linspace(-L/(2*Nx(i)), L + (L/(2*Nx(i))), Nx(i)+2);
    x_axis(1) = 0;
    x_axis(end) = L;
    y_axis = linspace(-H/(2*Ny(i)), H + (H/(2*Ny(i))), Ny(i)+2);
    y_axis(1) = 0;
    y_axis(end) = H;
    
    figure(1)
    plot(x_axis, u_vel_center)
    hold on
    figure(2)
    plot(y_axis, v_vel_center)
    hold on
end
figure(1)
xlabel('Distance along X [m]')
ylabel('Centerline U Velocity [m/s]')
lgd1 = legend('16X16', '32X32', '64X64', '128X128');
title(lgd1, 'Control Volumes')
hold off

figure(2)
xlabel('Distance from the bottom wall [m]')
ylabel('Centerline V Velocity [m/s]')
lgd1 = legend('16X16', '32X32', '64X64', '128X128');
title(lgd1, 'Control Volumes')
hold off

% expoert the results for comparison
filename = 'exam_2_comparison.xlsx';
writematrix('Y [m]', filename, 'Sheet', 1, 'Range', 'G3')
writematrix('U Velocity [m/s]', filename, 'Sheet', 1, 'Range', 'H3')
writematrix(x_axis', filename, 'Sheet', 1, 'Range', 'G4')
writematrix(u_vel_center', filename, 'Sheet', 1, 'Range', 'H4')
writematrix('X [m]', filename, 'Sheet', 1, 'Range', 'J3')
writematrix('V Velocity [m/s]', filename, 'Sheet', 1, 'Range', 'K3')
writematrix(y_axis', filename, 'Sheet', 1, 'Range', 'J4')
writematrix(v_vel_center, filename, 'Sheet', 1, 'Range', 'K4')
