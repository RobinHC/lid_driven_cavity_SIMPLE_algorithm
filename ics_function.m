% Initial condition setup function

function [u_vel] = ics_function(u_vel, u_top_ref, u_bottom_ref)
    u_vel(:, 1) = u_bottom_ref;
    u_vel(:, end) = u_top_ref;
    u_vel(2:end-1, 2:end-1) = 0.0001;
end
    