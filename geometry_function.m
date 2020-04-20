% Geometry function

function [IPCV, JPCV, dx, dy, L, H] = geometry_function(Nx, Ny)
    L = 0.5;
    H = 0.5;

    IPCV = Nx;
    JPCV = Ny;
    
    dx = L / IPCV;
    dy = H / JPCV;
end
    