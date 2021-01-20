function modes = zernike_order(N_modes)
% Author: Giulia Sereni
% Computes radial and azimuthal order of the first 'N_modes' Zernike modes.
% (see at 'zernfun_cart' function for a detailed description of the Zernike functions).
%
% INPUT
% N_modes = numbers of modes. e.g. N_modes = 2 generates the first two (n,m)
% couples. The corresponding names of these modes can be found in
% 'zernike_pyramid' function.
%
% OUTPUT
% modes = structure containing n and m. n is the radial order of the 
% zernike function, while m is the azimuthal order.

    mode.n = 0;
    mode.m = 0;
    modes(N_modes) = mode; % preallocate array of Nmodes mode structures

    for j = 1:N_modes                     
        modes(j).n = ceil(-1.5+0.5*sqrt(1+j*8));              %radial order
    end

    for j = 0:N_modes-1
        n = modes(j+1).n;
        modes(j+1).m = n-j+(n+1)*n/2;    %azimuthal order
    end
end