%% Demo script for zernike_fit
%
% This code generates a synthetic wavefront, adds noise, and
% fits the wavefront using the zernike_fit method
%
% Four steps are needed to perform the fit:
% 1) A region of interest is defined, using a 'Mask' object
% 2) A set of Zernike modes to be fitted is generated
% 3) A model matrix is constructed from the set of modes
% 4) The actual fitting is performed by calling zernike_fit
%
% Constructing the model matrix can be relatively slow if many modes are
% fitted. Therefore, the model matrix is pre-computed once, and used in
% all subsequent calls to zernike_fit.

%% Fitting options
% Define a coordinate system and region of interest by means of a
% 'mask' object. The mask object can be used to center and crop 
% experimental data. For synthetic data, we just use an ROI ranging 
% from 1 to 1000 in both dimensions.
% 
m = Mask([1, 1, 1000, 1000], 'Shape', 'circular', 'Smoothness', 0.3);  % coordinate system
N = 28; % number of zernike modes
noise_level = 0.5; % noise level to use in simulations

%% Construct the set of basis functions and fitting model
orders = zernike_order(N); % use up to order 28 Zernike mode
Z = zernfun_cart(m.x, m.y, [orders.n], [orders.m], false);

%% Construct synthetic data
a = randn(1, 1, N)*2-1; % random coefficients
phi = sum(Z .* a, 3);
noise = noise_level * sqrt(0.5) * (randn(size(phi)) + 1i * randn(size(phi)));
E = m.filter .* (noise+exp(1i * phi));    % smoothen the edges of the measured data

%% fit using the gradient-fit method
U = model(Z, m);
a = zernike_fit(E, Z, U, m);

%% calculate quality of fit, only use the pixels that are inside the ROI
E_rec = m.pack(abs(E)) .* exp(1i * m.pack(Z) * a);
E_sp = m.pack(E);
c2 = (abs(E_rec' * E_sp) / norm(E_rec) / norm(E_sp))^2;

%% print in the command window
fprintf('Quality of fit (correlation coefficient):\t%d\n', c2);

%% plot C2, original phase and residue
figure(1); imagesc(angle(E));
title('original'); colorbar; axis image;

figure(2); imagesc(angle(m.unpack(E_rec))); 
title(['reconstructed (', num2str(N), 'modes)']); 
colorbar; axis image;

figure(3); imagesc(angle(E .* conj(m.unpack(E_rec))));
title('residue'); 
colorbar; axis image;

