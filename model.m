function U = model(Z, m)
% Purpose: pre-compute model matrix that is used in zernike_fit
% INPUT:
% Z = Zernike polynomials built using zernfun_cart.m. 
%     Note: the first entry must be the piston distortion, which is treated
%     differently.
% m = Mask object. Typically, this is the same mask object as was used to
%     generate the Zernike polynomials. The mask object is used to pack
%     all pixels inside the unit disk into a single column vector. Pixels
%     outside the unit disk are not considered in the fit.
%
% OUTPUT: model matrix U for use in zernike_fit
    Z = m.crop(Z);
    [dZx, dZy] = diff2(Z);
    dZ = m.pack(m.apply(dZx), m.apply(dZy));
    U = pinv(dZ);
end
