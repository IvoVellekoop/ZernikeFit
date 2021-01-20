function best_a = zernike_fit(E, Z, U, m, a_initial)
%FIT_GRADIENT Fits the phase of field E to a set of modes Z, taking into
%             account mask m.
%
%   Since we don't know the absolute value phi (only the wrapped phase),
%   we cannot directly decompose phi into modes Z. However, we can recover
%   the gradient of phi with reasonable accuracy. Therefore, we decompose
%   the gradient of phi into a superposition of gradients of modes Z.
%   High values of dphi are very sensitive to noise (causing e.g phase
%   jumps in the wrong direction). Therefore, the fitting is performed
%   iteratively: first, a low pass filtered version of the field is fitted.
%   Then, this procedure is repeated iteratively for the residue (i.e. the phase 
%   difference between the fit and the actual, unfiltered field), choosing
%   a weaker low-pass filter every step.
%   INPUT: 
%       E   Complex values of the field to fit the wavefront to
%       Z   Matrix of Zernike modes as returned by zernfun_cart
%
%   TODO: implement weighted linear least squares to assign a lower weight
%   to pixels with lower intensity
%   TODO: find automatic way to choose sequence of low pass filters or
%   alllow the user to configure it.
%
%  
   N = size(U, 1); % number of modes to fit
   if nargin > 4
       % If starting values for a are provided, use those
       a = a_initial(:);
       if length(a) < N
           a(N) = 0;
       else
           a = a(1:N);
       end
   else
       a = zeros(N, 1);           
   end
   max_merit = 0;
   residue = E;
   E_ap = m.pack(abs(E)); %absolute value of field of pixels that are inside the mask aperture (used as weight in metric)
   
   for it=5:-0.2:0
        % low-pass filter the residue
        residue = imgaussfilt(real(residue), 2^it) + 1i * imgaussfilt(imag(residue), 2^it); 

        % calculate phase gradient and fit it with modes Z
        [dphix, dphiy] = diff2(angle(residue));
        dphix = m.apply(wrap(dphix));
        dphiy = m.apply(wrap(dphiy));
        dphi = m.pack(dphix, dphiy);
        a = a + U * dphi;

        % update the reconstruction and calculate the residue again
        rec = sum(Z .* reshape(a, 1, 1, N), 3);
        residue = E .* exp(-1i * rec);

        %% Keep best values of a. Leave out normalization of merit for performance reasons
        E_rec = m.pack(residue);
        inp = E_rec' * E_ap;
        merit = abs(inp);
        if merit > max_merit                
            best_a = a;
            best_a(1) = -angle(inp);
            max_merit = merit;
        end
        %% debugging only: show the current residue
        %imagesc(angle(residue)); title(['residue. |C|^2= ' num2str(merit)]); pause(0.1);
   end
end

function phi = wrap(phi)
    %% WRAP(phi) wraps phase phi to the range [-pi, +pi)
    phi = mod(phi + pi, 2*pi) - pi;
end
