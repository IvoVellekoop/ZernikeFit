# ZernikeFit
MATLAB script to fit Zernike modes to phase-wrapped data
(c) Giulia Sereni and Ivo Vellekoop, University of Twente

See demo.m for an example of how to use the zernike_fit script.

~~~
demo.m  - demo script
diff2.m - function to compute numerical gradient of a 2-d image
Mask.m  - mask object that contains information about coordinates and
          region of interest. Can be used to select what region of the
          input data is used in the fit.
          Also contains function to pack pixels to a vector (keeping only
          the pixels inside the circular ROI) or to reverse this operation
          (unpack).
model.m - prepares a model matrix for use in the zernike_fit function
zernfun_cart.m  - computes zernike modes in a 2-D cartesian basis
zernike_fit     - the function that performs the actual fitting algorithm
zernike_mode    - computes radial and azimuthal order of the first N Zernike modes
~~~

