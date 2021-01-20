function [dx, dy] = diff2(X)
%DIFF2 Calculate gradient of X along the first two dimensions (y and x)
% This function uses a backward finited difference approximation,
% which is different from the builtin gradient function,
% which uses centered differences. The advantage is that for wrapped
% phase data, diff2 returns +/- pi at wrapping discontinuities, whereas
% the builtin gradient function blurs the jump over two pixels 
    sz = size(X);
    sx = sz; sx(2) = 1;
    sy = sz; sy(1) = 1;    
    dx = [zeros(sx, 'like', X), diff(X, 1, 2)];
    dy = [zeros(sy, 'like', X); diff(X, 1, 1)];
end

