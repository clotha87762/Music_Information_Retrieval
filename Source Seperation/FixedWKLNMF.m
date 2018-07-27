function [ H ] = FixedWKLNMF(V , W  )
%FIXEDWKLNMF Summary of this function goes here
%   Detailed explanation goes here
F = size(V,1);
T = size(V,2);
rand('seed',0)
K = size(W,2);
H = 1+rand(K, T);
ONES = ones(F,T);
for i=1:100
% update activations
H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
% update dictionaries
end
% normalize W to sum to 1

%sumW = sum(W);
%H = diag(sumW)*H;

end

