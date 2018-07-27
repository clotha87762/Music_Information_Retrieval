function recursive_octave_accum
% RECURSIVE_OCTAVE_ACCUM performs a periodic acumulation sum of the columms
% of a matrix
%
% y = recursive_octave_accum(x,period)
% inputs : x,n_octave 
% outputs : y 
%       
%   y(i,:) = \sum_k x(i+k*period,:)
%
%  - x must be a real matrix of type: double
%  - n_octave be a real scalar of type: double