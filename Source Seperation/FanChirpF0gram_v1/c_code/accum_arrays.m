function accum_arrays
%  ACCUM_ARRAYS accumulation of concatenated arrays elements
%  S = accum_arrays(x,l)
%
%  inputs : x, l 
%  outputs : S
%  
%  x = [x_1 x_2 ... x_m]  where x_i are arrays 
%  of length l_i, and lengths l = [l_1 l_2 ... l_m] 
%  length(x) must be equal to sum(l)
%  The output S is given by:
%               S(i) = sum(x_i)
%
%  x must be a real or complex vector of type: single or double
%  lengths must be a real vector of type: double