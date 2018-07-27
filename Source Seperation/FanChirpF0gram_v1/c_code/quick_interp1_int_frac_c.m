%    QUICK_INTERP1_INT_FRAC_C quick linear interpolation
%    y = quick_interp1_int_frac_c(x,pos_int, pos_frac)
%    inputs : x, pos_int (NxM) ,pos_frac (NxM) 
%    outputs : y (NxM) 
%      
%    y(i,j) = x(1+pos_int(i,j)+pos_frac(i,j))  (linear interpolation)
%
%    - x must be a real vector of type: single or double
%    - pos_int must be a real matrix of size NxM and type: uint32
%    - pos_frac must be a real matriz of size NxM and type : single or double
%    - y is a real matrix of size NxM and type : single or double