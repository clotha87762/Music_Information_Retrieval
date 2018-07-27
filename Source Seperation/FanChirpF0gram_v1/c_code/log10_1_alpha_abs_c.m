function log10_1_alpha_abs_c
%  LOG10_1_ALPHA_ABS_c calculates log10(1+abs(\alpha*x)) efficiently
%     y = log10_1_alpha_abs_c(x,alpha)
%
%     input: x, alpha
%     output: y
%
%     y = log10(1+abs(\alpha*x)) 
%
%     - x can be real or complex and must be of type: double
%     - alpha should be a scalar of type: double