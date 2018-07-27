function iir_cqt_c
%  IIR_CQT_C Infinite Impulse Response Constant Q Transform
%  y = IIR_CQT_C(x, p, prepad, pospad)
%     input: x, p, prepad , pospad
%     output: y
%             
%     y = LTV filtering of x given the poles in p
%
%     - length(p) = length(x) > prepad, pospad
%     - x should be complex and of type double
%     - p should be real and of type double
%     - prepad , pospad are scalars of type double
%
% See article:
%  P. Cancela, M. Rocamora, and E. López, “An efficient multi-resolution
%  spectral transform for music analysis,” in International Society for
%  Music Information Retrieval Conference, 10th. ISMIR 2009. Kobe,
%  Japan, pp. 309–314, 2009. for details on the IIR-CQT