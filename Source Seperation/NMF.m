function [ W , H ] = NMF( V , R , epsilon )
%NMF Summary of this function goes here
%   Detailed explanation goes here
 W = rand(size(V,1),R);H = rand(R,size(V,2)); WW = W; HH = H;
 deltaW = 99999999; deltaH = 99999999;delta = 99999999;a = 0;
 while deltaH>epsilon || deltaW>epsilon,
    W = WW;
    H = HH;
    HH = H.*(((W)'*V)./ ((W)'*(W*H+eps)));
    WW = W.*((V*(HH)')./((W*HH+eps)*(HH)'));
    deltaW = norm(WW - W);
    deltaH = norm(HH - H);
 end
    W = WW;
    H = HH;
    
end

