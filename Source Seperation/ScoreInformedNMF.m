function [ H ] = ScoreInformedNMF( V, W ,initH  )
%SCOREINFORMEDNMF Summary of this function goes here
%   Detailed explanation goes here

    HH = initH;
    H = HH;
    delta = 999999999;
    mindelta = 999999999;
    minH = H;
    i = 0;

    while i<100,
         H = HH;
         HH = HH.*(((W)'*V)./ ((W)'*(W*HH+eps)));
         delta = norm(HH-H); 
         if delta<mindelta,
             mindelta = delta;
             minH = HH;
            % disp(mindelta);
         end
         i = i + 1;
         %disp(delta);
    end
    H = minH;


end

