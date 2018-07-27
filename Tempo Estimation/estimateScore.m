function [ P , ALOTC ] = estimateScore( t1 , t2 , G ,tempogram , t1i , t2i)
%ESTIMATESCORE Summary of this function goes here
%   Detailed explanation goes here
       %[ ~ , t1i ] =  min(abs(BPM(:) - double(t1)));
       %[ ~ , t2i ] = min(abs(BPM(:) - double(t2)));
       %t1i = int32(t1i);  t2i = int32(t2i);
         f1 = 0; f2 = 0;
         ALOTC = 0;
         P = 0;
            for i=1:size(tempogram,2),
               f1 = f1 + abs(tempogram(t1i,i));
               f2 = f2 + abs(tempogram(t2i,i));
            end
               f1 = f1 / size(tempogram,2);
               f2 = f2 / size(tempogram,2);
               
               %f1 = f1*f1;
               %f2 = f2*f2;
               
               s1 = f1 / (f1+f2);
               s2 = 1-s1;
               %disp(num2str(s1));
               %disp(num2str(s2));
               
               T1 = 0;
               T2 = 0;
               if abs(double(G-t1)/double(G)) <=0.08,
                  T1 = 1 ;
               end
               if abs(double(G-t2)/double(G)) <= 0.08,
                  T2 = 1; 
               end
               P = (s1 * T1) + (s2 * T2);
                if abs(double(G-t1)/G) <0.08 || abs(double(G-t2)/G) < 0.08,
                  ALOTC =  1; 
               end

end

