function [ P , ALOTC ] = estimateScore2(  t1 , t2 , G ,tempogram , t1i , t2i ,genre )
%ESTIMATESCORE2 Summary of this function goes here
%   Detailed explanation goes here
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
               
               %s1 = f1 / (f1+f2);
               %s2 = 1-s1;
               %disp(num2str(s1));
               %disp(num2str(s2));
               if strcmp(genre,'ChaChaCha'),
                  s1 = 1;
                  s2 = 0;
               elseif strcmp(genre,'Jive'),
                  s1 = 0;
                  s2 = 1;
               elseif strcmp(genre,'Quickstep'),
                  s1 = 0;
                  s2 = 1;
               elseif strcmp(genre,'Rumba'),
                  s1 = 1;
                  s2 = 0;
               elseif strcmp(genre,'Samba'),
                  s1 = 1;
                  s2 = 0;
               elseif strcmp(genre,'Tango'),
                 s1 =  0;
                  s2 = 1;
               elseif strcmp(genre,'VienneseWaltz'),
                      s1 = 0;
                      s2 = 1;
               elseif strcmp(genre,'Waltz'),
                      s1 = 1;
                      s2 = 0;
               end
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

