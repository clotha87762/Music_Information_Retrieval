function [BitMusic ] = Get_8_bit_instrument(valid , fs , argin )
%GET_8_BIT_INSTRUMENT Summary of this function goes here
%   Detailed explanation goes here
    
       win = 2048;
       hop = 1024;
       nfft = win;
  


       if size(valid,2)==2,
       valid(:,2) = [];
       end
       % valid = upsample(valid,int32(44100/8000));
       [vSpect ff tt] = spectrogram(valid,win,hop,nfft,fs);
       cangle = angle(vSpect);
       
       dictionary = [];
       for i=1:length(argin),
          if strcmp(argin{i},'triangle'),
             load('TriangleTemp',T);
             dictionary = [ dictionary T];
          elseif strcmp(argin{i},'square'),
             load('SquareTemp',T);
             dictionary = [ dictionary T];
          elseif strcmp(argin{i},'pulsenarrow'),
              load('PulsenarrowTemp',T);
             dictionary = [ dictionary T];
          elseif strcmp(argin{i},'doublesaw'),
              load('DoublesawTemp',T);
             dictionary = [ dictionary T];
          elseif srecmp(argin{i},'spiky'),
             load('SpikyTemp',T);
             dictionary = [ dictionary T];
          end
       end
       
       HH =FixedWKLNMF(abs(vSpect),vioDic2); 
      
      for a=1:size(HH,2),
         [mm mi] = max(HH(:,a)); 
         HH(find(HH(:,a)<(mm/5)),a) = 0;
      end
       
       figure();
       imagesc(t,1:120,abs(HH));
       colorbar;
       axis xy;
       
      % complexW = vioDic.*cos(angle) + i*vioDic.*sin(angle);
       Reconstructed = vioDic2 * HH ;
       complexR =Reconstructed.*cos(cangle) + i*Reconstructed.*sin(cangle);
       %rvio = istft(complexR,hop,nfft,fs);
       BitMusic = ispectrogram(complexR);
       BitMusic  =BitMusic';
end

