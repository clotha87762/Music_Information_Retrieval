

addpath(genpath('MIRtoolbox1.6.1'))

%% Q1
clear all;

a = audioread('audio/validation/01_vio.wav');
b = audioread('audio/validation/01_cla.wav');
c = audioread('audio/validation/01_mix.wav');
n = randn(length(a),1);

SDR = bss_eval_sources([c';c']/2,[a';b'])
SDR = bss_eval_sources([a';b'],[a';b'])
SDR = bss_eval_sources([b';b'],[a';b'])
SDR = bss_eval_sources([2*a';2*b'],[a';b'])
SDR = bss_eval_sources((a+0.01*n)',a')
SDR = bss_eval_sources((a+0.1*n)',a')
SDR = bss_eval_sources((a+n)',a')
SDR = bss_eval_sources((a+0.01*b)',a')
SDR = bss_eval_sources((a+0.1*b)',a')
SDR = bss_eval_sources((a+b)',a')


%% Q2
clear all;
clims = [-30 40];
[vio fs] = audioread('audio/train/vio/vio_64.wav');
vio = vio';
w = 2048;
h = 1024;
nfft = w;
%[spect f t ] = stft(vio,w,h,nfft,fs);
[spect f t] = spectrogram(vio,w,h,nfft,fs);
%spect = abs(spect(1:nfft/2+1,:));
%f = f * fs / 57.2958;
lspect = spect;
hspect = spect;
%lspect(f>1200,:) = 0;
hspect(f<1200,:) = 0;
%lvio = istft(lspect,h,nfft,fs);
%hvio = istft(hspect,h,nfft,fs);
lvio = ispectrogram(lspect);
lvio = lvio';
hvio = ispectrogram(hspect);
hvio = hvio';
audiowrite('vio_64_lp.wav',lvio,fs);
audiowrite('vio_64_hp.wav',hvio,fs);


figure();
%[s f t] = spectrogram(vio,w,h,'yaxis');
imagesc(t,f,(abs(spect)));
colorbar;
axis xy;



%{
figure();
imagesc(t,f,(abs(lspect)));
colorbar;
axis xy;

figure();
imagesc(t,f,(abs(hspect)));
colorbar;
axis xy;

%}


figure();
[sss f t ] = spectrogram(lvio,w,h,nfft,fs);
imagesc(t,f,(abs(sss)));
colorbar;
axis xy;

figure();
[sss f t ] = spectrogram(hvio,w,h,nfft,fs);
imagesc(t,f,(abs(sss)));
colorbar;
axis xy;



%{
figure();
spectrogram(vio,w,h,nfft,fs,'yaxis');
figure();
spectrogram(lvio,w,h,nfft,fs,'yaxis');
figure();
spectrogram(hvio,w,h,nfft,fs,'yaxis');
%}
l1 = length(vio);
l2 = length(lvio);
l3 = l1 - l2;
if l3~=0,
vio(l1-l3+1:l1) = [];
end
SDRL = bss_eval_sources(lvio,vio);
SDRH = bss_eval_sources(hvio,vio);

%% Q3
clear all;
win = 2048;
hop = win/2;
nfft = win;
R = 3;

[vio64 fs] = audioread('audio/train/vio/vio_64.wav');
[V f t] = spectrogram(vio64,win,hop,nfft,fs);
%f = f* fs  / 57.2958;
absV = abs(V);
[W H] = NMF(absV,R,0.1);

figure();
imagesc(1:3,f,((W)));
colorbar;
axis xy;

figure();
imagesc(t,f,(absV));
colorbar;
axis xy;

reV = W*H;
%reVio_64 = istft(reV,hop,nfft,fs);
reVio_64 = ispectrogram(reV);
reVio_64 = reVio_64';
audiowrite('reVio_64.wav',reVio_64,fs);
figure();
xxxx = 1:length(reVio_64);
plot(xxxx,reVio_64);

[vio88 fs] = audioread('audio/train/vio/vio_88.wav');
[V f t] = spectrogram(vio88,win,hop,nfft,fs);
%f = f* fs / 57.2958;
absV = abs(V);
[W H] = NMF(absV,R,0.1);

figure();
imagesc(1:3,f,(abs(W)));
colorbar;
axis xy;


[cla64 fs] = audioread('audio/train/cla/cla_64.wav');
[V f t] = spectrogram(cla64,win,hop,nfft,fs);
f = f* fs / 57.2958;
absV = abs(V);
[W H] = NMF(absV,R,0.1);

figure();
imagesc(1:3,f,(abs(W)));
colorbar;
axis xy;

%% Q4
    clear all;
   win = 2048;
   hop = 1024;
   nfft = win;
   vioDic = [];
   R = 3;
   
   alpha = 0.1;
    algs = {@nmf_kl_ns, @nmf_kl_sparse_es, @nmf_kl_sparse_v, @nmf_euc_sparse_es};
   
   
   for j = 13:49,
      name = strcat('audio/train/piano/piano_',num2str(j),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
       [W H] =NMF(absV,R,0.1);%nmf_alg(absV, R, 'alg',algs{2}, 'verb', 2, ...
               %                   'norm_w', 2, 'alpha', alpha, 'niter', 100);
       vioDic = [vioDic W];
   end
   save('Piano_Template.mat','vioDic');
 %{
      figure();
      imagesc(1:135,f,abs(vioDic));
      colorbar;
      axis xy;
   %}   
      %%
      [valid fs] = audioread('audio/validation/c1.wav_voice_rpca.wav');
    
      if size(valid,2)==2,
      valid(:,2) = [];
      end
     % valid = upsample(valid,int32(44100/8000));
      [vSpect ff tt] = spectrogram(valid,win,hop,nfft,fs);
      cangle = angle(vSpect);
    
      %[VHH] = FixedWNMF(abs(vSpect),vioDic);
      
      %nmf_alg(abs(vSpect), size(vioDic,2), 'alg',algs{2}, 'verb', 2, ...
       %                           'norm_w', 2, 'alpha', alpha, 'niter', 500,'W',vioDic);
      %%
      % temp = HH;
      %{
       HH = imgaussfilt(HH, 1);
       for a = 1: size(HH,2),
         % temp = medfilt1(HH(a,:),10);
          %HH(a,:) = temp;
          [~,mi] = max(HH(:,a));
          temp = HH(mi,a);
          HH(:,a) = 0;
          HH(mi,a) = temp;
       end
    %}
      % HH(find(HH<0.1)) = 0;
       %HH(find(HH>0.1)) = 0.3;
       %HH(1:30,:) = 0;
   vioDic2 = [];
   pulsenarrowTemp = [];
   squareTemp = [];
   triangleTemp = [];
   doublesawTemp = [];
   spikyTemp = [];
   for j = 30:92,
      name = strcat('audio/train/8bit/pulsenarrow/pulsenarrow_',num2str(j),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
       [W H] = nmf_alg(absV, R, 'alg', algs{2}, 'verb', 2, ...
                                  'norm_w', 2, 'alpha', alpha, 'niter', 500);
       vioDic2 = [vioDic2 W];
       pulsenarrowTemp = [pulsenarrowTemp W];
   end
   save('PulsenarrowTemp','pulsenarrowTemp');
   
   for j = 30:92,
      name = strcat('audio/train/8bit/spiky/spiky_',num2str(j),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
       [W H] = nmf_alg(absV, R, 'alg', algs{2}, 'verb', 2, ...
                                  'norm_w', 2, 'alpha', alpha, 'niter', 500);
       vioDic2 = [vioDic2 W];
       spikyTemp = [spikyTemp W];
   end
   save('SpikyTemp','spikyTemp');
   
   for j = 30:92,
      name = strcat('audio/train/8bit/square_',num2str(j),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
       [W H] = nmf_alg(absV, R, 'alg', algs{2}, 'verb', 2, ...
                                  'norm_w', 2, 'alpha', alpha, 'niter', 500);
       vioDic2 = [vioDic2 W];
       squareTemp = [squareTemp W];
   end
   save('SquareTemp','squareTemp');
   
    for j = 30:90,
      name = strcat('audio/train/8bit/triangle/triangle_',num2str(j),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
       [W H] = nmf_alg(absV, R, 'alg', algs{2}, 'verb', 2, ...
                                  'norm_w', 2, 'alpha', alpha, 'niter', 500);
       vioDic2 = [vioDic2 W];
       tringleTemp = [triangleTemp W];
    end
    
    save('TriangleTemp',triangleTemp);
    
    
    for j = 30:90,
      name = strcat('audio/train/8bit/doublesaw/doublesaw_',num2str(j),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
       [W H] = nmf_alg(absV, R, 'alg', algs{2}, 'verb', 2, ...
                                  'norm_w', 2, 'alpha', alpha, 'niter', 500);
       vioDic2 = [vioDic2 W];
       doublesawTemp = [doublesawTemp W];
    end
   save('DoublesawTemp','doublesawTemp');
   
   %%
   save('8Bit_Template','vioDic2');
   %%
  %{
   HH = zeros(111,size(VHH,2));
   for j=1:size(HH,1),
      HH(j,:) = VHH(j,:);
   end
   %}
 % HH(find(HH<max(max(HH))/10.0)) = 0;
  
 HH =FixedWKLNMF(abs(vSpect),vioDic2); 
 tempHH = HH;
 %[~, HH] = nmf_alg(abs(vSpect), size(vioDic2,2), 'alg',@nmf_beta, 'verb', 2, ...
  %                              'norm_w', 2, 'beta', 1, 'niter', 500,'W',vioDic2);
  %%
  
  mm = zeros(10,1);
  mi = zeros(5,1);
   for a=1:size(HH,2),
     for j = 1:5,
       [ mm(j)  mi(j)] = max(HH(:,a));
       HH(mi(j),a) = 0;
     end
     HH(:,a) = 0;
     for j=1:5,
        HH(mi(j),a) = mm(j); 
     end
  end
  
  
  for a=1:size(HH,2),
     [mm mi] = max(HH(:,a)); 
     HH(find(HH(:,a)<(mm/5)),a) = 0;
    
  end
  
 
  
  %%    
 
       figure();
       imagesc(tt,1:120,abs(HH));
       colorbar;
       axis xy;
       
      % complexW = vioDic.*cos(angle) + i*vioDic.*sin(angle);
       Reconstructed = vioDic2 * HH ;
       complexR =Reconstructed.*cos(cangle) + i*Reconstructed.*sin(cangle);
       %rvio = istft(complexR,hop,nfft,fs);
       rvio = ispectrogram(complexR);
       rvio = rvio';
       
       audiowrite('reconstructed_18_bit.wav',rvio,fs);
       
       
        l1 = length(valid);
        l2 = length(rvio);
        l3 = l1 - l2;
        if l3~=0,
        valid(l1-l3+1:l1) = [];
        end
      % SDR = bss_eval_sources(rvio,valid');
       
 %%  test Denoise
        rvio2 = audioread('Some_one_like_you(°t¼Ö).wav');
        %qq = wden(rvio2,'heursure','h','one',4,'sym8');
        sigden = cmddenoise(rvio2,'sym4',6,'s',3);
        audiowrite('test2.wav',sigden,fs);
%%  Q5

    [valid fss] = audioread('audio/validation/01_cla.wav');
    [vSpect ff tt] = spectrogram(valid,win,hop,nfft,fs);
    cangle = angle(vSpect);
    HH = FixedWNMF(abs(vSpect),vioDic);
    
       figure();
       imagesc(t,1:135,abs(HH));
       colorbar;
       axis xy;
       
       Reconstructed = vioDic * HH ;
       complexR = Reconstructed.*cos(cangle) + i*Reconstructed.*sin(cangle);
       %rcla = istft(complexR,hop,nfft,fs);
       rcla = ispectrogram(complexR);
       rcla = rcla';
       audiowrite('reconstructed_01_cla.wav',rcla,fs);
       
        l1 = length(valid);
        l2 = length(rcla);
        l3 = l1 - l2;
        if l3~=0,
        valid(l1-l3+1:l1) = [];
        end
       SDR = bss_eval_sources(rcla,valid');

%% Train Dictionary
   clear all;
   win = 2048;
   hop = 1024;
   nfft = win;
   vioDic = [];
   claDic = [];
   R = 3;

   for x = 55:99,
      name = strcat('audio/train/vio/vio_',num2str(x),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
     %[W H] = NMF(absV,R,0.05);
     [W H] = KLNMF(absV,R,200); 
     vioDic = [vioDic W];
   end
 
    for x=50:89,
      
      name = strcat('audio/train/cla/cla_',num2str(x),'.wav');
      [audio fs] = audioread(name);
      [V f t ] = spectrogram(audio,win,hop,nfft,fs);
       absV = abs(V);
      % [W H] = NMF(absV,R,0.05);
        [W H] = KLNMF(absV,R,200); 
       claDic = [claDic W];     
    end

   %% Q6

    isHard = 1;
    c = 1;
    SDR = zeros(2,5);
    SIR = zeros(2,5);
    SAR = zeros(2,5);
    perm = zeros(2,5);
   
    for x=1:5,
      
       name = strcat('audio/validation/0',num2str(x),'_mix.wav'); 
       [audio fs ] = audioread(name);
       [V f t] = spectrogram(audio,win,hop,nfft,fs);
       cangle = angle(V);
       absV = abs(V);
 
       vioH = FixedWNMF(absV,vioDic);
       claH = FixedWNMF(absV,claDic);
      % vioH = FixedWKLNMF(absV,vioDic);
       %claH = FixedWKLNMF(absV,claDic);
       rVio = vioDic * vioH;
       rCla = claDic * claH;
 
        if isHard==1,
           mask = rVio > rCla ;
           A = absV .* mask;
           mask = rCla > rVio;
           B = absV .* mask;
        else
            mask = (rVio.^c) ./ ( rVio.^c + rCla.^c);
            A = absV .* mask;
            mask = (rCla.^c) ./ ( rVio.^c + rCla.^c);
            B = absV .* mask;
        end
       AA = A.*cos(cangle) + 1i*A.*sin(cangle);
       BB = B.*cos(cangle) + 1i*B.*sin(cangle);
    
       %rv = istft(AA,hop,nfft,fs);
       %rc = istft(BB,hop,nfft,fs);
       rv = ispectrogram(AA);
       rv = rv';
       rc = ispectrogram(BB);
       rc = rc';
       name = strcat('0',num2str(x),'_vio_est.wav');
       audiowrite(name,rv,fs);
       name = strcat('0',num2str(x),'_cla_est.wav');
       audiowrite(name,rc,fs);
        
       name = strcat('audio/validation/0',num2str(x),'_vio.wav'); 
        vv = audioread(name);
        name = strcat('audio/validation/0',num2str(x),'_cla.wav'); 
        cc = audioread(name);
        
        l1 = length(vv);
        l2 = length(rv);
        l3 = l1 - l2;
        if l3~=0,
        vv(l1-l3+1:l1) = [];
        cc(l1-l3+1:l1) = [];
        end
        
       [SDR(:,x) , SIR(:,x), SAR(:,x) , perm(:,x)] = bss_eval_sources([rv;rc], [vv';cc']);
       effect = sum(sum(SDR));
    end




%% Q7

    isHard = 1;
    c = 1;
   
    for x=6:10,
       if x==10,
          name = strcat('audio/test/',num2str(x),'_mix.wav');  
       else
       name = strcat('audio/test/0',num2str(x),'_mix.wav'); 
       end
       [audio fs ] = audioread(name);
       [V f t] = spectrogram(audio,win,hop,nfft,fs);
       cangle = angle(V);
       absV = abs(V);
 
       vioH = FixedWNMF(absV,vioDic);
       claH = FixedWNMF(absV,claDic);
       rVio = vioDic * vioH;
       rCla = claDic * claH;
 
        if isHard==1,
           mask = rVio > rCla ;
           A = absV .* mask;
           mask = rCla > rVio;
           B = absV .* mask;
        else
            mask = (rVio.^c) ./ ( rVio.^c + rCla.^c);
            A = absV .* mask;
            mask = (rCla.^c) ./ ( rVio.^c + rCla.^c);
            B = absV .* mask;
        end
       AA = A.*cos(cangle) + 1i*A.*sin(cangle);
       BB = B.*cos(cangle) + 1i*B.*sin(cangle);
    
       %rv = istft(AA,hop,nfft,fs);
       %rc = istft(BB,hop,nfft,fs);
       rv = ispectrogram(AA);
       rv = rv';
       rc = ispectrogram(BB);
       rc = rc';
       if x==10,
       name = strcat('pred/',num2str(x),'_vio_est.wav');
       else
       name = strcat('pred/0',num2str(x),'_vio_est.wav');
       end
        audiowrite(name,rv,fs);
        
       if x==10,
       name = strcat('pred/',num2str(x),'_cla_est.wav');
       else
       name = strcat('pred/0',num2str(x),'_cla_est.wav');
       end

       audiowrite(name,rc,fs);
        
       
    end
%%  Adv 1
    
    isHard = 1;
    c = 1;
    SDR = zeros(2,5);
    SIR = zeros(2,5);
    SAR = zeros(2,5);
    perm = zeros(2,5);
     win = 2048;
   hop = 1024;
   nfft = win;
   R = 3;

    for x = 1:5,
         
       name = strcat('audio/validation/0',num2str(x),'_mix.wav'); 
       [audio fs ] = audioread(name);
       [V f t] = spectrogram(audio,win,hop,nfft,fs);
       cangle = angle(V);
       absV = abs(V);
 
       %%%%%%%%%%%%%
        onset= [];
        offset = [];
        kind = [];
        pitch = [];
        fid = fopen(strcat('score-info/0',num2str(x),'.txt'));
        while ~feof(fid),
             d = fscanf(fid,'%d',1);
             onset = [onset d];
             
             d = fscanf(fid,'%d',1);
             if d>5000,
                d = 5000; 
             end
             offset = [offset d];
             
             d = fscanf(fid,'%d',1);
             pitch = [pitch d];
             d = fscanf(fid,'%d',1);
             kind = [kind d];
        end
        l = size(absV,2);
        initVioH = zeros(R * 45,l);
        initClaH = zeros(R * 40,l);
        for y = 1:length(onset),
             if kind(y)==1,
              for xx = 1:R,
                initVioH(R*(pitch(y)-55)+xx, floor(onset(y)/5000 * l) : floor(offset(y)/5000 * l)) = 1;
              end
             else
               for xx=1:R,
                initClaH(R*(pitch(y)-50)+xx, floor(onset(y)/5000 * l) : floor(offset(y)/5000 * l)) = 1;
               end
           end
        end
       %%%%%%%%%%%%%
       
       
       vioH = ScoreInformedNMF(absV,vioDic,initVioH);
       claH = ScoreInformedNMF(absV,claDic,initClaH);

       rVio = vioDic * vioH;
       rCla = claDic * claH;
 
        if isHard==1,
           mask = rVio > rCla ;
           A = absV .* mask;
           mask = rCla > rVio;
           B = absV .* mask;
        else
            mask = (rVio.^c) ./ ( rVio.^c + rCla.^c);
            A = absV .* mask;
            mask = (rCla.^c) ./ ( rVio.^c + rCla.^c);
            B = absV .* mask;
        end
       AA = A.*cos(cangle) + 1i*A.*sin(cangle);
       BB = B.*cos(cangle) + 1i*B.*sin(cangle);
     
       %rv = istft(AA,hop,nfft,fs);
       %rc = istft(BB,hop,nfft,fs);
       rv = ispectrogram(AA);
       rv = rv';
       rc = ispectrogram(BB);
       rc = rc';
       name = strcat('informed/0',num2str(x),'_vio_est.wav');
       audiowrite(name,rv,fs);
       name = strcat('informed/0',num2str(x),'_cla_est.wav');
       audiowrite(name,rc,fs);
        
       name = strcat('audio/validation/0',num2str(x),'_vio.wav'); 
        vv = audioread(name);
        name = strcat('audio/validation/0',num2str(x),'_cla.wav'); 
        cc = audioread(name);
        
        l1 = length(vv);
        l2 = length(rv);
        l3 = l1 - l2;
        if l3~=0,
        vv(l1-l3+1:l1) = [];
        cc(l1-l3+1:l1) = [];
        end
        
       [SDR(:,x) , SIR(:,x), SAR(:,x) , perm(:,x)] = bss_eval_sources([rv;rc], [vv';cc']);
       effect = sum(sum(SDR));
        
        
            
    end






