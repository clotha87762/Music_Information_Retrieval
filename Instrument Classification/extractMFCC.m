function [feature, raw] = extractMFCC(file,w,h)
% file: file name
% w: window size (in samples)
% h: hop size (in samples)

% normalize with respect to RMS energy 
a = miraudio(file,'Normal');
% window length 1024 samples, hop size 512 samples
% w = 1024;
% h = 512;
f = mirframe(a,'Length',w,'sp','Hop',h,'sp');
% compute the spectrogram
S = mirspectrum(f);
% while S is an object, x is the values inside 
% x = mirgetdata(S); % size: 513 x number of frames
mfcc = mirmfcc(S,'Rank',1:13);
raw = mirgetdata(mfcc);
feature = [mean(raw,2); std(raw,0,2)]; % e.g. feature(21) = std(raw(1,:))
