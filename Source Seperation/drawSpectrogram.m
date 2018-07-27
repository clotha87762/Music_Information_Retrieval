function [  ] = drawSpectrogram( pathlist )
%DRAWSPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here
    addpath(genpath('MIRtoolbox1.6.1'))
    for i=1:size(pathlist,1),
       path = pathlist(i,:);
       audio = miraudio(path,'Normal');
       frame = mirframe(tempAudio,'Length',w,'sp','Hop',h,'sp');
    
    end
    rpath(genpath('MIRtoolbox1.6.1'))
end

