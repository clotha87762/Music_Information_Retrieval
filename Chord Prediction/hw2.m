%%

clear ;clc;close all;
%%

addpath(genpath('MATLAB-Chroma-Toolbox_2.0'));

%% Produce Pitch
genre = 'blues';

n = 99;
while n < 100,

    nstr = num2str(n);
    for i=1:(5-length(nstr)),
       nstr = strcat('0',nstr); 
    end
    
   filename = strcat('genres/',genre,'/',genre,'.',nstr,'.au');
  
   filename = strcat('piano_01.wav');
   %filename = 'Systematic_Chord-C-Major_Eight-Instruments.wav'; 
   
    [f_audio , fs ] = audioread(filename);
    clear parameter
    shiftFB = estimateTuning(f_audio);

     clear parameter
     parameter.winLenSTMSP = 4410;
     parameter.fs = 22050;
     parameter.save = 1;
     parameter.saveDir = 'mpitch/';
     parameter.saveFilename = strcat(genre,'_',nstr);
     parameter.shiftFB = shiftFB;
     parameter.saveAsTuned = 1;
     [f_pitch,sideinfo] = audio_to_pitch_via_FB(f_audio,parameter);
     
    parameter.usePitchNameLabels = 1;
    parameter.title = 'Logarithmic compression of amplitude';
    parameter.featureRate = sideinfo.pitch.featureRate;
    parameter.xlabel = 'Time [Seconds]';
    parameter.ylabel = 'Pitch';
    visualizePitch(log(5*f_pitch+1),parameter);
     %{
     parameter.usePitchNameLabels = 1;
     parameter.title = 'Logarithmic compression of amplitude';
     parameter.featureRate = sideinfo.pitch.featureRate;
     parameter.xlabel = 'Time [Seconds]';
     parameter.ylabel = 'Pitch';
     visualizePitch(log(5*f_pitch+1),parameter);
    %}
     n = n + 1;
     
     if n==100 && strcmp(genre,'blues'),
        genre = 'country';
        n = 0;
     elseif n==100 && strcmp(genre,'country'),
        genre = 'disco';
         n = 0;
     elseif n==100 && strcmp(genre,'disco'),
         genre = 'hiphop';
         n = 0;
     elseif n==100 && strcmp(genre,'hiphop'),
        genre = 'jazz';
         n = 0;
     elseif n==100 && strcmp(genre,'jazz'),
        genre = 'metal';
         n = 0;
     elseif n==100 && strcmp(genre,'metal'),
        genre = 'pop';
        n = 0;
     elseif n==100 && strcmp(genre,'pop'),
         genre = 'reggae';
          n = 0;
     elseif n==100&&strcmp(genre,'reggae'),
        genre = 'rock';
        n = 0;
     end
     
     
end

%% Produce Chroma

clear;
close all hidden;
genre = 'blues';
win_len = 4410;
directory = 'mpitch/';
n = 0;
gamma = 100;
while n<100,
  
    clear parameter;
    
    nstr = num2str(n);
    for i=1:(5-length(nstr)),
       nstr = strcat('0',nstr); 
    end
    
   pitch_name = strcat(genre,'_',nstr,'_pitch_',num2str(win_len));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loads pitch features (f_pitch) and computes chroma features (f_chroma)
    %
    % Note: feature filename is specified by WAV filename
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load(strcat(directory,pitch_name)); % load f_pitch and sideinfo;

    parameter.vis = 0;
    %parameter.save = 1;
    %parameter.save_dir = 'data_feature/';
    %parameter.save_filename = strcat(sideinfo.wav.filename(1:length(sideinfo.wav.filename)-4));
    [f_chroma_norm,sideinfo] = pitch_to_chroma(f_pitch,parameter,sideinfo);

    parameter.applyLogCompr = 1;
    
    parameter.factorLogCompr = gamma;
    
    f_logchroma_norm = pitch_to_chroma(f_pitch,parameter,sideinfo);
    
    save(strcat('mchroma/',genre,'_',nstr,'_logchroma_',num2str(parameter.factorLogCompr)),'f_logchroma_norm','parameter');
     
    n = n+ 1;
    
    if n==100 && strcmp(genre,'blues'),
        genre = 'country';
        n = 0;
     elseif n==100 && strcmp(genre,'country'),
        genre = 'disco';
         n = 0;
     elseif n==100 && strcmp(genre,'disco'),
         genre = 'hiphop';
         n = 0;
    elseif n==100 && strcmp(genre,'hiphop'),
        genre = 'jazz';
         n = 0;
     elseif n==100 && strcmp(genre,'jazz'),
        genre = 'metal';
         n = 0;
     elseif n==100 && strcmp(genre,'metal'),
        genre = 'pop';
        n = 0;
     elseif n==100 && strcmp(genre,'pop'),
         genre = 'reggae';
          n = 0;
     elseif n==100&&strcmp(genre,'reggae'),
        genre = 'rock';
        n = 0;
     end
end

%%  Algorithm 1 without weight

clear;

MIXED = 1;

if MIXED == 1,
    mixed = '_mixed';
else
    mixed = '';
end

genre = 'blues';
logScale = 100;


tempMajor = [1 0 1 0 1 1 0 1 0 1 0 1 ];
tempMinor = [1 0 1 1 0 1 0 1 1 0 1 0 ];
tempMajor = circshift(tempMajor,[0 5]);
tempMinor = circshift(tempMinor,[0 5]);
BTemp = zeros(24,12);


for i=1:12,
    BTemp(i,:) = tempMajor;
    tempMajor = circshift(tempMajor,[0,1]);
end
for i=13:24,
   BTemp(i,:) = tempMinor;
   tempMinor = circshift(tempMinor,[0,1]);
end


totalNum = 0;
totalMatched = 0;
genreMatched = 0;
finalKey = zeros(100,1);
flag = 1;
n = 0;
genreNum = 0;
while n<100,
    
  
    if flag==1,
        correctAns = ones(100,1) .* -1.0;
        for qq = 0:99,
             nstr = num2str(qq);
            for i=1:(5-length(nstr)),
               nstr = strcat('0',nstr); 
            end
             a = (textread(strcat('labels/',genre,'/',genre,'.',nstr,'.lerch.txt')));
             correctAns(qq+1) = a(1);
             if correctAns(qq+1)<0||correctAns(qq+1)>23,
                %disp(qq); 
             end
         end
        finalKey = zeros(100,1);
        flag = 0;
    end
    
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   
    
    totalNum = totalNum +1;
    genreNum = genreNum + 1;
    
    nstr = num2str(n);
    for i=1:(5-length(nstr)),
       nstr = strcat('0',nstr); 
    end
    
    load(strcat('mchroma/',genre,'_',nstr,'_logchroma_',num2str(logScale)));
    
    sumVec = sum(f_logchroma_norm,2);
    sumVec = circshift(sumVec,[3 0]);
    [m,index] = max(sumVec);
    xAvg = mean(sumVec);
    up = 0;down = 0;downX = 0;downY = 0;
    for i=1:12,
        up = up + (sumVec(i) - xAvg) * (BTemp(index,i) - (7/12));
        downX = downX +  (sumVec(i) - xAvg)^2;
        downY = downY +  (BTemp(index,i) - (7/12))^2;
    end
    down = (downX * downY)^(0.5);
    RMajor = up/down;
    
    up = 0;down = 0;downX = 0;downY = 0;
    for i=1:12,
        up = up + (sumVec(i) - xAvg) * (BTemp(index + 12,i) - (7/12));
        downX = downX +  (sumVec(i) - xAvg)^2;
        downY = downY +  (BTemp(index + 12,i) - (7/12))^2;
    end
    down = (downX * downY)^(0.5);
    RMinor = up/down;
    if RMajor>=RMinor,
        key = index - 1;
    elseif RMajor< RMinor,
        key = index + 11; 
    end
    
    if correctAns(n+1)~=-1 && key==correctAns(n+1),
        genreMatched = genreMatched + 1;
        totalMatched = totalMatched + 1;
    end
    %@@@@@@@@@@@@
    if correctAns(n+1)~=-1&&MIXED==1,
       if key <= 11 && key>=0,   
           if key==correctAns(n+1)+12,
              genreMatched = genreMatched + 0.2;
              totalMatched = totalMatched + 0.2;
           elseif key==mod(correctAns(n+1)+7,12),
               genreMatched = genreMatched + 0.5;
               totalMatched = totalMatched + 0.5;
           elseif key==(mod(correctAns(n+1)+9,12)+12),
              genreMatched = genreMatched + 0.3;
              totalMatched = totalMatched + 0.3;
           end
        elseif key<=23 && key>=12,
           if key == correctAns(n+1)-12,
              genreMatched = genreMatched + 0.2;
              totalMatched = totalMatched + 0.2;
           elseif key == mod(correctAns(n+1)-9,12),
              genreMatched = genreMatched + 0.3;
              totalMatched = totalMatched + 0.3;
           elseif key== mod(correctAns(n+1)-5,12) + 12,
               genreMatched = genreMatched + 0.5;
              totalMatched = totalMatched + 0.5;
           end
           
       end
       %@@@@@@@@@@@@@@@@@@@
       
    end
    
    finalKey(n+1) = key;
    if correctAns(n+1) == -1,
       totalNum = totalNum - 1;
       genreNum = genreNum - 1;
    end
    
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    
    n = n+1;
    if n==100 && strcmp(genre,'blues'),
        genreAccuracy = genreMatched / genreNum;
        save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
        fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'country';      
        n = 0;
        genreNum = 0;
        genreMatched = 0;
        flag = 1;
     elseif n==100 && strcmp(genre,'country'),
        genreAccuracy = genreMatched /genreNum;
        save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
         fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'disco';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
           flag = 1;
     elseif n==100 && strcmp(genre,'disco'),
         genreAccuracy = genreMatched / genreNum;
         save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
          fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
         genre = 'hiphop';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
           flag = 1;
    elseif n==100 && strcmp(genre,'hiphop'),
         genreAccuracy = genreMatched / genreNum;
         save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
          fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
         genre = 'jazz';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
           flag = 1;
     elseif n==100 && strcmp(genre,'jazz'),
         genreAccuracy = genreMatched / genreNum;
         save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
         fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy); 
        genre = 'metal';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
           flag = 1;
     elseif n==100 && strcmp(genre,'metal'),
        genreAccuracy = genreMatched / genreNum;
       save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
         fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'pop';
        n = 0;
        genreNum = 0;
        genreMatched = 0;
          flag = 1;
     elseif n==100 && strcmp(genre,'pop'),
         genreAccuracy = genreMatched / genreNum;
         save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
          fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
         genre = 'reggae';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
          flag = 1;
     elseif n==100&&strcmp(genre,'reggae'),
          genreAccuracy = genreMatched / genreNum;
         save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
          fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'rock';
        n = 0;
        genreNum = 0;
        genreMatched = 0;
        flag = 1;
     end
    
end
  genreAccuracy = genreMatched /genreNum;
  save(strcat('mKey/Bin_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
   fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);

totalAccuracy = totalMatched / totalNum;
save(strcat('mkey/Bin_Total_',genre,'_gamma',num2str(logScale),mixed),'totalAccuracy');
 fprintf('Total_gamma:%d_accuracy:%d\n\n\n',logScale,totalAccuracy);
  

%% Algorithm 2

clear;

MIXED = 1;

if MIXED == 1,
    mixed = 'mixed';
else 
    mixed = '';
end

genre = 'blues';
logScale = 100;


tempMajor = [6.35 2.23 3.48 2.33 4.38 4.09 2.52 5.19 2.39 3.66 2.29 2.88 ];
tempMinor = [6.33 2.68 3.52 5.38 2.60 3.53 2.54 4.75 3.98 2.69 3.34 3.17  ];
tempMajor = circshift(tempMajor,[0 0]);
tempMinor = circshift(tempMinor,[0 0]);
KSTemp = zeros(24,12);
tMajorAvg = mean(tempMajor);
tMinorAvg = mean(tempMinor);


for i=1:12,
    KSTemp(i,:) = tempMajor;
    tempMajor = circshift(tempMajor,[0,1]);
end
for i=13:24,
   KSTemp(i,:) = tempMinor;
   tempMinor = circshift(tempMinor,[0,1]);
end


totalNum = 0;
totalMatched = 0;
genreMatched = 0;
finalKey = zeros(100,1);
flag = 1;
n = 0;
genreNum = 0;
while n<100,
    
  
    if flag==1,
         correctAns = ones(100,1) .* -1.0;
        for qq = 0:99,
             nstr = num2str(qq);
            for i=1:(5-length(nstr)),
               nstr = strcat('0',nstr); 
            end
             a = (textread(strcat('labels/',genre,'/',genre,'.',nstr,'.lerch.txt')));
             correctAns(qq+1) = a(1);
             if correctAns(qq+1)<0||correctAns(qq+1)>23,
                %disp(n); 
             end
        end
        finalKey = zeros(100,1);
        flag = 0;
    end
    
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    
    
    totalNum = totalNum +1;
    genreNum = genreNum + 1;
    
    nstr = num2str(n);
    for i=1:(5-length(nstr)),
       nstr = strcat('0',nstr); 
    end
    
    load(strcat('mchroma/',genre,'_',nstr,'_logchroma_',num2str(logScale)));
    
    sumVec = sum(f_logchroma_norm,2);
    sumVec = circshift(sumVec,[3 0]);
    xAvg = mean(sumVec);
    RMax = -100;
    key = -1;
    for i=1:24,
         up = 0;down = 0;downX = 0;downY = 0;
             if i<=12&&i>=0,
                 t = tMajorAvg;
             else
                 t = tMinorAvg; 
             end
        for j=1:12,
             up = up + (sumVec(j) - xAvg) * (KSTemp(i,j) - (t));
             downX = downX +  (sumVec(j) - xAvg)^2;
             downY = downY +  (KSTemp(i,j) - (t))^2;
        end
         down = (downX * downY)^(0.5);
         R = up/down;
         if R>RMax,
            RMax = R;
            key = i-1;
         end
    end
    
   % fprintf('%d  %d \n',key,correctAns(n+1));
    
    if correctAns(n+1)~=-1 && key==correctAns(n+1),
        genreMatched = genreMatched + 1;
        totalMatched = totalMatched + 1;
    end
    %@@@@@@@@@@@@
    if correctAns(n+1)~=-1&&MIXED==1,
       if key <= 11 && key>=0,   
           if key==correctAns(n+1)+12,
              genreMatched = genreMatched + 0.2;
              totalMatched = totalMatched + 0.2;
           
           elseif key==mod(correctAns(n+1)+7,12),
               genreMatched = genreMatched + 0.5;
               totalMatched = totalMatched + 0.5;
           
           elseif key==(mod(correctAns(n+1)+9,12)+12),
              genreMatched = genreMatched + 0.3;
              totalMatched = totalMatched + 0.3;
           end
        elseif key<=23 && key>=12,
           if key == correctAns(n+1)-12,
              genreMatched = genreMatched + 0.2;
              totalMatched = totalMatched + 0.2;
           
           elseif key == mod(correctAns(n+1)-9,12),
              genreMatched = genreMatched + 0.3;
              totalMatched = totalMatched + 0.3;
           
           elseif key== mod(correctAns(n+1)-5,12) + 12,
               genreMatched = genreMatched + 0.5;
              totalMatched = totalMatched + 0.5;
           end
           
       end
       %@@@@@@@@@@@@@@@@@@@
       
    end
    
    finalKey(n+1) = key;
    if correctAns(n+1) == -1,
       totalNum = totalNum - 1;
       genreNum = genreNum - 1;
    end
    
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    
    n = n+1;
    if n==100 && strcmp(genre,'blues'),
        genreAccuracy = genreMatched / genreNum;
      save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
        fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'country';      
        n = 0;
        genreNum = 0;
        genreMatched = 0;
        flag = 1;
     elseif n==100 && strcmp(genre,'country'),
        genreAccuracy = genreMatched /genreNum;
       save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
         fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'disco';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
         flag = 1;
     elseif n==100 && strcmp(genre,'disco'),
         genreAccuracy = genreMatched / genreNum;
       save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
          fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
         genre = 'hiphop';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
         flag = 1;
    elseif n==100 && strcmp(genre,'hiphop'),
         genreAccuracy = genreMatched / genreNum;
         save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
          fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
         genre = 'jazz';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
         flag = 1;
     elseif n==100 && strcmp(genre,'jazz'),
         genreAccuracy = genreMatched / genreNum;
       save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
         fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy); 
        genre = 'metal';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
         flag = 1;
     elseif n==100 && strcmp(genre,'metal'),
        genreAccuracy = genreMatched / genreNum;
        save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
         fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'pop';
        n = 0;
        genreNum = 0;
        genreMatched = 0;
        flag = 1;
     elseif n==100 && strcmp(genre,'pop'),
         genreAccuracy = genreMatched /genreNum;
         save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
         fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
         genre = 'reggae';
         n = 0;
         genreNum = 0;
         genreMatched = 0;
         flag = 1;
     elseif n==100&&strcmp(genre,'reggae'),
          genreAccuracy = genreMatched / genreNum;
          save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
          fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);
        genre = 'rock';
        n = 0;
        genreNum = 0;
        genreMatched = 0;
        flag = 1;
     end
    
end
  genreAccuracy = genreMatched / genreNum;
  save(strcat('mKey/KS_',genre,'_',num2str(logScale),mixed),'genreAccuracy','finalKey');
   fprintf('%s_gamma:%d_accuracy:%d\n',genre,logScale,genreAccuracy);

totalAccuracy = totalMatched / totalNum;
save(strcat('mkey/KS_Total_',genre,'_gamma',num2str(logScale),mixed),'totalAccuracy');
  fprintf('Total_gamma:%d_accuracy:%d\n\n',logScale,totalAccuracy);
 
 
 
 
 