%% 
    clear all;
    addpath(genpath('MATLAB-Tempogram-Toolbox_1.0'));
    
    %% Q1
    clear all;
    genreList = cell(8,1);
    %eraseList = [];
    %{
    for i=1:length(genreList),
       if strcmp(genreList(i).name , '.') ||strcmp( genreList(i).name,'..') || strcmp(genreList(i).name,'allBallroomFiles')||strcmp(genreList(i).name , 'nada'),
           eraseList = [eraseList; i];
       end
    end
    genreList(eraseList) = [];
    %}
    genreList{1} = 'ChaChaCha';
    genreList{2} = 'Jive' ;
    genreList{3} = 'Quickstep' ;
    genreList{4} = 'Rumba' ;
    genreList{5} = 'Samba' ;
    genreList{6} = 'Tango' ;
    genreList{7} = 'VienneseWaltz' ;
    genreList{8} = 'Waltz' ;
    
    pause(5);
    disp('Fourier Tempogram');
    fileID = fopen('fourier tempogram.txt','w');
    
    
    for k=1:length(genreList),
        
      
       filenames = dir(strcat('BallroomData/',genreList{k},'/'));
       filenames(1:2) = [];
       dirWav = strcat('BallroomData/',genreList{k},'/');
        
       genreT21 = 0;
        genreT1Q = 0;
        genreT2Q = 0;
        count = 0;
       
       
       
       genreP = zeros(4,1);
       genreALOTC = zeros(4,1) ;
       for j=1:length(filenames),
           
          [audio, fs ] = audioread(strcat(dirWav,filenames(j).name));
          Fs = fs;
           
           parameterNovelty = [];
           [noveltyCurve,featureRate] = audio_to_noveltyCurve(audio, Fs, parameterNovelty);
           
            parameterTempogram = [];
            parameterTempogram.featureRate = featureRate/2;
            parameterTempogram.tempoWindow = 6;         % window length in sec
            parameterTempogram.BPM = 30:1:600;          % tempo values

            [tempogram, T, BPM] = noveltyCurve_to_tempogram_via_DFT(noveltyCurve,parameterTempogram);
            T = T./2;
            %figure();
           % visualize_tempogram(tempogram,T,BPM);
            
            %compSum = sum(tempogram,2);
            
            
            realTempogram = abs(tempogram);
            realSum = sum(realTempogram,2);
            
            %realSum(find(BPM<50))=0;
            %realSum(find(BPM>300))=0;
            
            [MAX Index1] = max(realSum);
            t1 = int32(BPM(Index1));
            range = double(t1 * 0.25);
            downLimit = double(t1 - range);
            upLimit = double(t1 + range);
            [ ~ ,li] = min(abs(BPM(:) - (downLimit)));
            [ ~ ,ui] = min(abs(BPM(:) - (upLimit)));
            li = int32(li); ui = int32(ui);
            realSum(li:ui) = 0;
           % realSum(MAXI) = 0;
            [SUB Index2] = max(realSum);
            t2 = int32(BPM(Index2));
            
           
            if t1> t2,
               tt = t1;
               t1 = t2;
               t2 = tt;
               tt = Index1;
               Index1 = Index2;
               Index2 = tt;
            end
            %}
            
            [~,FName ~]  = fileparts(filenames(j).name);
            G = textread(strcat('BallroomAnnotations/ballroomGroundTruth/',FName,'.bpm'));
            genreT21 = genreT21 +  (double(t2) / double(t1));
            genreT1Q = genreT1Q + (double(t1) / double(G));
            genreT2Q = genreT2Q + (double(t2) / double(G));
           % disp(strcat('t1 = ',num2str(t1),' t2= ',num2str(t2),' G= ',num2str(G)));
             %{
             if ((abs((double(t1)/3) - G)/G) <0.08) || ((abs((double(t2)/3) - G)/G) < 0.08),
                fprintf('%d   %d   %d   %f   %f  \n',t1/3,t2/3,G,abs((t1/3) - G)/G), ((abs((t2/3) - G)/G)  ); 
                count = count +1; 
             end
             %}
             [ P ALOTC] = estimateScore(t1,t2,G,tempogram, Index1 , Index2);
             [ P2 ALOTC2] = estimateScore(t1/2,t2/2,G,tempogram,Index1 , Index2);
             [ P3 ALOTC3] = estimateScore(t1/3,t2/3,G,tempogram,Index1 , Index2);
             [ P4 ALOTC4] = estimateScore(t1/4,t2/4,G,tempogram,Index1 , Index2);

               genreP(1) = genreP(1) + P;
               genreP(2) = genreP(2) + P2;
               genreP(3) = genreP(3) + P3;
               genreP(4) = genreP(4) + P4;
               genreALOTC(1) = genreALOTC(1) + ALOTC ;
               genreALOTC(2) = genreALOTC(2) + ALOTC2 ;
               genreALOTC(3) = genreALOTC(3) + ALOTC3 ;
               genreALOTC(4) = genreALOTC(4) + ALOTC4 ;
             
               
       end
       
       genreP(:) = genreP(:) ./ length(filenames);
       genreALOTC(:) = genreALOTC(:) ./ length(filenames);
       genreT21 = genreT21 / length(filenames);
       genreT1Q = genreT1Q / length(filenames);
       genreT2Q = genreT2Q / length(filenames);
       
        fprintf(fileID,strcat(genreList{k},' P-Score: ',num2str(genreP(1)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/2): ',num2str(genreP(2)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/3): ',num2str(genreP(3)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/4): ',num2str(genreP(4)),'\n'));
        fprintf(fileID,strcat(genreList{k},' ALOTC-Score: ',num2str(genreALOTC(1)),'\n'));
        fprintf(fileID,strcat(genreList{k},' T2/T1 : ',num2str(genreT21),'\n'));
        fprintf(fileID,strcat(genreList{k},' T1/Q : ',num2str(genreT1Q),'\n'));
        fprintf(fileID,strcat(genreList{k},' T2/Q : ',num2str(genreT2Q),'\n'));
      
        fprintf(fileID,'----------------------------- \n');
    end
    
    
       fprintf('\n\n');
    
 %% Q4
   
  clear all;
    genreList = cell(8,1);
    %eraseList = [];
    %{
    for i=1:length(genreList),
       if strcmp(genreList(i).name , '.') ||strcmp( genreList(i).name,'..') || strcmp(genreList(i).name,'allBallroomFiles')||strcmp(genreList(i).name , 'nada'),
           eraseList = [eraseList; i];
       end
    end
    genreList(eraseList) = [];
    %}
    genreList{1} = 'ChaChaCha';
    genreList{2} = 'Jive' ;
    genreList{3} = 'Quickstep' ;
    genreList{4} = 'Rumba' ;
    genreList{5} = 'Samba' ;
    genreList{6} = 'Tango' ;
    genreList{7} = 'VienneseWaltz' ;
    genreList{8} = 'Waltz' ;
    
    pause(5);
    disp('Auto Correlation Function');
    fileID = fopen('auto correlation function.txt','w');
    for k=1:length(genreList),
        
      
       filenames = dir(strcat('BallroomData/',genreList{k},'/'));
       filenames(1:2) = [];
       dirWav = strcat('BallroomData/',genreList{k},'/');
        
       genreT21 = 0;
        genreT1Q = 0;
        genreT2Q = 0;
        count = 0;
       
      
       
       genreP = zeros(7,1);
       genreALOTC = zeros(7,1) ;
       for j=1:length(filenames),
           
          [audio, fs ] = audioread(strcat(dirWav,filenames(j).name));
          Fs = fs;
           
           parameterNovelty = [];
           [noveltyCurve,featureRate] = audio_to_noveltyCurve(audio, Fs, parameterNovelty);
           
           parameterTempogram = [];
            parameterTempogram.featureRate = featureRate;
            parameterTempogram.tempoWindow = 6;                  % window length in sec
            parameterTempogram.maxLag = 2;                    % corresponding to 30 bpm
            parameterTempogram.minLag = 0.1;                  % corresponding to 600 bpm

            [tempogram, T, Lag] = noveltyCurve_to_tempogram_via_ACF(noveltyCurve,parameterTempogram);
            %T = T./2;
            BPM_in = 60./Lag;
            BPM_out = 30:1:600;
            [tempogram,BPM] = rescaleTempoAxis(tempogram,BPM_in,BPM_out);
            %figure();
           % visualize_tempogram(tempogram,T,BPM);
            
            %compSum = sum(tempogram,2);
            realTempogram = abs(tempogram);
            realSum = sum(realTempogram,2);
            [MAX Index1] = max(realSum);
            t1 = int32(BPM(Index1));
            range = double(t1 * 0.25);
            downLimit = double(t1 - range);

            upLimit = double(t1 + range);
            [ ~ ,li] = min(abs(BPM(:) - (downLimit)));
            [ ~ ,ui] = min(abs(BPM(:) - (upLimit)));
            li = int32(li); ui = int32(ui);
            realSum(li:ui) = 0;
           % realSum(MAXI) = 0;
            
            [SUB Index2] = max(realSum);
            t2 = int32(BPM(Index2));
            
            if t1> t2,
               tt = t1;
               t1 = t2;
               t2 = tt;
               tt = Index1;
               Index1 = Index2;
               Index2 = tt;
            end
            
            [~,FName ~]  = fileparts(filenames(j).name);
            G = textread(strcat('BallroomAnnotations/ballroomGroundTruth/',FName,'.bpm'));
            genreT21 = genreT21 +  (double(t2) / double(t1));
            genreT1Q = genreT1Q + (double(t1) / double(G));
            genreT2Q = genreT2Q + (double(t2) / double(G));
           % disp(strcat('t1 = ',num2str(t1),' t2= ',num2str(t2),' G= ',num2str(G)));
             %{
             if ((abs((double(t1)/3) - G)/G) <0.08) || ((abs((double(t2)/3) - G)/G) < 0.08),
                fprintf('%d   %d   %d   %f   %f  \n',t1/3,t2/3,G,abs((t1/3) - G)/G), ((abs((t2/3) - G)/G)  ); 
                count = count +1; 
             end
             %}
             [ P ALOTC] = estimateScore(t1,t2,G,tempogram, Index1 , Index2);
             [ P2 ALOTC2] = estimateScore(t1/2,t2/2,G,tempogram,Index1 , Index2);
             [ P3 ALOTC3] = estimateScore(t1/3,t2/3,G,tempogram,Index1 , Index2);
             [ P4 ALOTC4] = estimateScore(t1/4,t2/4,G,tempogram,Index1 , Index2);
             [ P5 ALOTC5] = estimateScore(t1*2,t2*2,G,tempogram,Index1 , Index2);
             [ P6 ALOTC6] = estimateScore(t1*3,t2*3,G,tempogram,Index1 , Index2);
             [ P7 ALOTC7] = estimateScore(t1*4,4*t2,G,tempogram,Index1 , Index2);

               genreP(1) = genreP(1) + P;
               genreP(2) = genreP(2) + P2;
               genreP(3) = genreP(3) + P3;
               genreP(4) = genreP(4) + P4;
               genreP(5) = genreP(5) + P5;
               genreP(6) = genreP(6) + P6;
               genreP(7) = genreP(7) + P7;
               genreALOTC(1) = genreALOTC(1) + ALOTC ;
               genreALOTC(2) = genreALOTC(2) + ALOTC2 ;
               genreALOTC(3) = genreALOTC(3) + ALOTC3 ;
               genreALOTC(4) = genreALOTC(4) + ALOTC4 ;
               genreALOTC(5) = genreALOTC(5) + ALOTC5 ;
               genreALOTC(6) = genreALOTC(6) + ALOTC6 ;
               genreALOTC(7) = genreALOTC(7) + ALOTC7 ;
             
               
       end
       
       genreP(:) = genreP(:) ./ length(filenames);
       genreALOTC(:) = genreALOTC(:) ./ length(filenames);
       genreT21 = genreT21 / length(filenames);
       genreT1Q = genreT1Q / length(filenames);
       genreT2Q = genreT2Q / length(filenames);
       
        fprintf(fileID,strcat(genreList{k},' P-Score: ',num2str(genreP(1)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/2): ',num2str(genreP(2)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/3): ',num2str(genreP(3)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/4): ',num2str(genreP(4)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(2*T): ',num2str(genreP(5)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(3*T): ',num2str(genreP(6)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(4*T): ',num2str(genreP(7)),'\n'));
        fprintf(fileID,strcat(genreList{k},' ALOTC-Score: ',num2str(genreALOTC(1)),'\n'));
        fprintf(fileID,strcat(genreList{k},' T2/T1 : ',num2str(genreT21),'\n'));
        fprintf(fileID,strcat(genreList{k},' T1/Q : ',num2str(genreT1Q),'\n'));
        fprintf(fileID,strcat(genreList{k},' T2/Q : ',num2str(genreT2Q),'\n'));
      
        fprintf(fileID,'----------------------------- \n');
    end
    
       fprintf('\n\n');
       
 %%    Q6
 
         clear all;
    genreList = cell(8,1);
  
    genreList{1} = 'ChaChaCha';
    genreList{2} = 'Jive' ;
    genreList{3} = 'Quickstep' ;
    genreList{4} = 'Rumba' ;
    genreList{5} = 'Samba' ;
    genreList{6} = 'Tango' ;
    genreList{7} = 'VienneseWaltz' ;
    genreList{8} = 'Waltz' ;
   
    genreList{1} = 'Tango';
    genreList{2} = 'Waltz';
    pause(5);
   fileID = fopen('mixed3.txt','w');
    for k=1:length(genreList),
        
      
       filenames = dir(strcat('BallroomData/',genreList{k},'/'));
       filenames(1:2) = [];
       dirWav = strcat('BallroomData/',genreList{k},'/');
        
       genreT21 = 0;
        genreT1Q = 0;
        genreT2Q = 0;
        count = 0;
       
       
       
       genreP = zeros(7,1);
       genreALOTC = zeros(7,1) ;
       for j=1:length(filenames),
           
          [audio, fs ] = audioread(strcat(dirWav,filenames(j).name));
          Fs = fs;
           
           parameterNovelty = [];
           [noveltyCurve,featureRate] = audio_to_noveltyCurve(audio, Fs, parameterNovelty);
           
           parameterTempogram = [];
            parameterTempogram.featureRate = featureRate/2;
            parameterTempogram.tempoWindow = 6;                  % window length in sec
            parameterTempogram.maxLag = 2;                    % corresponding to 30 bpm
            parameterTempogram.minLag = 0.1;                  % corresponding to 600 bpm

            [tempogram, T, Lag] = noveltyCurve_to_tempogram_via_ACF(noveltyCurve,parameterTempogram);
            %T = T./2;
            BPM_in = 60./Lag;
            BPM_out = 30:1:600;
            [ACFtempogram,~] = rescaleTempoAxis(tempogram,BPM_in,BPM_out);
            %figure();
           % visualize_tempogram(tempogram,T,BPM);
           
            parameterTempogram = [];
            parameterTempogram.featureRate = featureRate;
            parameterTempogram.tempoWindow = 6;         % window length in sec
            parameterTempogram.BPM = 30:1:600;          % tempo values

            [Ftempogram, ~, BPM] = noveltyCurve_to_tempogram_via_DFT(noveltyCurve,parameterTempogram);
             realF = abs(Ftempogram);
             realACF = abs(ACFtempogram);
           %{
            realF = abs(Ftempogram);
            realF = realF./norm(realF);
            realF(:,2:2:size(realF,2)) = (realF(:,2:2:size(realF,2)) + realF(:,1:2:size(realF,2)-1))/2;
            realF(:,1:2:size(realF,2)) = [];
            realACF = abs(ACFtempogram);
            realACF = realACF ./norm(realACF);
            if size(realF,2) > size(realACF,2),
               sss = (size(realF,2) - size(realACF,2))-1;
               realF(:, (size(realF,2)-sss):size(realF,2)) = [];
            elseif size(realF,2) < size(realACF,2),
                sss = (size(realACF,2) - size(realF,2)) -1;
               realACF(:, (size(realACF,2)-sss):size(realACF,2)) = [];
            end
            %}
            
            realTempogram = realF.*realACF;
            
            
             realSum = sum(realTempogram,2);
            [MAX Index1] = max(realSum);
            t1 = int32(BPM(Index1));
            range = double(t1 * 0.25);
            downLimit = double(t1 - range);

            upLimit = double(t1 + range);
            [ ~ ,li] = min(abs(BPM(:) - (downLimit)));
            [ ~ ,ui] = min(abs(BPM(:) - (upLimit)));
            li = int32(li); ui = int32(ui);
            realSum(li:ui) = 0;
           % realSum(MAXI) = 0;
            
            [SUB Index2] = max(realSum);
            t2 = int32(BPM(Index2));
            
            
            if t1> t2,
               tt = t1;
               t1 = t2;
               t2 = tt;
               tt = Index1;
               Index1 = Index2;
               Index2 = tt;
            end
            
            [~,FName ~]  = fileparts(filenames(j).name);
            G = textread(strcat('BallroomAnnotations/ballroomGroundTruth/',FName,'.bpm'));
            genreT21 = genreT21 +  (double(t2) / double(t1));
            genreT1Q = genreT1Q + (double(t1) / double(G));
            genreT2Q = genreT2Q + (double(t2) / double(G));
           % disp(strcat('t1 = ',num2str(t1),' t2= ',num2str(t2),' G= ',num2str(G)));
             %{
             if ((abs((double(t1)/3) - G)/G) <0.08) || ((abs((double(t2)/3) - G)/G) < 0.08),
                fprintf('%d   %d   %d   %f   %f  \n',t1/3,t2/3,G,abs((t1/3) - G)/G), ((abs((t2/3) - G)/G)  ); 
                count = count +1; 
             end
             %}
             [ P ALOTC] = estimateScore2(t1,t2,G,realTempogram, Index1 , Index2,genreList{k});
             [ P2 ALOTC2] = estimateScore2(t1/2,t2/2,G,realTempogram,Index1 , Index2,genreList{k});
             [ P3 ALOTC3] = estimateScore2(t1/3,t2/3,G,realTempogram,Index1 , Index2,genreList{k});
             [ P4 ALOTC4] = estimateScore2(t1/4,t2/4,G,realTempogram,Index1 , Index2,genreList{k});
             [ P5 ALOTC5] = estimateScore2(t1*2,t2*2,G,realTempogram,Index1 , Index2,genreList{k});
             [ P6 ALOTC6] = estimateScore2(t1*3,t2*3,G,realTempogram,Index1 , Index2,genreList{k});
             [ P7 ALOTC7] = estimateScore2(t1*4,4*t2,G,realTempogram,Index1 , Index2,genreList{k});

               genreP(1) = genreP(1) + P;
               genreP(2) = genreP(2) + P2;
               genreP(3) = genreP(3) + P3;
               genreP(4) = genreP(4) + P4;
               genreP(5) = genreP(5) + P5;
               genreP(6) = genreP(6) + P6;
               genreP(7) = genreP(7) + P7;
               genreALOTC(1) = genreALOTC(1) + ALOTC ;
               genreALOTC(2) = genreALOTC(2) + ALOTC2 ;
               genreALOTC(3) = genreALOTC(3) + ALOTC3 ;
               genreALOTC(4) = genreALOTC(4) + ALOTC4 ;
               genreALOTC(5) = genreALOTC(5) + ALOTC5 ;
               genreALOTC(6) = genreALOTC(6) + ALOTC6 ;
               genreALOTC(7) = genreALOTC(7) + ALOTC7 ;
             
               
       end
       
       genreP(:) = genreP(:) ./ length(filenames);
       genreALOTC(:) = genreALOTC(:) ./ length(filenames);
       genreT21 = genreT21 / length(filenames);
       genreT1Q = genreT1Q / length(filenames);
       genreT2Q = genreT2Q / length(filenames);
       
         fprintf(fileID,strcat(genreList{k},' P-Score: ',num2str(genreP(1)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/2): ',num2str(genreP(2)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/3): ',num2str(genreP(3)),'\n'));
        fprintf(fileID,strcat(genreList{k},' P-Score(T/4): ',num2str(genreP(4)),'\n'));
        fprintf(fileID,strcat(genreList{k},' ALOTC-Score: ',num2str(genreALOTC(1)),'\n'));
        fprintf(fileID,strcat(genreList{k},' T2/T1 : ',num2str(genreT21),'\n'));
        fprintf(fileID,strcat(genreList{k},' T1/Q : ',num2str(genreT1Q),'\n'));
        fprintf(fileID,strcat(genreList{k},' T2/Q : ',num2str(genreT2Q),'\n'));
      
        fprintf(fileID,'----------------------------- \n');
       
       
    end
    
    fprintf('\n\n');
 
 
 %%  Beat Tracking
 clear all;
    genreList = cell(8,1);
    
    genreList{1} = 'ChaChaCha';
    genreList{2} = 'Jive' ;
    genreList{3} = 'Quickstep' ;
    genreList{4} = 'Rumba' ;
    genreList{5} = 'Samba' ;
    genreList{6} = 'Tango' ;
    genreList{7} = 'VienneseWaltz' ;
    genreList{8} = 'Waltz' ;
    
    fileID = fopen('beatTracking.txt','w');
    for k=1:length(genreList),
        
       filenames = dir(strcat('BallroomData/',genreList{k},'/'));
       filenames(1:2) = [];
       dirWav = strcat('BallroomData/',genreList{k},'/');
       genrePrecision = 0;
       genreRecall = 0;
       genreFScore = 0;
       for j=1:length(filenames),
                
            [audio, fs ] = audioread(strcat(dirWav,filenames(j).name));
            Fs = fs;
            
            parameterNovelty = [];
            [noveltyCurve,featureRate] = audio_to_noveltyCurve(audio, Fs, parameterNovelty);
            
            parameterTempogram = [];
            parameterTempogram.featureRate = featureRate;
            parameterTempogram.tempoWindow = 6;         % window length in sec
            parameterTempogram.BPM = 30:1:600;          % tempo values

            [tempogram, T, BPM] = noveltyCurve_to_tempogram_via_DFT(noveltyCurve,parameterTempogram);
            
            
            
            parameterPLP = [];
            parameterPLP.featureRate = featureRate;
            parameterPLP.tempoWindow = parameterTempogram.tempoWindow;

            [PLP,featureRate] = tempogram_to_PLPcurve(tempogram, T, BPM, parameterPLP);
            PLP = PLP(1:length(noveltyCurve));  % PLP curve will be longer (zero padding)
           
            pos = [PLP,PLP(end)] > [PLP(1),PLP];
            neg = ~pos;
            peaks = find(pos(1:end-1).*neg(2:end));
            estimateBeatTime = [];
            
            for idx = 1:length(peaks)
                start = (floor(peaks(idx)./featureRate*Fs)+1) / Fs;
                estimateBeatTime = [estimateBeatTime ; start ]; 
            end
            
             [~,FName ~]  = fileparts(filenames(j).name);
            beatID = fopen(strcat('BallroomAnnotations/beats/',FName,'.beats'),'r');
            beats = [];
            while ~feof(beatID),
             beatTime = fscanf(beatID,'%f',1);
             [~] = fscanf(beatID,'%d',1);
             beats = [beats ; beatTime];
            end
            TP = 0; FP = 0; FN = 0;
            %disp(num2str((length(estimateBeatTime) - length(beats))));
            for idx = 1 : length(estimateBeatTime),
              if length(beats)>=1,
                [minValue , minIndex ] = min(abs(beats(:) - estimateBeatTime(idx)));
                if minValue<0.07,
                   TP = TP + 1;
                   beats(minIndex) = [];
                else
                   FP = FP + 1;
                end
              else
                   FP = FP + 1;
              end
            end
            FN = length(beats);
            Precision = double(TP) / double(TP+FP);
            Recall = double(TP) / double(TP+FN);
            FScore = (2 * Precision * Recall) / (Precision + Recall);
            genrePrecision = genrePrecision + Precision;
            genreRecall = genreRecall + Recall;
            genreFScore = genreFScore + FScore;
       end
        genrePrecision = genrePrecision / length(filenames);
        genreRecall = genreRecall / length(filenames);
        genreFScore = genreFScore / length(filenames);
        fprintf(fileID,strcat(genreList{k},'  Precision: ',num2str(genrePrecision),'\n'));
        fprintf(fileID,strcat(genreList{k},'  Recall: ',num2str(genreRecall),'\n'));
        fprintf(fileID,strcat(genreList{k},'  F-Score: ',num2str(genreFScore),'\n'));
        
    end
 
    fprintf('\n\n');