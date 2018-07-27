%%
clear; clc; close all

%% first, creat a structure array for iterating thru the audio files
%

PATH_AUDIO = 'audio';

addpath(genpath('MIRtoolbox1.6.1'))

listOfSongsGuitar = listfile(fullfile(PATH_AUDIO,'guitar'))';
listOfSongsViolin = listfile(fullfile(PATH_AUDIO,'violin'))';
listOfSongsPiano = listfile(fullfile(PATH_AUDIO,'piano'))';
listOfSongsVoice = listfile(fullfile(PATH_AUDIO,'voice'))';
listOfSongsTest = listfile(fullfile(PATH_AUDIO,'test'))';

listOfSongsTrain = [listOfSongsGuitar;listOfSongsViolin;listOfSongsPiano;listOfSongsVoice];

%%  Q1

audio1 = miraudio(listOfSongsGuitar(2),'Normal');
a1 = mirgetdata(audio1);
time = (1:length(a1)) / 16;
frame1 = mirframe(audio1,'Length',1024,'sp','Hop',512,'sp');
spectrum = mirspectrum(frame1,'Max',2000,'Cents');
spectrum
frame2 = mirframe(audio1,'Length',2048,'sp','Hop',1024,'sp');
spectrum2 = mirspectrum(frame2,'Max',2000,'Cents');
spectrum2

%% Q2
audio21 = miraudio(listOfSongsViolin(88),'Normal');
audio22 = miraudio(listOfSongsVoice(88),'Normal');
frame21 = mirframe(audio21,'Length',1024,'sp','Hop',512,'sp');
frame22 = mirframe(audio22,'Length',1024,'sp','Hop',512,'sp');
spectrum21 = mirspectrum(frame21,'Max',2000);
spectrum22 = mirspectrum(frame22,'Max',2000);
spectrum21
spectrum22
%% Q2 - Compare
mircentroid(spectrum21)
mircentroid(spectrum22)
mirregularity(spectrum21)
mirregularity(spectrum22)
mirbrightness(spectrum21)
mirbrightness(spectrum22)
mirspread(spectrum21)
mirspread(spectrum22)
mirflatness(spectrum21)
mirflatness(spectrum22)
%mirzerocross(spectrum21)
%mirzerocross(spectrum22)
%[ mfccMean21  mfccSTD21] = extractMFCC(listOfSongsViolin(3),1024,512);
%[ mfccMean22  mfccSTD22] = extractMFCC(listOfSongsVoice(3),1024,512);
mirmfcc(spectrum21)
mirmfcc(spectrum22)
%%  Q3 - Extract Feature
addpath(genpath('MIRtoolbox1.6.1'))
w = 1024;
h = 512;

PATH_FEAT = 'feature';
X = [];

for i=1:800
    
    disp(i)
    fn = listOfSongsTrain{i};
    
    x = extractMFCC(fn,w,h);
    x = x';
    tempAudio = miraudio(fn,'Normal');
    tempFrame = mirframe(tempAudio,'Length',w,'sp','Hop',h,'sp');
    tempSpect = mirspectrum(tempFrame);
    
    temp  = mirgetdata(mirspread(tempSpect));
    %disp(temp);
    x = [x mean(temp,2) std(temp,0,2)];
   % disp(x);
    temp  = mirgetdata(mirflatness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mircentroid(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirbrightness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirroughness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirentropy(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirregularity(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirrolloff(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirattacktime(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirkurtosis(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirskewness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirattackslope(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirattackleap(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    %temp  = mirgetdata(mirzerocross(tempSpect));
    %x = [x mean(temp,2) std(temp,0,2)];
    
    
    % we can either save the extracted features as files
    [pathstr,fn,ext] = fileparts(fn);
    if i<=200
        fn = fullfile(PATH_FEAT,'guitar',[fn '.csv']);
    elseif i<=400
        fn = fullfile(PATH_FEAT,'violin',[fn '.csv']);
    elseif i<=600
        fn = fullfile(PATH_FEAT,'piano',[fn '.csv']);
    else
        fn = fullfile(PATH_FEAT,'voice',[fn '.csv']);
    end
    csvwrite(fn,x);
    
    % or we simply use a big matrix to keep them
    X = [X; x];
end
% save X X    % save it if you want


% extract the features for the test set by the way
% make sure you use the same features for the training and test sets

Xtest = [];
for i=1:200    
    disp(i)
    fn = listOfSongsTest{i};
    x = extractMFCC(fn,w,h);
    
    x = x';
    tempAudio = miraudio(fn,'Normal');
    tempFrame = mirframe(tempAudio,'Length',w,'sp','Hop',h,'sp');
    tempSpect = mirspectrum(tempFrame);
    
    temp  = mirgetdata(mirspread(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirflatness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mircentroid(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirbrightness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirroughness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirentropy(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirregularity(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirrolloff(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirattacktime(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirkurtosis(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirskewness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
     temp  = mirgetdata(mirattackslope(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirattackleap(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    %temp  = mirgetdata(mirzerocross(tempSpect));
    %x = [x mean(temp,2) std(temp,0,2)];
    
    Xtest = [Xtest; x];
end
% save Xtest Xtest  


% set the groundtruth
% 1: guitar 
% 2: violin
% 3: piano
% 4: voice
Y = [ones(200,1);2*ones(200,1);3*ones(200,1);4*ones(200,1)];
% save Y Y 

% clear some variables to keep the work space neat
 % clearvars -except X Y Xtest
 
 %% Transfer the NaNs into mean
 
 
for i = 1 : size(X,1),
    for j = 1 : size(X,2),
        if isnan(X(i,j)),
            X(i,j) =  mean(X(:,j),'omitnan');
        end 
    end
end
%%
for i = 1 : size(Xtest,1),
    for j = 1 : size(Xtest,2),
        if isnan(Xtest(i,j)),
            Xtest(i,j) =  mean(Xtest(:,j),'omitnan');
        end 
    end
end
%%
 MeanGuitar = mean(X (1:200,:));
 MeanViolin = mean(X (201:400,:));
 MeanPiano = mean(X (401:600,:));
 MeanVoice = mean(X (601:800,:));

%% split the data into training and validation (note: this is not cross validation)

nValidation = floor(size(X,1)/5); % use 20% of the data for validation
nTrain = size(X,1) - nValidation;
nTest = size(Xtest,1);

rp = randperm(800); % generate a random permutation
indexValidation = rp(1:nValidation);
indexTrain = rp(nValidation+1:end);

Xvalidation = X(indexValidation,:);
Yvalidation = Y(indexValidation,:);
Xtrain = X(indexTrain,:);
Ytrain = Y(indexTrain,:);



%% feature normalization   -- In Q4 , Skip this step

% the values for normalization are computed from the training set
featMean = mean(Xtrain);
featSTD = std(Xtrain);


% then apply to the training, validation, and test sets
Xtrain = (Xtrain - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
Xvalidation = (Xvalidation - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

Xtest2 = (Xtest - repmat(featMean,nTest,1))./(repmat(featSTD,nTest,1)+eps);


%Xtrain(isnan(Xtrain)) = 0;
%Xvalidation(isnan(Xvalidation)) = 0;
%Xtest(isnan(Xtest)) = 0;

%% train the classifier
%addpath('./libsvm-3.21/matlab');
% svmtrain
% svmpredict



c0 = 1;
g0 = 1/size(X,2);

Cs = [c0 c0*10 c0*100 c0*1000]; % possible range of the parameter C
Gs = [g0 g0/10 g0/100]; % possible range of the parameter gamma

% default
bestAccuraccy = 0.25; 
bestModel = {};
bestC = nan;
bestG = nan;

for c=1:length(Cs)
    for g=1:length(Gs)
        
        %model = svmtrain(Ytrain,Xtrain,sprintf('-t 2 -c %f -g %f',Cs(c),Gs(g)));
        model = svmtrain(Ytrain,Xtrain,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g))); % quiet mode
        % actually, you can also use svmtrain(...,'-v 5') to implement 5-fold 
        % cross validation, but we are not using that in this code
        
        % Yvaliation is the groundtruth
        % Ypred is the prediction result
        [Ypred, accuracy, ~] = svmpredict(Yvalidation, Xvalidation, model, '-q');
        accuracy = accuracy(1); % the first one correponds to classification accuracy
                                % accuracy = sum(Ypred==Yvalidation)/length(Yvalidation)
        disp(sprintf('c=%f g=%f accuracy=%f',Cs(c),Gs(g),accuracy))

        if accuracy > bestAccuraccy
            bestAccuraccy = accuracy;
            bestModel = model;
            bestC = Cs(c);
            bestG = Gs(g);
        end
    end
end

% from Ypred and Yvalidation you can create the confusion tlabe
[Ypred, accuracy, ~] = svmpredict(Yvalidation, Xvalidation, bestModel);
ConfusionTable = zeros(4,4);

totalGuess = zeros(4,1);
totalWhat = zeros(4,1);
Precision = zeros(4,1);
Recall = zeros(4,1);
FScore = zeros(4,1);
for i=1:length(Ypred),
ConfusionTable(Yvalidation(i),Ypred(i)) = ConfusionTable(Yvalidation(i),Ypred(i)) +1 ;  
totalGuess(Ypred(i)) = totalGuess(Ypred(i)) + 1;
totalWhat(Yvalidation(i)) = totalWhat(Yvalidation(i)) + 1;
end

for i=1:4,
  Precision(i) =  ConfusionTable(i,i) / totalGuess(i);
  Recall(i) = ConfusionTable(i,i) / totalWhat(i);
  FScore(i) = 2* Precision(i) * Recall(i) / (Precision(i) + Recall(i));
end

showConfusionTable(ConfusionTable,Recall,Precision , FScore);


%%
YtestPred = svmpredict(zeros(nTest,1), Xtest2, bestModel,'-q');
csvwrite('YtestPred.csv',YtestPred);



%%  Q6 - Add other Feature
addpath(genpath('MIRtoolbox1.6.1'))
w = 1024;
h = 512;

PATH_FEAT = 'mFeature';
X = [];
XSpread = [];  %40~50
XIrregular = [];  % 30~40
XBrightness = []; % 35 ~ 45
XFlatness = []; % 40~50
XCentroid = []; %40~50
XRoughness = [];
XBcdf = [];
XRolloff = [];
XKurtosis = [];
XZerocross = [];
for i=1:800
    
    disp(i)
    fn = listOfSongsTrain{i};
    
    x = extractMFCC(fn,w,h);
    x = x';
    tempAudio = miraudio(fn,'Normal');
    tempFrame = mirframe(tempAudio,'Length',w,'sp','Hop',h,'sp');
    tempSpect = mirspectrum(tempFrame);
   
    
    temp  = mirgetdata(mirentropy(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirspread(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
    temp  = mirgetdata(mirroughness(tempSpect));
   x = [x mean(temp,2) std(temp,0,2)];
     temp  = mirgetdata(mirbrightness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
     temp  = mirgetdata(mirroughness(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
     temp  = mirgetdata(mirskewness(tempSpect));
   x = [x mean(temp,2) std(temp,0,2)];
     temp  = mirgetdata(mirrolloff(tempSpect));
     x = [x mean(temp,2) std(temp,0,2)];
     temp  = mirgetdata(mirzerocross(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
     temp  = mirgetdata(mirkurtosis(tempSpect));
     x = [x mean(temp,2) std(temp,0,2)];
      temp  = mirgetdata(mirregularity(tempSpect));
    x = [x mean(temp,2) std(temp,0,2)];
     
     
    [pathstr,fn,ext] = fileparts(fn);
    if i<=200
        fn = fullfile(PATH_FEAT,'guitar',[fn '.csv']);
    elseif i<=400
        fn = fullfile(PATH_FEAT,'violin',[fn '.csv']);
    elseif i<=600
        fn = fullfile(PATH_FEAT,'piano',[fn '.csv']);
    else
        fn = fullfile(PATH_FEAT,'voice',[fn '.csv']);
    end
    csvwrite(fn,x);
    
    % or we simply use a big matrix to keep them
    X = [X; x];
  
    XCentroid = [XCentroid; x(27:28)];
    XSpread = [XSpread; x(29:30)];
    XFlatness = [XFlatness; x(31:32)];
    XBrightness = [ XBrightness ; x(33:34)];
    XRoughness = [ XRoughness;x(35:36)];
    XBcdf = [ XBcdf;x(37:38)];
    XRolloff = [ XRolloff;x(39:40)];
    XZerocross = [XZerocross;x(41:42)];
    XKurtosis = [XKurtosis;x(43:44)];
    XIrregular = [XIrregular;x(45:46)];
    
    %{
    XCentroid = [XCentroid; x(2)];
    XSpread = [XSpread; x(3)];
    XFlatness = [XFlatness; x(4)];
    XBrightness = [ XBrightness ; x(5)];
    %}
end
% save X X    % save it if you want



for i = 1 : size(X,1),
    for j = 1 : size(X,2),
        if isnan(X(i,j)),
            X(i,j) =  mean(X(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XCentroid,1),
    for j = 1 : size(XCentroid,2),
        if isnan(XCentroid(i,j)),
            XCentroid(i,j) =  mean(XCentroid(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XSpread,1),
    for j = 1 : size(XSpread,2),
        if isnan(XSpread(i,j)),
            XSpread(i,j) =  mean(XSpread(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XFlatness,1),
    for j = 1 : size(XFlatness,2),
        if isnan(XFlatness(i,j)),
            XFlatness(i,j) =  mean(XFlatness(:,j),'omitnan');
        end 
    end
end

for i = 1 : size(XRoughness,1),
    for j = 1 : size(XRoughness,2),
        if isnan(XRoughness(i,j)),
            XRoughness(i,j) =  mean(XRoughness(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XBcdf,1),
    for j = 1 : size(XBcdf,2),
        if isnan(XBcdf(i,j)),
            XBcdf(i,j) =  mean(XBcdf(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XRolloff,1),
    for j = 1 : size(XRolloff,2),
        if isnan(XRolloff(i,j)),
            XRolloff(i,j) =  mean(XRolloff(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XBrightness,1),
    for j = 1 : size(XBrightness,2),
        if isnan(XBrightness(i,j)),
            XBrightness(i,j) =  mean(XBrightness(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XKurtosis,1),
    for j = 1 : size(XKurtosis,2),
        if isnan(XKurtosis(i,j)),
            XKurtosis(i,j) =  mean(XKurtosis(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XZerocross,1),
    for j = 1 : size(XZerocross,2),
        if isnan(XZerocross(i,j)),
            XZerocross(i,j) =  mean(XZerocross(:,j),'omitnan');
        end 
    end
end
for i = 1 : size(XIrregular,1),
    for j = 1 : size(XIrregular,2),
        if isnan(XIrregular(i,j)),
            XIrregular(i,j) =  mean(XIrregular(:,j),'omitnan');
        end 
    end
end
% extract the features for the test set by the way
% make sure you use the same features for the training and test sets

%{

Xtest = [];
XtestCentroid = [];
XtestSpread = [];
XFlatness = [];
for i=1:200    
    
    disp(i)
    fn = listOfSongsTest{i};
    %x = extractMFCC(fn,w,h);
    
    tempAudio = miraudio(fn,'Normal');
    tempFrame = mirframe(tempAudio,'Length',w,'sp','Hop',h,'sp');
    tempSpect = mirspectrum(tempFrame);
    temp  = mirgetdata(mircentroid(tempSpect));
    x = [x;temp'];
    temp  = mirgetdata(mirspread(tempSpect));
    x = [x;temp'];
    temp  = mirgetdata(mirflatness(tempSpect));
    x = [x;temp'];
     temp  = mirgetdata(mirbrightness(tempSpect));
     x = [x;temp'];
    
    Xtest = [Xtest; x'];
    
end
% save Xtest Xtest  

%}
% set the groundtruth
% 1: guitar 
% 2: violin
% 3: piano
% 4: voice
Y = [ones(200,1);2*ones(200,1);3*ones(200,1);4*ones(200,1)];
% save Y Y 

%%
X = X(:,1:26);
%% split the data into training and validation (note: this is not cross validation)

nValidation = floor(size(X,1)/5); % use 20% of the data for validation
nTrain = size(X,1) - nValidation;
%nTest = size(Xtest,1);

rp = randperm(800); % generate a random permutation
indexValidation = rp(1:nValidation);
indexTrain = rp(nValidation+1:end);

Xvalidation = X(indexValidation,:);
Xtrain = X(indexTrain,:);

XvalidationCentroid = XCentroid(indexValidation,:);
XtrainCentroid = XCentroid(indexTrain,:);

XvalidationSpread = XSpread(indexValidation,:);
XtrainSpread = XSpread(indexTrain,:);

XvalidationFlatness = XFlatness(indexValidation,:);
XtrainFlatness = XFlatness(indexTrain,:);

XvalidationBrightness = XBrightness(indexValidation,:);
XtrainBrightness = XBrightness(indexTrain,:);

XvalidationRoughness = XRoughness(indexValidation,:);
XtrainRoughness = XRoughness(indexTrain,:);

XvalidationBcdf = XBcdf(indexValidation,:);
XtrainBcdf = XBcdf(indexTrain,:);

XvalidationRolloff = XRolloff(indexValidation,:);
XtrainRolloff = XRolloff(indexTrain,:);

XvalidationKurtosis = XKurtosis(indexValidation,:);
XtrainKurtosis = XKurtosis(indexTrain,:);

XvalidationZerocross = XZerocross(indexValidation,:);
XtrainZerocross = XZerocross(indexTrain,:);

XvalidationIrregular = XIrregular(indexValidation,:);
XtrainIrregular = XIrregular(indexTrain,:);

Yvalidation = Y(indexValidation,:);
Ytrain = Y(indexTrain,:);

%% feature normalization   -- In Q4 , Skip this step

% the values for normalization are computed from the training set
featMean = mean(Xtrain);
featSTD = std(Xtrain);
% then apply to the training, validation, and test sets
Xtrain = (Xtrain - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
Xvalidation = (Xvalidation - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

%Xtest = (Xtest - repmat(featMean,nTest,1))./(repmat(featSTD,nTest,1)+eps);


featMean = mean(XtrainCentroid);
featSTD = std(XtrainCentroid);
XtrainCentroid = (XtrainCentroid - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationCentroid = (XvalidationCentroid - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

featMean = mean(XtrainSpread);
featSTD = std(XtrainSpread);
XtrainSpread = (XtrainSpread - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationSpread = (XvalidationSpread - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

featMean = mean(XtrainFlatness);
featSTD = std(XtrainFlatness);
XtrainFlatness = (XtrainFlatness - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationFlatness = (XvalidationFlatness - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);


featMean = mean(XtrainBrightness);
featSTD = std(XtrainBrightness);
XtrainBrightness = (XtrainBrightness - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationBrightness = (XvalidationBrightness - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

featMean = mean(XtrainRoughness);
featSTD = std(XtrainRoughness);
XtrainRoughness = (XtrainRoughness - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationRoughness = (XvalidationRoughness - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

featMean = mean(XtrainBcdf);
featSTD = std(XtrainBcdf);
XtrainBcdf = (XtrainBcdf - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationBcdf = (XvalidationBcdf - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

featMean = mean(XtrainRolloff);
featSTD = std(XtrainRolloff);
XtrainRolloff = (XtrainRolloff - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationRolloff = (XvalidationRolloff - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);


featMean = mean(XtrainKurtosis);
featSTD = std(XtrainKurtosis);
XtrainKurtosis = (XtrainKurtosis - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationKurtosis = (XvalidationKurtosis - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);

featMean = mean(XtrainZerocross);
featSTD = std(XtrainZerocross);
XtrainZerocross = (XtrainZerocross - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationZerocross = (XvalidationZerocross - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);


featMean = mean(XtrainIrregular);
featSTD = std(XtrainIrregular);
XtrainIrregular = (XtrainIrregular - repmat(featMean,nTrain,1))./(repmat(featSTD,nTrain,1)+eps);
XvalidationIrregular = (XvalidationIrregular - repmat(featMean,nValidation,1))./(repmat(featSTD,nValidation,1)+eps);
%XtrainSpread(isnan(Xtrain)) = 0;
%Xvalidation(isnan(Xvalidation)) = 0;
%Xtest(isnan(Xtest)) = 0;

%% train the classifier  with own features
%addpath('./libsvm-3.21/matlab');
% svmtrain
% svmpredict



c0 = 1;
g0 = 1/size(X,2);

Cs = [c0 c0*10 c0*100 c0*1000]; % possible range of the parameter C
Gs = [g0 g0/10 g0/100]; % possible range of the parameter gamma

% default
bestAccuraccy = 0.25; 
bestModel = {};
bestC = nan;
bestG = nan;
bestPred = zeros(200,1);

for c=1:length(Cs)
    for g=1:length(Gs)
        
        %model = svmtrain(Ytrain,Xtrain,sprintf('-t 2 -c %f -g %f',Cs(c),Gs(g)));
        model = svmtrain(Ytrain,Xtrain,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g))); % quiet mode
        model2 = svmtrain(Ytrain,XtrainCentroid, sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model3 = svmtrain(Ytrain,XtrainSpread, sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model4 = svmtrain(Ytrain,XtrainFlatness,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model5 = svmtrain(Ytrain,XtrainBrightness,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model6 = svmtrain(Ytrain,XtrainRoughness,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model7 = svmtrain(Ytrain,XtrainBcdf,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model8 = svmtrain(Ytrain,XtrainRolloff,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model9 = svmtrain(Ytrain,XtrainKurtosis,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model10 = svmtrain(Ytrain,XtrainZerocross,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        model11 = svmtrain(Ytrain,XtrainIrregular,sprintf('-t 2 -c %f -g %f -q',Cs(c),Gs(g)));
        % actually, you can also use svmtrain(...,'-v 5') to implement 5-fold 
        % cross validation, but we are not using that in this code
        
        cumPred = zeros(160,11);
        
        % Yvaliation is the groundtruth
        % Ypred is the prediction result
        [Ypred, accuracy, ~] = svmpredict(Yvalidation, Xvalidation, model, '-q');
        accuracy = accuracy(1); % the first one correponds to classification accuracy
                                % accuracy = sum(Ypred==Yvalidation)/length(Yvalidation) 
         disp(sprintf('c=%f g=%f accuracy=%f',Cs(c),Gs(g),accuracy))
         [Ypred2, accuracy2, ~] = svmpredict(Yvalidation, XvalidationCentroid, model2, '-q');
        accuracy2 = accuracy2(1);
        disp(sprintf('c=%f g=%f accuracy2=%f',Cs(c),Gs(g),accuracy2))
        
        [Ypred3, accuracy3, ~] = svmpredict(Yvalidation, XvalidationSpread, model3, '-q');
        accuracy3 = accuracy3(1);
        disp(sprintf('c=%f g=%f accuracy3=%f',Cs(c),Gs(g),accuracy3))
        
        [Ypred4, accuracy4, ~] = svmpredict(Yvalidation, XvalidationFlatness, model4, '-q');
        accuracy4 = accuracy4(1);
        disp(sprintf('c=%f g=%f accuracy4=%f',Cs(c),Gs(g),accuracy4))
        
        [Ypred5, accuracy5, ~] = svmpredict(Yvalidation, XvalidationBrightness, model5, '-q');
        accuracy5 = accuracy5(1);
        disp(sprintf('c=%f g=%f accuracy5=%f',Cs(c),Gs(g),accuracy5))
        
        [Ypred6, accuracy6, ~] = svmpredict(Yvalidation, XvalidationRoughness, model6, '-q');
        accuracy6 = accuracy6(1);
        disp(sprintf('c=%f g=%f accuracy6=%f',Cs(c),Gs(g),accuracy6))
        
        [Ypred7, accuracy7, ~] = svmpredict(Yvalidation, XvalidationBcdf, model7, '-q');
        accuracy7 = accuracy7(1);
        disp(sprintf('c=%f g=%f accuracy7=%f',Cs(c),Gs(g),accuracy7))
        
        [Ypred8, accuracy8, ~] = svmpredict(Yvalidation, XvalidationRolloff, model8, '-q');
        accuracy8 = accuracy8(1);
        disp(sprintf('c=%f g=%f accuracy8=%f',Cs(c),Gs(g),accuracy8))
        
         [Ypred9, accuracy9, ~] = svmpredict(Yvalidation, XvalidationKurtosis, model9, '-q');
        accuracy9 = accuracy9(1);
        disp(sprintf('c=%f g=%f accuracy9=%f',Cs(c),Gs(g),accuracy9))
        
         [Ypred10, accuracy10, ~] = svmpredict(Yvalidation, XvalidationZerocross, model10, '-q');
        accuracy10 = accuracy10(1);
        disp(sprintf('c=%f g=%f accuracy10=%f',Cs(c),Gs(g),accuracy10))
        
         [Ypred11, accuracy11, ~] = svmpredict(Yvalidation, XvalidationIrregular, model11, '-q');
        accuracy11 = accuracy11(1);
        disp(sprintf('c=%f g=%f accuracy11=%f',Cs(c),Gs(g),accuracy11))
        
        finalPred = zeros(length(Ypred),1);
        
        for i=1:length(Ypred),
           cumPred(i,Ypred(i)) =   cumPred(i,Ypred(i)) + accuracy;
           cumPred(i,Ypred2(i)) =   cumPred(i,Ypred2(i)) + accuracy2;
           cumPred(i,Ypred3(i)) =   cumPred(i,Ypred3(i)) + accuracy3;
           cumPred(i,Ypred4(i)) =   cumPred(i,Ypred4(i)) + accuracy4;
           cumPred(i,Ypred5(i)) =   cumPred(i,Ypred5(i)) + accuracy5;
           cumPred(i,Ypred6(i)) =   cumPred(i,Ypred6(i)) + accuracy6;
           cumPred(i,Ypred7(i)) =   cumPred(i,Ypred7(i)) + accuracy7;
           cumPred(i,Ypred8(i)) =   cumPred(i,Ypred8(i)) + accuracy8;
           cumPred(i,Ypred9(i)) =   cumPred(i,Ypred9(i)) + accuracy9;
           cumPred(i,Ypred10(i)) =   cumPred(i,Ypred10(i)) + accuracy10;
           cumPred(i,Ypred11(i)) =   cumPred(i,Ypred11(i)) + accuracy11;
           m = 0;
           for j=1:11,
               if cumPred(i,j) > m ,
                  m = cumPred(i,j);
                  finalPred(i) = j ;
               end
           end
           
        end
        count = 0.0;
        for i=1:length(Ypred),
            if finalPred(i) == Yvalidation(i),
               count = count + 1; 
            end
        end
        finalAccuracy = count /length(Ypred);

        if finalAccuracy > bestAccuraccy
            bestAccuraccy = finalAccuracy;
            bestModel = model;
            bestC = Cs(c);
            bestG = Gs(g);
            bestPred = finalPred;
        end
    end
end
disp('QAQ');
% from Ypred and Yvalidation you can create the confusion tlabe
%[Ypred, accuracy, ~] = svmpredict(Yvalidation, Xvalidation, bestModel);
%ConfusionTable = zeros(4,4);

