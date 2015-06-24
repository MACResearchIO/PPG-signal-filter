clear; clc; close all;
 
addpath(genpath('PPG'));

% Test Dataset IDs
testID = { 'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
       'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
       'TEST_S07_T02', 'TEST_S08_T01'};  
DataID = { 'DATA_01_TYPE01', 'DATA_02_TYPE02', 'DATA_03_TYPE02', 'DATA_04_TYPE02', ...
   'DATA_05_TYPE02', 'DATA_06_TYPE02', 'DATA_07_TYPE02', 'DATA_08_TYPE02',...
   'DATA_09_TYPE02', 'DATA_10_TYPE02', 'DATA_11_TYPE02', 'DATA_12_TYPE02'};  
BPMID = { 'Trace1', 'Trace2', 'Trace3', 'Trace4', ...
   'Trace5', 'Trace6', 'Trace7', 'Trace8',...
   'Trace9', 'Trace10', 'Trace11', 'Trace12'}; 
resultID = { 'Result_S01_T01', 'Result_S02_T02', 'Result_S03_T02', 'Result_S04_T02', ...
   'Result_S05_T02', 'Result_S06_T02', 'Result_S07_T02', 'Result_S08_T02',...
   'Result_S09_T02', 'Result_S10_T02','Result_S11_T02','Result_S12_T02'};    
         
   
for idnb = 1:12
    windowMaxFFTAmplitudeArray = zeros;
    windowMaxFFTindexArray = zeros;
    load(DataID{idnb});                      % load test dataset
    display(idnb);
    srate = 125;                             % 125 Hz
    window   = 8 * srate;                    % window length is 8 seconds
    currentWindow = 1;                      % initialize previous window
    step     = 2 * srate;                    % step size is 2 seconds
    windowNb = floor((length(sig)-window)/step)+1;  % total number of windows(estimates)
    %**********************************************************************
    % Please write your codes as follows (i.e.,inputing data window by window)
    BPM = zeros(windowNb,1);
    for i =   1  :  windowNb
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        [BPM(i,1), windowMaxFFTindexArray, windowMaxFFTAmplitudeArray] = MACVSSA(sig(1,curSegment), sig(2,curSegment), windowMaxFFTindexArray, windowMaxFFTAmplitudeArray, i);
    end
        save(resultID{idnb},'BPM');          % load test dataset 
end

for i = 1:12
   EKG = load(BPMID{i});
   PPG = load(resultID{i}); 
   EKG = EKG.BPM0;
   PPG = PPG.BPM;
   bs=(abs(PPG(1:numel(EKG),1)-EKG(1:numel(EKG),1))./EKG(1:numel(EKG),1))*100;
   Error(i) = mean(bs);
end

totalAvgError = mean(Error);

