clear
close all
clc


Param.PreStim           = 1;   % (ms) a negatove integer, how long before stimulus onset you want to start decoding
Param.PostStim          = 300; % (ms) a positive integer, how long before stimulus onset you want to start decoding
Param.MovWind           = 30;  % (ms)a positive integer, the width of time window for counting spikes
Param.NumOfTrailGroup   = 100; % number of trials for each oriention for decoding
Param.NumOfCells        = 40;  % how many neurons do you want to use for decoding
Param.NumCrossVal       = 10;  % how many cross validatoin you need
Param.PercentOfTest     = 0.2; % the percentation of test data
Param.NumUniqStim       = 6;   % number of orientation you want to decode (12 Orientation)
Param.TimResolution     = 2;   % temporal resolution for decoder
Param.Classifier        = 'lda'; % LDA or NBC (Naive Bayes Classifier), case insensitive
Param.FileName          = 'F:\CJ194\Data\CJ194_datafile028.nev';  % select a data file name and dicrctor
Param.StimulusName      = 'F:\CJ194\Stimulus\Paired_Stimulus_File_CJ194_0005.mat';
% Param.SelectedChannels  = [1,2,3,4,5,6,7,8,9,10,11,12,14,19,21,22,26,27,29,31,32,37,40,41,42,44,46,50,51,52,53,54,55,56,57,58,62,63,65,66,67,73,75,76,81,83,84,85,86,87,94,95];
% Param.SelectedChannels  = [1:14 16:19 21 22 26 27 29 31 32 37 40:42 44 46 47 50:58 61:70 73 75 76 81 83:88 91 93:96 ];
Param.SelectedChannels    = [1:12 17 19 21 23 26 27 29 32 32 37 40 41 42 44 46 50 51:57 66 73 75 76 83 85 86 87];
Param.select_Alignments = 1;
Param.NormalizeRF       = true;  % true will zscore the spike count
Param.SDF               = true;  % feed the classifer with FR instead of spike count
Param.SDF_binSize       = 15;    % ms
%% read cells and spikes
rng(1); % For reproducibility

% make some zore array to stor the results
Group_Cnt                 = 1;
accyracy_Matrix           = zeros(Param.NumCrossVal ,length(Param.PreStim : Param.TimResolution : Param.PostStim- Param.MovWind), 400/Param.NumOfTrailGroup);
raw_SpkieCount            = zeros(Param.NumCrossVal ,length(Param.PreStim : Param.TimResolution : Param.PostStim- Param.MovWind), 400/Param.NumOfTrailGroup);
raw_Mean_SpkieCount       = zeros(Param.NumCrossVal ,length(Param.PreStim : Param.TimResolution : Param.PostStim- Param.MovWind), 400/Param.NumOfTrailGroup);
normalize_SpkieCount      = zeros(Param.NumCrossVal ,length(Param.PreStim : Param.TimResolution : Param.PostStim- Param.MovWind), 400/Param.NumOfTrailGroup);
normalize_Mean_SpkieCount = zeros(Param.NumCrossVal ,length(Param.PreStim : Param.TimResolution : Param.PostStim- Param.MovWind), 400/Param.NumOfTrailGroup);

% read cells and spike
Param  = construct_FR_Date(Param);

% make Clsss labels
Labels = reshape(repmat([1:Param.NumUniqStim]',[1 Param.NumOfTrailGroup])', [], 1);
N      = length(Labels);

clc
% Decoding ...
for iGroup = 1 : Param.NumOfTrailGroup : 400
    
    time_Ind = 1;
    
    for iTime = Param.PreStim : Param.TimResolution : Param.PostStim - Param.MovWind
        
        fprintf(['\n Decoding trials: ' num2str(iGroup) '-' num2str(iGroup + Param.NumOfTrailGroup-1) ', Decoding  data at time: ' num2str(iTime) ' (ms)...'])
        fprintf( '\nCross Validation #: ' )
        
        % perform some corss validation
        for cr = 1 : Param.NumCrossVal
            
            fprintf([ num2str(cr) ' '])
            randomized_Cell_Order = randperm(Param.NumOfCells);  % randomize cell order
            
            spike_Count = [];
            for st = 1 : Param.NumUniqStim
                
                spike_Count = [spike_Count; squeeze(sum(Param.all_Resp{st}(iGroup : iGroup + Param.NumOfTrailGroup -1, iTime : Param.MovWind + iTime - 1, randomized_Cell_Order(1:Param.NumOfCells)), 2))];
                
            end
            raw_SpkieCount(cr, time_Ind, Group_Cnt)      = sum(spike_Count(:));
            raw_Mean_SpkieCount(cr, time_Ind, Group_Cnt) = mean(spike_Count(:));
            if Param.NormalizeRF == true
                spike_Count = zscore(spike_Count, [], 2);
            end
            normalize_SpkieCount(cr, time_Ind, Group_Cnt)      = sum(spike_Count(:));
            normalize_Mean_SpkieCount(cr, time_Ind, Group_Cnt) = mean(spike_Count(:));
            
            cvp     = cvpartition(N, 'Holdout', Param.PercentOfTest);
            idxTrn  = training(cvp); % Training set indices
            idxTest = test(cvp);    % Test set indices
            fprintf(['(' num2str(size(spike_Count, 1)) ', ' num2str(size(spike_Count, 2)) '), '])
            if strcmpi(Param.Classifier, 'lda')
                
                Mdl                                      = fitcdiscr(spike_Count(idxTrn, :), Labels(idxTrn, :));
                labels                                   = predict(Mdl, spike_Count(idxTest, :));
                accyracy_Matrix(cr, time_Ind, Group_Cnt) = 100*mean(labels==Labels(idxTest, :));
                  
            elseif strcmpi(Param.Classifier, 'nbc')
                
                Mdl                                      = fitcnb(spike_Count(idxTrn, :), Labels(idxTrn, :));
                labels                                   = predict(Mdl, spike_Count(idxTest, :));
                accyracy_Matrix(cr, time_Ind, Group_Cnt) = 100*mean(labels==Labels(idxTest, :));
           
            elseif strcmpi(Param.Classifier, 'svm')
                
                Mdl                                      = fitcecoc(spike_Count(idxTrn, :), Labels(idxTrn, :));
                labels                                   = predict(Mdl, spike_Count(idxTest, :));
                accyracy_Matrix(cr, time_Ind, Group_Cnt) = 100*mean(labels==Labels(idxTest, :));
                
            else
                error('**** The classifier should be either *LDA* or *NBC* (Naive Bayes Classifier)')
            end
            
            
        end
        
        time_Ind = time_Ind + 1;
        
    end
    
    Group_Cnt = Group_Cnt +1;
end
%%
c = colormap;

subplot(231)
color_Ind = floor(size(c, 1)/size(accyracy_Matrix, 3));
ci = 1;
for iGroup = 1 : size(accyracy_Matrix, 3)
    
    plot(mean(accyracy_Matrix(:, :, iGroup)), 'color',c(ci, :) )
    ci = ci + color_Ind;
    hold on
    
end
xlabel('Time (ms)')
ylabel('accuracy')

subplot(232)
ci = 1;
for iGroup = 1 : size(raw_SpkieCount, 3)
    
    plot(mean(raw_SpkieCount(:, :, iGroup)), 'color',c(ci, :) )
    ci = ci + color_Ind;
    hold on
    
end
xlabel('Time (ms)')
ylabel('spike count (raw)')

subplot(233)
ci = 1;
for iGroup = 1 : size(raw_Mean_SpkieCount, 3)
    
    plot(mean(raw_Mean_SpkieCount(:,:,iGroup)), 'color',c(ci, :) )
    ci = ci + color_Ind;
    hold on
    
end
xlabel('Time (ms)')
ylabel('mean of spike count (raw)')

subplot(234)
ci = 1;
for iGroup = 1 : size(normalize_Mean_SpkieCount, 3)
    
    plot(mean(normalize_SpkieCount(:, :, iGroup)), 'color',c(ci, :) )
    ci = ci + color_Ind;
    hold on
    
end
xlabel('Time (ms)')
ylabel('spike count (normalized)')

subplot(235)
ci = 1;
for iGroup = 1 : size(normalize_Mean_SpkieCount, 3)
    
    plot(mean(normalize_Mean_SpkieCount(:, :, iGroup)), 'color',c(ci, :) )
    ci = ci + color_Ind;
    hold on
    
end
xlabel('Time (ms)')
ylabel('mean spike count (normalized)')
