clear
close all
clc


Param.PreStim           = 1;   % (ms) a negatove integer, how long before stimulus onset you want to start decoding
Param.PostStim          = 300; % (ms) a positive integer, how long before stimulus onset you want to start decoding
Param.MovWind           = 30;  % (ms)a positive integer, the width of time window for counting spikes
Param.NumOfTrailGroup   = 100; % number of trials for each oriention for decoding
Param.NumOfCells        = 50;  % how many neurons do you want to use for decoding
Param.NumCrossVal       = 10;  % how many cross validatoin you need
Param.PercentOfTest     = 0.2; % the percentation of test data
Param.NumUniqStim       = 6;   % number of orientation you want to decode (12 Orientation)
Param.TimResolution     = 2;   % temporal resolution for decoder
Param.Classifier        = 'LDA'; % LDA or NBC (Naive Bayes Classifier), case insensitive
Param.FileName          = 'F:\CJ194\Data\CJ194_datafile025.nev';  % select a data file name and dicrctor
Param.StimulusName      = 'F:\CJ194\Stimulus\Paired_Stimulus_File_CJ194_0001.mat';
Param.SelectedChannels  = [1,2,3,4,5,6,7,8,9,10,11,12,14,19,21,22,26,27,29,31,32,37,40,41,42,44,46,50,51,52,53,54,55,56,57,58,62,63,65,66,67,73,75,76,81,83,84,85,86,87,94,95];
% Param.SelectedChannels  = [1:14 16:19 21 22 26 27 29 31 32 37 40:42 44 46 47 50:58 61:70 73 75 76 81 83:88 91 93:96 ];
Param.select_Alignments = 1;
Param.NormalizeRF       = true;  % true will zscore the spike count
Param.SDF               = true;  % feed the classifer with FR instead of spike count
Param.SDF_binSize       = 15;    % ms
Param.TarinTestCellRandomization = false;  % true will randomize test and train cells
%% read cells and spikes
rng(1); % For reproducibility

% make some zore array to stor the results
Group_Cnt                 = 1;
accuracy_Matrix           = zeros(length(Param.PreStim : Param.TimResolution : Param.PostStim- Param.MovWind) ,length(Param.PreStim : Param.TimResolution : Param.PostStim- Param.MovWind), 400/Param.NumOfTrailGroup);

% read cells and spike
Param  = construct_FR_Date(Param);

% make Clsss labels
Labels = reshape(repmat([1:Param.NumUniqStim]',[1 Param.NumOfTrailGroup])', [], 1);
N      = length(Labels);

clc
% Decoding ...
for iGroup = 1 : Param.NumOfTrailGroup : 400
    
    time_Train_Ind = 1;
    
    for iTimeTrain = Param.PreStim : Param.TimResolution : Param.PostStim - Param.MovWind
        
        time_Test_Ind  = 1;
        
        for iTimeTest = Param.PreStim : Param.TimResolution : Param.PostStim - Param.MovWind
            
            fprintf(['\n Decoding trials: ' num2str(iGroup) '-' num2str(iGroup + Param.NumOfTrailGroup-1) ', Decoding  data at time: ' num2str(iTimeTrain) ' (ms)...'])
            fprintf( '\nCross Validation #: ' )
            
            % perform some corss validation
            for cr = 1 : Param.NumCrossVal
                
                fprintf([ num2str(cr) ' '])
                if Param.TarinTestCellRandomization == true
                    randomized_Cell_Order_Train = randperm(Param.NumOfCells);  % randomize cell order for train
                    randomized_Cell_Order_Test  = randperm(Param.NumOfCells);  % randomize cell order for test
                else
                    
                    randomized_Cell_Order_Train = randperm(Param.NumOfCells);  % randomize cell order for train
                    randomized_Cell_Order_Test  = randomized_Cell_Order_Train;  % randomize cell order for test
                end
                spike_Count_Train = [];
                spike_Count_Test = [];
                for st = 1 : Param.NumUniqStim
                    
                    spike_Count_Train = [spike_Count_Train; squeeze(sum(Param.all_Resp{st}(iGroup : iGroup + Param.NumOfTrailGroup -1, iTimeTrain : Param.MovWind + iTimeTrain - 1, randomized_Cell_Order_Train(1:Param.NumOfCells)), 2))];                   
                    spike_Count_Test = [spike_Count_Test;   squeeze(sum(Param.all_Resp{st}(iGroup : iGroup + Param.NumOfTrailGroup -1, iTimeTest : Param.MovWind + iTimeTest - 1, randomized_Cell_Order_Test(1:Param.NumOfCells)), 2))];
                    
                end
                
                if Param.NormalizeRF == true
                    spike_Count_Train = zscore(spike_Count_Train, [], 2);
                    spike_Count_Test  = zscore(spike_Count_Test,  [], 2);
                end
                
                cvp     = cvpartition(N, 'Holdout', Param.PercentOfTest);
                idxTrn  = training(cvp); % Training set indices
                idxTest = test(cvp);    % Test set indices
                fprintf(['(' num2str(size(spike_Count_Train, 1)) ', ' num2str(size(spike_Count_Train, 2)) '), '])
                if strcmpi(Param.Classifier, 'lda')
                    
                    Mdl                                      = fitcdiscr(spike_Count_Train(idxTrn, :), Labels(idxTrn, :));
                    labels                                   = predict(Mdl, spike_Count_Test(idxTest, :));
                    accuracy_Matrix(time_Train_Ind, time_Test_Ind, Group_Cnt) = accuracy_Matrix(time_Train_Ind, time_Test_Ind, Group_Cnt) + 100*mean(labels==Labels(idxTest, :));
                    
                elseif strcmpi(Param.Classifier, 'nbc')
                    
                    Mdl                                      = fitcnb(spike_Count_Train(idxTrn, :), Labels(idxTrn, :));
                    labels                                   = predict(Mdl, spike_Count_Test(idxTest, :));
                    accuracy_Matrix(time_Train_Ind, time_Test_Ind, Group_Cnt) = accuracy_Matrix(time_Train_Ind, time_Test_Ind, Group_Cnt) + 100*mean(labels==Labels(idxTest, :));
                    
                elseif strcmpi(Param.Classifier, 'svm')
                    
                    Mdl                                      = fitcecoc(spike_Count_Train(idxTrn, :), Labels(idxTrn, :));
                    labels                                   = predict(Mdl, spike_Count_Test(idxTest, :));
                    accuracy_Matrix(time_Train_Ind, time_Test_Ind, Group_Cnt) = accuracy_Matrix(time_Train_Ind, time_Test_Ind, Group_Cnt) + 100*mean(labels==Labels(idxTest, :));
                    
                else
                    error('**** The classifier should be either *LDA* or *NBC* (Naive Bayes Classifier)')
                end
                
                
            end
            
            time_Test_Ind = time_Test_Ind + 1;
        end
        time_Train_Ind = time_Train_Ind + 1;
    end
    
    Group_Cnt = Group_Cnt +1;
end
save(['accuracy_Matrix_Generalization_' date '.mat'], 'accuracy_Matrix', 'Param')
%%
colormap;
figure(1)
for iGroup = 1 : size(accuracy_Matrix, 3)
    subplot(2,2,iGroup)
    imagesc(accuracy_Matrix(:, end:-1:1, iGroup)./Param.NumCrossVal);
    aX = gca;
    aX.TickDir = 'out';
    aX.Box    = 'off';
    aX.TickLength = 3*aX.TickLength;
    aX.XLabel.String = 'Training Window';
    aX.YLabel.String = 'Testing Window';
    aX.YTick  = 2 : 26 : size(accuracy_Matrix, 1);
    aX.XTick  = 1 : 26 : size(accuracy_Matrix, 1);
    aX.XTickLabel = linspace(0, Param.PostStim - Param.MovWind, length(aX.XTick));
    aX.YTickLabel = linspace(Param.PostStim - Param.MovWind, 0, length(aX.YTick));
    caxis([1/Param.NumUniqStim 90])
    
    h = colorbar;
    h.Box = 'off';
    h.TickDirection = 'out';
    h.TickLength = 3*h.TickLength;
    h.Label.String = 'Accuracy';
    
end

