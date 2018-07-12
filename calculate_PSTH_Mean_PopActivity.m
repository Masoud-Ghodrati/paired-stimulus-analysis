clear
close all
clc

% load the NEV file and do some pre-processing.
% data_Path = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Marmoset data\CJ194\';
% data_FileName = 'CJ194_datafile025.nev';
%
% stimulus_Path = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Marmoset data\CJ194 Stimulus Files\';
% stimulus_FileName = 'Paired_Stimulus_File_CJ194_0001.mat';
data_Path = 'F:\CJ194\Data\';
data_FileName = 'CJ194_datafile030.nev';

stimulus_Path = 'F:\CJ194\Stimulus\';
stimulus_FileName = 'Paired_Stimulus_File_CJ194_0007.mat';

load([stimulus_Path stimulus_FileName])

if exist([data_Path data_FileName(1:end-3) 'mat'], 'file')
    load([data_Path data_FileName(1:end-3) 'mat'])
else
    openNEV([data_Path data_FileName], 'read', 'nosave');
    NEV.Data.Spikes.Waveform = [];
    save([data_Path data_FileName(1:end-3) 'mat'], 'NEV')
end

if exist('NEVdata', 'var')
    NEV   = NEVdata;
    clear   NEVdata;
end

%% Extract some event information and timing
dat       = cbmex_Parse_data(NEV);
clear   NEV;
tRes      = dat.MetaTags.TimeRes;  % sampling resolution
spikes    = double(dat.Data.Spikes.TimeStamp)/tRes*1000;  % spike times (ms)

% Digital Timings
RawDIO        = dat.Data.SerialDigitalIO.UnparsedData;  % DIO tags
RawTimes      = double(dat.Data.SerialDigitalIO.TimeStamp)/tRes*1000;  % DIO digital time(ms)
DIO           = mod(RawDIO, 128);   % digital line without photodiode
stim_OnTime1  = RawTimes(DIO == 2);  % stim 1 onset
stim_OffTime1 = RawTimes(DIO == 3);  % stim 1 offset
stim_OnTime2  = RawTimes(DIO == 4);  % stim 2 onset
stim_OffTime2 = RawTimes(DIO == 5);  % stim 1 offset


% Channels information
electrodes        = unique(dat.Data.Spikes.Electrode);  % electrode numbers
if strcmp(data_FileName, 'CJ194_datafile025.nev')
    [stim_OnTime1, stim_OffTime1, stim_OffTime2, comments] = correct_Timing_CJ194_datafile025(stim_OnTime1, stim_OffTime1, stim_OffTime2, dat, stim);
    select_Electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,14,19,21,22,26,27,29,31,32,37,40,41,42,44,46,50,51,52,53,54,55,56,57,58,62,63,65,66,67,73,75,76,81,83,84,85,86,87,94,95]; % 25
elseif strcmp(data_FileName, 'CJ194_datafile026.nev')
    [stim_OnTime1, stim_OffTime1, stim_OffTime2, comments] = correct_Timing_CJ194_datafile026(stim_OnTime1, stim_OffTime1, stim_OffTime2, dat, stim);
    select_Electrodes = [1:14 16:19 21 22 26 27 29 31 32 37 40:42 44 46 47 50:58 61:70 73 75 76 81 83:88 91 93:96 ]; % 26
elseif strcmp(data_FileName, 'CJ194_datafile028.nev')
    [stim_OnTime1, stim_OffTime1, stim_OffTime2, comments] = correct_Timing_CJ194_datafile028(stim_OnTime1, stim_OffTime1, stim_OffTime2, dat, stim);
    select_Electrodes = [1:12 17 19 21 23 26 27 29 32 32 37 40 41 42 44 46 50 51:57 66 73 75 76 83 85 86 87]; % 28
elseif strcmp(data_FileName, 'CJ194_datafile030.nev')
    [stim_OnTime1, stim_OffTime1, stim_OffTime2, comments] = correct_Timing_CJ194_datafile030(stim_OnTime1, stim_OffTime1, stim_OffTime2, dat, stim);
    select_Electrodes = [1:12 14:17 19 21 22 26 27 29 31 32 37 40:42 44 46 50:57 63 64 66 73 75 76 83 84:87 93:96];
else
    
    cStruct   = dat.Data.Comments;  % comments
    comments1 = double([cStruct.TimeStamp])/tRes*1000;  % comment times (ms)
    % txt = reshape(NEV.Data.Comments.Text,[],92);
    % comment_txt        = reshape(cStruct.Comments,[],92);
    % [match, noMatch]   = regexp(cellstr(comment_txt(:, 1:22)),'\d','match','forceCellOutput');
    % trial_NumCellArray = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
    % trial_NumArray     = cellfun(@str2num, trial_NumCellArray, 'UniformOutput', false);
    % % find(diff(cell2mat(trial_NumArray))>1)+1;
    % comments           = comments1(~[1 0 diff(cell2mat(trial_NumArray))'>1]);
    
    comments      = comments1;
    select_Electrodes = 1:96;
end
% Photodiode
PDTimes = double(dat.Data.Spikes.TimeStamp(dat.Data.Spikes.Electrode == 129))/tRes*1000;
PDTimes = PDTimes(PDTimes > comments(1));

%% extract stimulus information
if strcmpi(stim.textureType, 'texture')
    stim_LeadStim  = stim.TextFamilies(1:length(stim.TextFamilies)/2);  % leading stimulus names/indexes
    stim_TrailStim = stim.TextFamilies(1+(length(stim.TextFamilies)/2):end);  % trailing stimulus names/indexes
else
    stim_LeadStim  = 1:length(stim.oriList)/2;  % leading stimulus names/indexes
    stim_TrailStim = (1+(length(stim.oriList)/2)):length(stim.oriList);  % trailing stimulus names/indexes
end
stim_Train     = stim.allStimTrain;  % stimulus train. This should be a matrix of 3*n. 1st row: leading stim name/ind, 2nd trailing stim name/ind, last sample number
stim_Images    = stim.allStimFile;  % presented image file

if (size([stim_OnTime1; stim_OffTime1; stim_OffTime2; comments], 2) - length(stim_Train)) ~= 0
    
    error('DIO/Comments length does not match with stim length')
    
end
%% make a spike train for each selected channel

sTrain = zeros(length(select_Electrodes), ceil(max(spikes)));
for iElectrode = 1 : length(select_Electrodes)
    sTrain(iElectrode, round(spikes(dat.Data.Spikes.Electrode == select_Electrodes(iElectrode)))) = 1;
end

%% PSHT for a n*n pairing matrix for each selected channel individually
close all

SDF_binSize       = 35;  % ms
leadStimDuration  = 1000*stim.durationLeadStim;  % presentation time of leading stimulus
trailStimDuration = 1000*stim.durationTrailStim;  % presentation time of trailing stimulus
ISIDurartin       = 1000*stim.durationISI;  % ISI duration
pre_Stim          = 100;  % time before stimulus onset
post_Stim         = 200;  % time after stimulus onset
winSize           = leadStimDuration + trailStimDuration + ISIDurartin + pre_Stim + post_Stim;  % ms (PSTH length)

% there can be 4 Alignments: start of 1st stim, end of 1st stim, start of
% ISI, comments
other_Alignments  = [-pre_Stim; -(leadStimDuration + pre_Stim); -(leadStimDuration + trailStimDuration + ISIDurartin + pre_Stim)];
group_Trials      = 100;  % group every "group_Trails" trials to see the effect of learning
FigureTab         = true;  % if ture it plots a figure with mutiple tabs, otheriwse, multiple figure
select_Alignments = 1;

%% plot predicted and unpredicted responses for single channels
figure(1)
line_Color = colormap('parula');
line_width = 1;
predicted_Resp             = cell(1, 6);
predicted_Resp_SingleTrial = cell(1, 6);
unpredicted_resp           = cell(1, 6);
if FigureTab
    %     figure('units','normalized','outerposition',[0 0 1 1]);
    tab_group = uitabgroup; % tabgroup
end
for iElectrode = 1 : length(select_Electrodes)
    
    if FigureTab
        thistab = uitab(tab_group);  % build a tab
        axes('Parent', thistab);  % somewhere to plot
    else
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    sdf = conv(ones(1, SDF_binSize), sTrain(iElectrode,:))*(1/(SDF_binSize/1000));
    
    for iTrailStim = 1 : length(stim_TrailStim)
        temp_Resp = [];
        subplot(1, 6, iTrailStim)
        
        for iLeadStim = 1 : length(stim_LeadStim)
            
            this_Pair = find(stim_Train(1, :) == stim_LeadStim(iLeadStim) &  stim_Train(2, :) == stim_TrailStim(iTrailStim));
            switch select_Alignments
                case 1
                    % PSTH aligned to the start of first event
                    this_Stim_OnTime1  = round(stim_OnTime1(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_Stim_OnTime1  + other_Alignments(1)-1]';
                    resps = sdf(this_Epochs);
                case 2
                    % PSTH aligned to the end of first event
                    this_stim_OffTime1 = round(stim_OffTime1(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime1 + other_Alignments(2)-1]';
                    resps = sdf(this_Epochs);
                case 3
                    % PSTH aligned to the start of ISI
                    this_stim_OffTime2 = round(stim_OffTime2(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime2 + other_Alignments(3)-1]';
                    resps = sdf(this_Epochs);
                case 4
                    % PSTH aligned to the start of comments
                    this_comments      = round(comments(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_comments      + other_Alignments(3)-1]';
                    resps = sdf(this_Epochs);
            end
            
            if iTrailStim == iLeadStim
                predicted_Resp{1, iTrailStim}(iElectrode, :) =  mean(resps);
                predicted_Resp_SingleTrial{1, iTrailStim}(:, :, iElectrode) = resps;
            else
                temp_Resp = [temp_Resp; resps];
                
            end
            
        end
        unpredicted_resp{1, iTrailStim}(iElectrode, :) = mean(temp_Resp);
        
        plot(predicted_Resp{1, iTrailStim}(iElectrode, :), 'r'), hold on
        h = plot(unpredicted_resp{1, iTrailStim}(iElectrode, :), 'b');
        
        h.Parent.Box = 'off';
        h.Parent.TickDir = 'out';
        if iTrailStim == 1
            h.Parent.XLabel.String = 'Time (ms)';
            h.Parent.YLabel.String = 'Firing Rate (spk/s)';
        end
    end
    
    legend('Predicted', 'Nonpredicted')
    legend boxoff
    if FigureTab
        thistab.Title = ['Chn ' num2str(iElectrode)];
    else
        suptitle = ['Chn ' num2str(iElectrode)];
    end
    
end

%% plot population average for predicted and unpredicted responses
figure(2)
for iTrailStim = 1 : length(stim_TrailStim)
    
    subplot(1, 6, iTrailStim)
    plot(mean(predicted_Resp{1, iTrailStim}), 'r'), hold on
    h = plot(mean(unpredicted_resp{1, iTrailStim}), 'b');
    
    h.Parent.Box = 'off';
    h.Parent.TickDir = 'out';
    if iTrailStim == 1
        h.Parent.XLabel.String = 'Time (ms)';
        h.Parent.YLabel.String = 'Firing Rate (spk/s)';
    end
end
legend('Predicted', 'Nonpredicted')
legend boxoff

%% plot 6 pairs every group_Trials (100) trials, averaged over all channels
figure(3)
sub_Ind = 1;
for iTrial = 1 : group_Trials: size(predicted_Resp_SingleTrial{1}, 1)
    
    color_Ind = 1;
    subplot(1,round(size(predicted_Resp_SingleTrial{1}, 1)/group_Trials), sub_Ind)
    
    for iPair = 1 : 6
        if (iTrial + group_Trials-1)<=length(this_Pair)
            this_Mean = mean(mean(predicted_Resp_SingleTrial{iPair}(iTrial : iTrial + group_Trials-1, :, : )), 3);
            h = plot(1:winSize, this_Mean, 'color', line_Color(color_Ind, :)); hold on
        else
            this_Mean = mean(mean(predicted_Resp_SingleTrial{iPair}(iTrial : end, :, : )), 3);
            h = plot(1:winSize, this_Mean, 'color', line_Color(color_Ind, :)); hold on
        end
        h.LineWidth = line_width;
        h.Parent.Box = 'off';
        h.Parent.TickDir = 'out';
        color_Ind = color_Ind  + 10;
    end
    
    if sub_Ind == 1
        h.Parent.XLabel.String = 'Time (ms)';
        h.Parent.YLabel.String = 'Firing Rate (spk/s)';
        
    end
    sub_Ind = sub_Ind + 1;
    
end
legend('Pair 1','Pair 2','Pair 3','Pair 4','Pair 5','Pair 6')
legend boxoff


%% plot 6 pairs every group_Trials (100) trials, for every single channels
figure(4)
if FigureTab
    %     figure('units','normalized','outerposition',[0 0 1 1]);
    tab_group = uitabgroup; % tabgroup
end
for iElectrode = 1 : length(select_Electrodes)
    
    if FigureTab
        thistab = uitab(tab_group);  % build a tab
        axes('Parent', thistab); % somewhere to plot
    else
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    sub_Ind = 1;
    for iTrial = 1 : group_Trials: size(predicted_Resp_SingleTrial{1}, 1)
        
        color_Ind = 1;
        subplot(1,round(size(predicted_Resp_SingleTrial{1}, 1)/group_Trials),sub_Ind)
        
        for iPair = 1 : 6
            if (iTrial + group_Trials-1)<=length(this_Pair)
                this_Mean = mean(mean(predicted_Resp_SingleTrial{iPair}(iTrial : iTrial + group_Trials-1, :, iElectrode )), 3);
                h = plot(1:winSize, this_Mean, 'color', line_Color(color_Ind, :)); hold on
            else
                this_Mean = mean(mean(predicted_Resp_SingleTrial{iPair}(iTrial : end, :, iElectrode )), 3);
                h = plot(1:winSize, this_Mean, 'color', line_Color(color_Ind, :)); hold on
            end
            h.LineWidth = line_width;
            h.Parent.Box = 'off';
            h.Parent.TickDir = 'out';
            color_Ind = color_Ind  + 10;
        end
        
        if sub_Ind == 1
            h.Parent.XLabel.String = 'Time (ms)';
            h.Parent.YLabel.String = 'Firing Rate (spk/s)';
            
        end
        sub_Ind = sub_Ind + 1;
        
    end
    legend('Pair 1','Pair 2','Pair 3','Pair 4','Pair 5','Pair 6')
    legend boxoff
    if FigureTab
        thistab.Title = ['Chn ' num2str(select_Electrodes(iElectrode))];
    else
        suptitle = ['Chn ' num2str(select_Electrodes(iElectrode))];
    end
end

%% show images of all 6 pairs
figure(5)

for iPair = 1 : 6
    color_Ind = 1;
    subplot(1, 6, iPair)
    for iTrial = 1 : group_Trials: size(predicted_Resp_SingleTrial{1}, 1)
        
        if (iTrial + group_Trials-1)<=length(this_Pair)
            this_Mean = mean(mean(predicted_Resp_SingleTrial{iPair}(iTrial : iTrial + group_Trials-1, :, : )), 3);
            h = plot(1:winSize, this_Mean, 'color', line_Color(color_Ind, :)); hold on
        else
            this_Mean = mean(mean(predicted_Resp_SingleTrial{iPair}(iTrial : end, :, : )), 3);
            h = plot(1:winSize, this_Mean, 'color', line_Color(color_Ind, :)); hold on
        end
        h.LineWidth = line_width;
        color_Ind = color_Ind  + 15;
        
        h.Parent.Box = 'off';
        h.Parent.TickDir = 'out';
        
    end
    if iPair == 1
        h.Parent.XLabel.String = 'Time (ms)';
        h.Parent.YLabel.String = 'Firing Rate (spk/s)';
        
    end
end
legend('0-100','101-200','201-300','301-400')
legend boxoff