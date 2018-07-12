clear
close all
clc

% load the NEV file and do some pre-processing.
data_Path = 'F:\CJ194\Data\';
data_FileName = 'CJ194_datafile028.nev';

stimulus_Path = 'F:\CJ194\Stimulus\';
stimulus_FileName = 'Paired_Stimulus_File_CJ194_0005.mat';

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
end

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

%% make a spike train for each selected channel

sTrain = zeros(length(select_Electrodes), ceil(max(spikes)));
for iElectrode = 1 : length(select_Electrodes)
    sTrain(iElectrode, round(spikes(dat.Data.Spikes.Electrode == select_Electrodes(iElectrode)))) = 1;
end

%% PSHT for a n*n pairing matrix for each selected channel individually
close all

SDF_binSize       = 20;  % ms
leadStimDuration  = 1000*stim.durationLeadStim;  % presentation time of leading stimulus
trailStimDuration = 1000*stim.durationTrailStim;  % presentation time of trailing stimulus
ISIDurartin       = 1000*stim.durationISI;  % ISI duration

winSize           = leadStimDuration + trailStimDuration + ISIDurartin +200;  % ms (PSTH length)

% there can be 4 Alignments: start of 1st stim, end of 1st stim, start of
% ISI, comments
other_Alignments = [0 winSize;
    -leadStimDuration winSize-leadStimDuration;
    -(leadStimDuration + trailStimDuration + ISIDurartin) winSize-(leadStimDuration + trailStimDuration + ISIDurartin)];

group_Trials = 100;  % group every "group_Trails" trials to see the effect of learning
line_Color   = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
FigureTab   = true;  % if ture it plots a figure with mutiple tabs, otheriwse, multiple figure

if FigureTab
    figure('units','normalized','outerposition',[0 0 1 1]);
    tab_group = uitabgroup; % tabgroup
end

for iElectrode = 1 : length(select_Electrodes)

    if FigureTab
        thistab = uitab(tab_group);  % build a tab
        axes('Parent', thistab); % somewhere to plot
    else
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    sdf = conv(ones(1, SDF_binSize), sTrain(iElectrode,:))*(1/(SDF_binSize/1000));

    iPanel = 1;
    for iLeadStim = 1 : length(stim_LeadStim)
        for iTrailStim = 1 : length(stim_TrailStim)

            subplot(6, 6, iPanel)
            this_Pair = find(stim_Train(1, :) == stim_LeadStim(iLeadStim) &  stim_Train(2, :) == stim_TrailStim(iTrailStim));

            % get some timing event for different alignments
            this_Stim_OnTime1  = round(stim_OnTime1(this_Pair));
            this_stim_OffTime1 = round(stim_OffTime1(this_Pair));
            this_stim_OffTime2 = round(stim_OffTime2(this_Pair));
            this_comments      = round(comments(this_Pair));

            % PSTH aligned to the start of first event
            this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_Stim_OnTime1  + other_Alignments(1, 1)-1]';
            resps1 = sdf(this_Epochs);
            % PSTH aligned to the end of first event
            this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime1 + other_Alignments(2, 1)-1]';
            resps2 = sdf(this_Epochs);
            % PSTH aligned to the start of ISI
            this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime2 + other_Alignments(3, 1)-1]';
            resps3 = sdf(this_Epochs);
            % PSTH aligned to the start of comments
            this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_comments      + other_Alignments(3, 1)-1]';
            resps4 = sdf(this_Epochs);

            plot(1:winSize, mean(resps1), 'color', line_Color(1, :)), hold on
            plot(1:winSize, mean(resps2), 'color', line_Color(2, :)), hold on
            plot(1:winSize, mean(resps3), 'color', line_Color(3, :)), hold on
            h = plot(1:winSize, mean(resps4), 'color', line_Color(4, :));
            h.Parent.Box = 'off';
            h.Parent.TickDir = 'out';
            if iPanel == 31
                h.Parent.XLabel.String = 'Time (ms)';
                h.Parent.YLabel.String = 'Firing Rate (spk/s)';
            end
            iPanel = iPanel + 1;
        end
    end

    if FigureTab
        thistab.Title = ['Chn ' num2str(iElectrode)];
    else
        suptitle = ['Chn ' num2str(iElectrode)];
    end
end

%%

if FigureTab
    figure('units','normalized','outerposition',[0 0 1 1]);
    tab_group = uitabgroup; % tabgroup
end
select_Alignments = 1;
line_Color2 = colormap('parula');
line_width  = 1;
for iElectrode = 1 : length(select_Electrodes)
    
    if FigureTab
        thistab = uitab(tab_group);  % build a tab
        axes('Parent', thistab); % somewhere to plot
    else
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    sdf = conv(ones(1, SDF_binSize), sTrain(iElectrode,:))*(1/(SDF_binSize/1000));
    
    iPanel = 1;
    for iLeadStim = 1 : length(stim_LeadStim)
        for iTrailStim = 1 : length(stim_TrailStim)
            
            subplot(6, 6, iPanel)
            this_Pair = find(stim_Train(1, :) == stim_LeadStim(iLeadStim) &  stim_Train(2, :) == stim_TrailStim(iTrailStim));
             
            switch select_Alignments
                case 1
                    % PSTH aligned to the start of first event
                    this_Stim_OnTime1  = round(stim_OnTime1(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_Stim_OnTime1  + other_Alignments(1, 1)-1]';
                    resps = sdf(this_Epochs);
                    %                     resps = zeros(length(this_Stim_OnTime1), winSize);
                    %                     for iStim = 1 : length(this_Stim_OnTime1)
                    %                         resps(iStim, :) = sdf(this_Stim_OnTime1(iStim) : this_Stim_OnTime1(iStim) + winSize-1);
                    %
                    %                     end
                    
                case 2
                    
                    % PSTH aligned to the end of first event

                    this_stim_OffTime1 = round(stim_OffTime1(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime1 + other_Alignments(2, 1)-1]';
                    resps = sdf(this_Epochs);
                case 3
                    
                    % PSTH aligned to the start of ISI
                                        this_stim_OffTime2 = round(stim_OffTime2(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime2 + other_Alignments(3, 1)-1]';
                    resps = sdf(this_Epochs);
                case 4
                    % PSTH aligned to the start of comments
                    this_comments      = round(comments(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_comments       + other_Alignments(3, 1)-1]';
                    resps = sdf(this_Epochs);
            end
            
            if  length(this_Pair) > stim.trialsOtherStim
                remaining_Trials = mod(length(this_Pair),group_Trials);
                color_Ind = 1;
                for iTrial = 1 : group_Trials: length(this_Pair)
                    if (iTrial + group_Trials-1)<=length(this_Pair)
                        h = plot(1:winSize, mean(resps(iTrial : iTrial + group_Trials-1, : )), 'color', line_Color2(color_Ind, :)); hold on
                    else
                        h = plot(1:winSize, mean(resps(iTrial : end, : )), 'color', line_Color2(color_Ind, :)); hold on
                    end
                    h.LineWidth = line_width;
                    color_Ind = color_Ind  + 15;
                end
            else
                h = plot(1:winSize, mean(resps), 'color', line_Color(select_Alignments, :));
                h.LineWidth = line_width;
                hold on
            end
            
            h.Parent.Box = 'off';
            h.Parent.TickDir = 'out';
            if iPanel == 31
                h.Parent.XLabel.String = 'Time (ms)';
                h.Parent.YLabel.String = 'Firing Rate (spk/s)';
            end
            iPanel = iPanel + 1;
        end
    end
    
    
    if FigureTab
        thistab.Title = ['Chn ' num2str(iElectrode)];
    else
        suptitle = ['Chn ' num2str(iElectrode)];
    end
    
end

%% here we do some normalization with the hope to slove latency issue, Nic's suggestion

sdf_Normalized = zeros(size(sdf));
for iElectrode = 1 : length(select_Electrodes)
    
    sdf_Normalized = sdf_Normalized + zscore(conv(ones(1, SDF_binSize), sTrain(iElectrode,:))*(1/(SDF_binSize/1000)));
    
end
sdf_Normalized = sdf_Normalized./length(select_Electrodes);



figure('units','normalized','outerposition',[0 0 1 1]);
select_Alignments = 2;
line_Color2 = colormap('parula');
line_width  = 1;

iPanel = 1;
for iLeadStim = 1 : length(stim_LeadStim)
    for iTrailStim = 1 : length(stim_TrailStim)
        
        subplot(6, 6, iPanel)
        this_Pair = find(stim_Train(1, :) == stim_LeadStim(iLeadStim) &  stim_Train(2, :) == stim_TrailStim(iTrailStim));
        
        % get some timing event for different alignments
        this_Stim_OnTime1  = round(stim_OnTime1(this_Pair));
        this_stim_OffTime1 = round(stim_OffTime1(this_Pair));
        this_stim_OffTime2 = round(stim_OffTime2(this_Pair));
        this_comments      = round(comments(this_Pair));
        
        switch select_Alignments
            case 1
                % PSTH aligned to the start of first event
                this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_Stim_OnTime1  + other_Alignments(1, 1)-1]';
                resps = sdf_Normalized(this_Epochs);
            case 2
                
                % PSTH aligned to the end of first event
                this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime1 + other_Alignments(2, 1)-1]';
                resps = sdf_Normalized(this_Epochs);
            case 3
                
                % PSTH aligned to the start of ISI
                this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime2 + other_Alignments(3, 1)-1]';
                resps = sdf_Normalized(this_Epochs);
            case 4
                % PSTH aligned to the start of comments
                this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_comments       + other_Alignments(3, 1)-1]';
                resps = sdf_Normalized(this_Epochs);
        end
        
        if  length(this_Pair) > stim.trialsOtherStim
            remaining_Trials = mod(length(this_Pair),group_Trials);
            color_Ind = 1;
            for iTrial = 1 : group_Trials: length(this_Pair)
                if (iTrial + group_Trials-1)<=length(this_Pair)
                    h = plot(1:winSize, mean(resps(iTrial : iTrial + group_Trials-1, : )), 'color', line_Color2(color_Ind, :)); hold on
                else
                    h = plot(1:winSize, mean(resps(iTrial : end, : )), 'color', line_Color2(color_Ind, :)); hold on
                end
                h.LineWidth = line_width;
                color_Ind = color_Ind  + 15;
            end
        else
            h = plot(1:winSize, mean(resps), 'color', line_Color(select_Alignments, :));
            h.LineWidth = line_width;
            hold on
        end
        
        h.Parent.Box = 'off';
        h.Parent.TickDir = 'out';
        if iPanel == 31
            h.Parent.XLabel.String = 'Time (ms)';
            h.Parent.YLabel.String = 'Firing Rate (spk/s)';
        end
        iPanel = iPanel + 1;
    end
end


