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
data_FileName = 'CJ194_datafile026.nev';

stimulus_Path = 'F:\CJ194\Stimulus\';
stimulus_FileName = 'Paired_Stimulus_File_CJ194_0002.mat';

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
dat = cbmex_Parse_data(NEV);
clear   NEV;
tRes = dat.MetaTags.TimeRes;  % sampling resolution
cStruct = dat.Data.Comments;  % comments
comments1 = double([cStruct.TimeStamp])/tRes*1000;  % comment times (ms)

spikes = double(dat.Data.Spikes.TimeStamp)/tRes*1000;  % spike times (ms)


% txt = reshape(NEV.Data.Comments.Text,[],92);
% comment_txt        = reshape(cStruct.Comments,[],92);
% [match, noMatch]   = regexp(cellstr(comment_txt(:, 1:22)),'\d','match','forceCellOutput');
% trial_NumCellArray = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
% trial_NumArray     = cellfun(@str2num, trial_NumCellArray, 'UniformOutput', false);
% % find(diff(cell2mat(trial_NumArray))>1)+1;
% comments           = comments1(~[1 0 diff(cell2mat(trial_NumArray))'>1]);

comments           = comments1;
% Digital Timings
RawDIO        = dat.Data.SerialDigitalIO.UnparsedData;  % DIO tags
RawTimes      = double(dat.Data.SerialDigitalIO.TimeStamp)/tRes*1000;  % DIO digital time(ms)
DIO           = mod(RawDIO, 128);   % digital line without photodiode
stim_OnTime1  = RawTimes(DIO == 2);  % stim 1 onset
stim_OffTime1 = RawTimes(DIO == 3);  % stim 1 offset
stim_OnTime2  = RawTimes(DIO == 4);  % stim 2 onset
stim_OffTime2 = RawTimes(DIO == 5);  % stim 1 offset

% Photodiode
PDTimes = double(dat.Data.Spikes.TimeStamp(dat.Data.Spikes.Electrode == 129))/tRes*1000;
PDTimes = PDTimes(PDTimes > comments(1));

% Channels information
electrodes = unique(dat.Data.Spikes.Electrode);  % electrode numbers
% select_Electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,14,19,21,22,26,27,29,31,32,37,40,41,42,44,46,50,51,52,53,54,55,56,57,58,62,63,65,66,67,73,75,76,81,83,84,85,86,87,94,95]; % 25
select_Electrodes = [1:14 16:19 21 22 26 27 29 31 32 37 40:42 44 46 47 50:58 61:70 73 75 76 81 83:88 91 93:96 ]; % 26
% select_Electrodes = [1:12 17 19 21 23 26 27 29 32 32 37 40 41 42 44 46 50 51:57 66 73 75 76 83 85 86 87]; % 28
%% extract stimulus information
if strcmpi(stim.textureType, 'texture')
    stim_LeadStim  = stim.TextFamilies(1:length(stim.TextFamilies)/2);  % leading stimulus names/indexes
    stim_TrailStim = stim.TextFamilies(1+(length(stim.TextFamilies)/2):end);  % trailing stimulus names/indexes
else
    stim_LeadStim  = 1:length(stim.oriList)/2;  % leading stimulus names/indexes
    stim_TrailStim = (1+(length(stim.oriList)/2)):length(stim.oriList);  % trailing stimulus names/indexes
    %     stim_LeadStim = stim_TrailStim;
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

SDF_binSize       = 5;  % ms
leadStimDuration  = 1000*stim.durationLeadStim;  % presentation time of leading stimulus
trailStimDuration = 1000*stim.durationTrailStim;  % presentation time of trailing stimulus
ISIDurartin       = 1000*stim.durationISI;  % ISI duration

winSize           = leadStimDuration + trailStimDuration + ISIDurartin + 100;  % ms (PSTH length)

% there can be 4 Alignments: start of 1st stim, end of 1st stim, start of
% ISI, comments
other_Alignments = [0 winSize;
    -leadStimDuration winSize-leadStimDuration;
    -(leadStimDuration + trailStimDuration + ISIDurartin) winSize-(leadStimDuration + trailStimDuration + ISIDurartin)];

group_Trials    = 100;  % group every "group_Trails" trials to see the effect of learning
line_Color      = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
baseline_Window = round(1000*[stim.durationITI(1)-0.1 stim.durationITI(1)]);
FigureTab       = true;  % if ture it plots a figure with mutiple tabs, otheriwse, multiple figure
GeoMean         = true;

%% find the time window of analysis

if FigureTab
    tab_group = uitabgroup; % tabgroup
end
select_Alignments = 1;
for iElectrode = 1 : length(select_Electrodes)
    
    if FigureTab
        thistab = uitab(tab_group);  % build a tab
        axes('Parent', thistab); % somewhere to plot
    else
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    sdf = conv(ones(1, SDF_binSize), sTrain(iElectrode,:))*(1/(SDF_binSize/1000));
    all_Resp = [];
    for iTrailStim = 1 : length(stim_TrailStim)
        
        for iLeadStim = 1 : length(stim_LeadStim)
            
            this_Pair = find(stim_Train(1, :) == stim_LeadStim(iLeadStim) &  stim_Train(2, :) == stim_TrailStim(iTrailStim));
            
            switch select_Alignments
                case 1
                    % PSTH aligned to the start of first event
                    this_Stim_OnTime1  = round(stim_OnTime1(this_Pair));
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_Stim_OnTime1  + other_Alignments(1, 1)-1]';
                    resps = sdf(this_Epochs);
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
            
            all_Resp = [all_Resp; resps];
            
        end
        
        
    end
    
    plot(mean(all_Resp)), hold on
    mean_Resp              = mean(all_Resp);
    time_Wind_1            = [1 100];
    time_Wind_2            = [150 winSize];
    [max_Val_1, max_Ind_1] = max(mean_Resp(time_Wind_1(1):time_Wind_1(2)));
    [max_Val_2, max_Ind_2] = max(mean_Resp(time_Wind_2(1):time_Wind_2(2)));
    
    plot(max_Ind_1,max_Val_1, 'ro');
    hold on
    h2 = plot(max_Ind_2 + time_Wind_2(1), max_Val_2, 'ro');
    h2.Parent.Box = 'off';
    h2.Parent.TickDir = 'out';
    h2.Parent.XLabel.String = 'Time (ms)';
    h2.Parent.YLabel.String = 'Firing Rate (spk/s)';
    if FigureTab
        thistab.Title = ['Chn ' num2str(iElectrode)];
    else
        suptitle = ['Chn ' num2str(iElectrode)];
    end
    
    cell_Analysis_Wind(:, iElectrode) = [max_Ind_1; max_Ind_2 + time_Wind_2(1)];
end


%% calculate some spike count stat (e.g., FR change ration over time)
for iElectrode = 1 : length(select_Electrodes)
    
    sdf      = conv(ones(1, SDF_binSize), sTrain(iElectrode,:))*(1/(SDF_binSize/1000));
    pair_Ind = 1;
    for iTrailStim = 1 : length(stim_TrailStim)
        
        for iLeadStim = 1 : length(stim_LeadStim)
            
            if iTrailStim == iLeadStim
                
                this_Pair = find(stim_Train(1, :) == stim_LeadStim(iLeadStim) &  stim_Train(2, :) == stim_TrailStim(iTrailStim));
                
                switch select_Alignments
                    case 1
                        % PSTH aligned to the start of first event
                        this_Stim_OnTime1  = round(stim_OnTime1(this_Pair));
                        this_Epochs        = repmat(1:winSize, [length(this_Pair) 1]) + [this_Stim_OnTime1  + other_Alignments(1, 1)-1]';
                        resps              = sdf(this_Epochs);
                    case 2
                        % PSTH aligned to the end of first event
                        this_stim_OffTime1 = round(stim_OffTime1(this_Pair));
                        this_Epochs        = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime1 + other_Alignments(2, 1)-1]';
                        resps              = sdf(this_Epochs);
                    case 3
                        % PSTH aligned to the start of ISI
                        this_stim_OffTime2 = round(stim_OffTime2(this_Pair));
                        this_Epochs        = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime2 + other_Alignments(3, 1)-1]';
                        resps              = sdf(this_Epochs);
                    case 4
                        % PSTH aligned to the start of comments
                        this_comments      = round(comments(this_Pair));
                        this_Epochs        = repmat(1:winSize, [length(this_Pair) 1]) + [this_comments       + other_Alignments(3, 1)-1]';
                        resps              = sdf(this_Epochs);
                end
                
                remaining_Trials = mod(length(this_Pair), group_Trials);
                
                groupd_Ind              = 1;
                peak_Time_1             = [];
                peak_Time_2             = [];
                peak_rate               = [];
                mean_Group_FR           = [];
                analysis_Wind_HalfWidth = 15;
                for iTrial = 1 : group_Trials : length(this_Pair)
                    if (iTrial + group_Trials-1) <= length(this_Pair)
                        
                        mean_Group_FR(groupd_Ind, :)  = mean(resps(iTrial : iTrial + group_Trials-1, :));
                        peak_Time_1(groupd_Ind)       = mean(mean(resps(iTrial : iTrial + group_Trials-1, cell_Analysis_Wind(1, iElectrode)-analysis_Wind_HalfWidth : cell_Analysis_Wind(1, iElectrode)+analysis_Wind_HalfWidth ), 2));
                        peak_Time_2(groupd_Ind)       = mean(mean(resps(iTrial : iTrial + group_Trials-1, cell_Analysis_Wind(2, iElectrode)-analysis_Wind_HalfWidth : cell_Analysis_Wind(2, iElectrode)+analysis_Wind_HalfWidth ), 2));
                        peak_rate(groupd_Ind)         = peak_Time_1(groupd_Ind)/peak_Time_2(groupd_Ind);
                        
                    else
                        mean_Group_FR(groupd_Ind, :) = mean(resps(iTrial : end, :));
                        peak_Time_1(groupd_Ind)      = mean(mean(resps(iTrial : end, cell_Analysis_Wind(1, iElectrode)-analysis_Wind_HalfWidth : cell_Analysis_Wind(1, iElectrode)+analysis_Wind_HalfWidth ), 2));
                        peak_Time_2(groupd_Ind)      = mean(mean(resps(iTrial : end, cell_Analysis_Wind(2, iElectrode)-analysis_Wind_HalfWidth : cell_Analysis_Wind(2, iElectrode)+analysis_Wind_HalfWidth ), 2));
                        peak_rate(groupd_Ind)        = peak_Time_1(groupd_Ind)/peak_Time_2(groupd_Ind);
                        mnn
                    end
                    
                    groupd_Ind = groupd_Ind  + 1;
                    
                end
                
                
                fr_Stat_Matrix{pair_Ind}(:,:, iElectrode) = [peak_Time_1; peak_Time_2; peak_rate];
                fr_mean_Pop{pair_Ind}(:,:, iElectrode)    = mean_Group_FR;
                pair_Ind = pair_Ind + 1;
                
            end
            
        end
        
    end
    
end
%% plot the results
close all

figure(1)
line_Color   = colormap('parula');
select_Param = 1;
x_Scaling    = 10;
for iPair = 1 : length(fr_Stat_Matrix)
    subplot(1,6, iPair)
    color_Ind    = 1;
    for iElectrode = 1 : length(select_Electrodes)
        
        xAXIS_Rand_Pos     = x_Scaling*randn;
        h1                 = plot([group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, fr_Stat_Matrix{iPair}(select_Param, : , iElectrode));
        h1.Color           = 0.7*[1 1 1];
        h1.LineWidth       = 0.5;
        hold on
        h2                 = scatter([group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, fr_Stat_Matrix{iPair}(select_Param, : , iElectrode), 20);
        h2.MarkerEdgeColor = 'none';
        h2.MarkerFaceColor = line_Color(color_Ind, :);
        h2.MarkerFaceAlpha = 0.5;
        
        %         color_Ind          = color_Ind + 1;
    end
    this_Pop_Data = fr_Stat_Matrix{iPair}(select_Param, : ,:);
    this_Pop_Data = squeeze(this_Pop_Data)';
    this_Pop_Data = this_Pop_Data(sum(isinf(this_Pop_Data),2) < 1, :);
    if GeoMean
        h2 = errorbar([group_Trials : group_Trials : length(this_Pair)], geomean(this_Pop_Data), std(this_Pop_Data)./sqrt(size(this_Pop_Data, 1)), 'or');
    else
        h2 = errorbar([group_Trials : group_Trials : length(this_Pair)], nanmean(this_Pop_Data), std(this_Pop_Data)./sqrt(size(this_Pop_Data, 1)), 'or');
    end
    h2.CapSize   = 0;
    h2.LineWidth = 2;
    h2.Color     = [1 0 0];
    aX          = gca;
    aX.TickDir  = 'out';
    aX.XTick    = group_Trials : group_Trials : length(this_Pair);
    aX.XLim     = [50 450];
    aX.Box      = 'off';
    
    if iPair == 1
        
        aX.XLabel.String = '# Trials';
        aX.YLabel.String = 'Firing Rate';
        
    end
end

figure(2)
select_Param = 2;
for iPair = 1 : length(fr_Stat_Matrix)
    subplot(1,6, iPair)
    color_Ind    = 1;
    for iElectrode = 1 : length(select_Electrodes)
        
        xAXIS_Rand_Pos     = x_Scaling*randn;
        h1                 = plot([group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, fr_Stat_Matrix{iPair}(select_Param, : , iElectrode));
        h1.Color           = 0.7*[1 1 1];
        h1.LineWidth       = 0.5;
        hold on
        h2                 = scatter([group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, fr_Stat_Matrix{iPair}(select_Param, : , iElectrode), 20);
        h2.MarkerEdgeColor = 'none';
        h2.MarkerFaceColor = line_Color(color_Ind, :);
        h2.MarkerFaceAlpha = 0.5;
        
        %         color_Ind          = color_Ind + 1;
    end
    this_Pop_Data = fr_Stat_Matrix{iPair}(select_Param, : ,:);
    this_Pop_Data = squeeze(this_Pop_Data)';
    this_Pop_Data = this_Pop_Data(sum(isinf(this_Pop_Data),2) < 1, :);
    if GeoMean
        h2 = errorbar([group_Trials : group_Trials : length(this_Pair)], geomean(this_Pop_Data), std(this_Pop_Data)./sqrt(size(this_Pop_Data, 1)), 'or');
    else
        h2 = errorbar([group_Trials : group_Trials : length(this_Pair)], nanmean(this_Pop_Data), std(this_Pop_Data)./sqrt(size(this_Pop_Data, 1)), 'or');
    end
    h2.CapSize   = 0;
    h2.LineWidth = 2;
    h2.Color     = [1 0 0];
    aX          = gca;
    aX.TickDir  = 'out';
    aX.XTick    = group_Trials : group_Trials : length(this_Pair);
    aX.XLim     = [50 450];
    aX.Box      = 'off';
    
    if iPair == 1
        
        aX.XLabel.String = '# Trials';
        aX.YLabel.String = 'Firing Rate';
        
    end
end


figure(3)
select_Param = 3;
for iPair = 1 : length(fr_Stat_Matrix)
    subplot(1,6, iPair)
    color_Ind    = 1;
    for iElectrode = 1 : length(select_Electrodes)
        
        xAXIS_Rand_Pos     = x_Scaling*randn;
        h1                 = plot([group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, fr_Stat_Matrix{iPair}(select_Param, : , iElectrode));
        h1.Color           = 0.7*[1 1 1];
        h1.LineWidth       = 0.5;
        hold on
        h2                 = scatter([group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, fr_Stat_Matrix{iPair}(select_Param, : , iElectrode), 20);
        h2.MarkerEdgeColor = 'none';
        h2.MarkerFaceColor = line_Color(color_Ind, :);
        h2.MarkerFaceAlpha = 0.5;
        
        %         color_Ind          = color_Ind + 1;
    end
    this_Pop_Data = fr_Stat_Matrix{iPair}(select_Param, : ,:);
    this_Pop_Data = squeeze(this_Pop_Data)';
    this_Pop_Data = this_Pop_Data(sum(isinf(this_Pop_Data),2) < 1, :);
    if GeoMean
        h2 = errorbar([group_Trials : group_Trials : length(this_Pair)], geomean(this_Pop_Data), nanstd(this_Pop_Data)./sqrt(size(this_Pop_Data, 1)), 'or');
    else
        h2 = errorbar([group_Trials : group_Trials : length(this_Pair)], nanmedian(this_Pop_Data), nanstd(this_Pop_Data)./sqrt(size(this_Pop_Data, 1)), 'or');
    end
    h2.CapSize   = 0;
    h2.LineWidth = 2;
    h2.Color     = [1 0 0];
    aX          = gca;
    aX.TickDir  = 'out';
    aX.XTick    = group_Trials : group_Trials : length(this_Pair);
    aX.XLim     = [50 450];
    aX.Box      = 'off';
    
    if iPair == 1
        
        aX.XLabel.String = '# Trials';
        aX.YLabel.String = 'Ratio (Peak 1 / Peak 2)';
        
    end
end

figure(4)
select_Param = 1;
YLIM         = [-1 1];
for iPair = 1 : length(fr_Stat_Matrix)
    subplot(1,6, iPair)
    color_Ind    = 1;
    for iElectrode = 1 : length(select_Electrodes)
        
        xAXIS_Rand_Pos     = x_Scaling*randn;
        
        this_Cell_Data     = fr_Stat_Matrix{iPair}(select_Param, : , iElectrode);
        this_Cell_FR_Change= (this_Cell_Data(2:end) - this_Cell_Data(1))./this_Cell_Data(1);
        h1                 = plot([2*group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, this_Cell_FR_Change);
        h1.Color           = 0.7*[1 1 1];
        h1.LineWidth       = 0.5;
        hold on
        h2                 = scatter([2*group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, this_Cell_FR_Change, 20);
        h2.MarkerEdgeColor = 'none';
        h2.MarkerFaceColor = line_Color(color_Ind, :);
        h2.MarkerFaceAlpha = 0.5;
        
        %         color_Ind          = color_Ind + 1;
    end
    this_Pop_Data   = fr_Stat_Matrix{iPair}(select_Param, : ,:);
    this_Pop_Data   = squeeze(this_Pop_Data)';
    this_FR_Diff    = this_Pop_Data(:,2:end) - this_Pop_Data(:, 1);
    this_Pop_Change = this_FR_Diff./this_Pop_Data(:, 1);
    this_Pop_Change = this_Pop_Change(sum(isinf(this_Pop_Change),2) < 1, :);
    h2 = errorbar([2*group_Trials : group_Trials : length(this_Pair)], nanmedian(this_Pop_Change), nanstd(this_Pop_Change)./sqrt(size(this_Pop_Change, 1)), 'or');
    h2.CapSize   = 0;
    h2.LineWidth = 2;
    h2.Color     = [1 0 0];
    aX          = gca;
    aX.TickDir  = 'out';
    aX.XTick    = 2*group_Trials : group_Trials : length(this_Pair);
    aX.XLim     = [150 450];
    aX.YLim     = YLIM;
    aX.Box      = 'off';
    
    if iPair == 1
        
        aX.XLabel.String = '# Trials';
        aX.YLabel.String = 'Index / Proportion of change';
        
    end
end

figure(5)
select_Param = 2;
for iPair = 1 : length(fr_Stat_Matrix)
    subplot(1,6, iPair)
    color_Ind    = 1;
    for iElectrode = 1 : length(select_Electrodes)
        
        xAXIS_Rand_Pos     = x_Scaling*randn;
        
        this_Cell_Data     = fr_Stat_Matrix{iPair}(select_Param, : , iElectrode);
        this_Cell_FR_Change= (this_Cell_Data(2:end) - this_Cell_Data(1))./this_Cell_Data(1);
        h1                 = plot([2*group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, this_Cell_FR_Change);
        h1.Color           = 0.7*[1 1 1];
        h1.LineWidth       = 0.5;
        hold on
        h2                 = scatter([2*group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, this_Cell_FR_Change, 20);
        h2.MarkerEdgeColor = 'none';
        h2.MarkerFaceColor = line_Color(color_Ind, :);
        h2.MarkerFaceAlpha = 0.5;
        
        %         color_Ind          = color_Ind + 1;
    end
    this_Pop_Data   = fr_Stat_Matrix{iPair}(select_Param, : ,:);
    this_Pop_Data   = squeeze(this_Pop_Data)';
    this_FR_Diff    = this_Pop_Data(:,2:end) - this_Pop_Data(:, 1);
    this_Pop_Change = this_FR_Diff./this_Pop_Data(:, 1);
    this_Pop_Change = this_Pop_Change(sum(isinf(this_Pop_Change),2) < 1, :);
    h2 = errorbar([2*group_Trials : group_Trials : length(this_Pair)], nanmedian(this_Pop_Change), nanstd(this_Pop_Change)./sqrt(size(this_Pop_Change, 1)), 'or');
    h2.CapSize   = 0;
    h2.LineWidth = 2;
    h2.Color     = [1 0 0];
    aX          = gca;
    aX.TickDir  = 'out';
    aX.XTick    = 2*group_Trials : group_Trials : length(this_Pair);
    aX.XLim     = [150 450];
    aX.YLim     = YLIM;
    aX.Box      = 'off';
    
    if iPair == 1
        
        aX.XLabel.String = '# Trials';
        aaX.YLabel.String = 'Index / Proportion of change';
        
    end
end

figure(6)
select_Param = 3;
for iPair = 1 : length(fr_Stat_Matrix)
    subplot(1,6, iPair)
    color_Ind    = 1;
    for iElectrode = 1 : length(select_Electrodes)
        
        xAXIS_Rand_Pos     = x_Scaling*randn;
        
        this_Cell_Data     = fr_Stat_Matrix{iPair}(select_Param, : , iElectrode);
        this_Cell_FR_Change= (this_Cell_Data(2:end) - this_Cell_Data(1))./this_Cell_Data(1);
        h1                 = plot([2*group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, this_Cell_FR_Change);
        h1.Color           = 0.7*[1 1 1];
        h1.LineWidth       = 0.5;
        hold on
        h2                 = scatter([2*group_Trials : group_Trials : length(this_Pair)]+xAXIS_Rand_Pos, this_Cell_FR_Change, 20);
        h2.MarkerEdgeColor = 'none';
        h2.MarkerFaceColor = line_Color(color_Ind, :);
        h2.MarkerFaceAlpha = 0.5;
        
        %         color_Ind          = color_Ind + 1;
    end
    this_Pop_Data   = fr_Stat_Matrix{iPair}(select_Param, : ,:);
    this_Pop_Data   = squeeze(this_Pop_Data)';
    this_FR_Diff    = this_Pop_Data(:,2:end) - this_Pop_Data(:, 1);
    this_Pop_Change = this_FR_Diff./this_Pop_Data(:, 1);
    this_Pop_Change = this_Pop_Change(sum(isinf(this_Pop_Change),2) < 1, :);
    h2 = errorbar([2*group_Trials : group_Trials : length(this_Pair)], nanmedian(this_Pop_Change), nanstd(this_Pop_Change)./sqrt(size(this_Pop_Change, 1)), 'or');
    h2.CapSize   = 0;
    h2.LineWidth = 2;
    h2.Color     = [1 0 0];
    aX          = gca;
    aX.TickDir  = 'out';
    aX.XTick    = 2*group_Trials : group_Trials : length(this_Pair);
    aX.XLim     = [150 450];
    aX.YLim     = YLIM;
    aX.Box      = 'off';
    
    if iPair == 1
        
        aX.XLabel.String = '# Trials';
        aX.YLabel.String = 'Index / Proportion of change';
        
    end
end

figure(7)
subInd = 1;
for iPair = 1 : length(stim_LeadStim)
    
    subplot(6,1,iPair)
    imshow([stim.allStimFile{stim_LeadStim(iPair)}{1}.im 255*ones(size(stim.allStimFile{stim_LeadStim(iPair)}{1}.im, 1), 20) stim.allStimFile{stim_TrailStim(iPair)}{1}.im], [])
    
end


