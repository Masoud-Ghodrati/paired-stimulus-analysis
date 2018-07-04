clear
close all
clc

% load the NEV file and do some pre-processing.
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

%% psth setting
leadStimDuration  = 1000*stim.durationLeadStim;  % presentation time of leading stimulus
trailStimDuration = 1000*stim.durationTrailStim;  % presentation time of trailing stimulus
ISIDurartin       = 1000*stim.durationISI;  % ISI duration
pre_Stim          = 1000;  % time before stimulus onset
post_Stim         = 1000;  % time after stimulus onset
winSize           = leadStimDuration + trailStimDuration + ISIDurartin + pre_Stim + post_Stim;  % ms (PSTH length)

% there can be 4 Alignments: start of 1st stim, end of 1st stim, start of
% ISI, comments
other_Alignments  = [-pre_Stim; -(leadStimDuration + pre_Stim); -(leadStimDuration + trailStimDuration + ISIDurartin + pre_Stim)];
group_Trials      = 200;  % group every "group_Trails" trials to see the effect of learning
select_Alignments = 1;

%% make a spike train for each selected channel
close all
SDF_binSize        = 50;  % ms
time_Window        = 1000*[0 2522];  % time window ro plot spikes (sec)
line_Color         = [1 0 0; 0 1 0; 0 0 1; 0 0 0]; % line colors for tag
% Raster plot properties
YAXIS_LIM          = [0, 60];
RASTER_AXIS_INCRIM = 0.1;
RASTER_START_POINT = 40;
RASTER_LINE        = [YAXIS_LIM(1) + RASTER_START_POINT, YAXIS_LIM(1) + RASTER_START_POINT + RASTER_AXIS_INCRIM];
RASTER_COLOR       = 0.0*[1 1 1];
RASTER_LINE_WIDTH  = 0.1;
PDF_RESOLUTION     = '-r300';
FILE_NAME          = 'PSTH_AND_RASTER';
raster_Increment   = 0;

figure
sdf = 0;
for iElectrode = 1 : length(electrodes)
    
    this_Electrode_Spikes = spikes(dat.Data.Spikes.Electrode == electrodes(iElectrode));
    this_Electrode_Spikes = this_Electrode_Spikes(this_Electrode_Spikes >= time_Window(1) &  this_Electrode_Spikes <= time_Window(2));
    this_Electrode_Spikes = this_Electrode_Spikes - time_Window(1);
    
    if sum(this_Electrode_Spikes <= 0.5) > 0
        this_Electrode_Spikes = this_Electrode_Spikes(this_Electrode_Spikes >= 0.5);
    end
    sTrain = zeros(1, diff(time_Window));
    sTrain(round(this_Electrode_Spikes)) = 1;
    sdf = sdf + conv(sTrain, ones(1, SDF_binSize), 'same')*(1/(SDF_binSize/1000));
    
    this_Raster = repmat(this_Electrode_Spikes, [2 1]);
    %     line(this_Raster, RASTER_LINE + raster_Increment, 'Color', RASTER_COLOR, 'LineWidth', RASTER_LINE_WIDTH); hold on
    plot(this_Raster(1,:), RASTER_LINE(1) + raster_Increment + ones(size(this_Raster(1,:))), '.',...
        'markersize',1,'Color', RASTER_COLOR, 'LineWidth', RASTER_LINE_WIDTH); hold on
    
    raster_Increment = raster_Increment + RASTER_AXIS_INCRIM;
    
end

% plot(1:diff(time_Window), sdf(1:diff(time_Window))./length(electrodes))
plot(1:diff(time_Window), sdf./length(electrodes))

% stim 1 onset
this_Tag_On1 = stim_OnTime1(stim_OnTime1 >= time_Window(1) &  stim_OnTime1 <= time_Window(2));
this_Tag_On1 = this_Tag_On1 - time_Window(1);
this_Tag_On1 = repmat(this_Tag_On1, [2 1]);
line(this_Tag_On1', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(1, :), 'LineWidth', RASTER_LINE_WIDTH)

% stim 1 offset
this_Tag_Off1 = stim_OffTime1(stim_OffTime1 >= time_Window(1) &  stim_OffTime1 <= time_Window(2));
this_Tag_Off1 = repmat(this_Tag_Off1, [2 1]);
this_Tag_Off1 = this_Tag_Off1 - time_Window(1);
line(this_Tag_Off1', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(2, :), 'LineWidth', RASTER_LINE_WIDTH)

% stim 1 offset
this_Tag_Off2 = stim_OffTime2(stim_OffTime2 >= time_Window(1) &  stim_OffTime2 <= time_Window(2));
this_Tag_Off2 = repmat(this_Tag_Off2, [2 1]);
this_Tag_Off2 = this_Tag_Off2 - time_Window(1);
line(this_Tag_Off2', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(3, :), 'LineWidth', RASTER_LINE_WIDTH)

% comment
this_Tag_Cmt = comments(comments >= time_Window(1) &  comments <= time_Window(2));
this_Tag_Cmt = repmat(this_Tag_Cmt, [2 1]);
this_Tag_Cmt = this_Tag_Cmt - time_Window(1);
line(this_Tag_Cmt', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(4, :), 'LineWidth', RASTER_LINE_WIDTH);

aX = gca;
aX.Box = 'off';
aX.TickDir = 'out';
aX.XTick = 0:1000: diff(time_Window);
aX.XTickLabel = [time_Window(1):1000: time_Window(2)]/1000;
aX.XLabel.String = 'Time (s)';
% aX.YTick = RASTER_AXIS_INCRIM*select_Electrodes;
aX.YTickLabel = [];
aX.YLabel.String = 'Channels';
aX.YLim = [YAXIS_LIM(1) RASTER_START_POINT + raster_Increment];

%% plot population psth

figure
line_Color = colormap('parula');
color_Ind = 1;
sTrain = zeros(1, ceil(max(spikes)));
sdf    = 0;
for iElectrode = 1 : length(select_Electrodes)
    
    sTrain(round(spikes(dat.Data.Spikes.Electrode == select_Electrodes(iElectrode)))) = 1;
    sdf = sdf + conv(ones(1, SDF_binSize), sTrain)*(1/(SDF_binSize/1000));
    
end

sdf = sdf./length(select_Electrodes);
inc_Ind = floor(length(line_Color)/(length(stim.allStimTrain)/group_Trials));
legend_Labels = [];
for iTrial = 1 : group_Trials: length(stim.allStimTrain)
    
    this_Trials = iTrial : iTrial + group_Trials-1;
    legend_Labels = [ legend_Labels {['Tr: ' num2str(this_Trials(1)) ' - ' num2str(this_Trials(end))]}];
    
    resps = [];
    switch select_Alignments
        case 1
            % PSTH aligned to the start of first event
            this_Stim_OnTime1  = round(stim_OnTime1(this_Trials));
            this_Epochs = repmat(1:winSize, [length(this_Trials) 1]) + [this_Stim_OnTime1  + other_Alignments(1)-1]';
            resps = sdf(this_Epochs);
        case 2
            % PSTH aligned to the end of first event
            this_stim_OffTime1 = round(stim_OffTime1(this_Trials));
            this_Epochs = repmat(1:winSize, [length(this_Trials) 1]) + [this_stim_OffTime1 + other_Alignments(2)-1]';
            resps = sdf(this_Epochs);
        case 3
            % PSTH aligned to the start of ISI
            this_stim_OffTime2 = round(stim_OffTime2(this_Trials));
            this_Epochs = repmat(1:winSize, [length(this_Trials) 1]) + [this_stim_OffTime2 + other_Alignments(3)-1]';
            resps = sdf(this_Epochs);
        case 4
            % PSTH aligned to the start of comments
            this_comments      = round(comments(this_Trials));
            this_Epochs = repmat(1:winSize, [length(this_Trials) 1]) + [this_comments      + other_Alignments(3)-1]';
            resps = sdf(this_Epochs);
    end
    
    
    h = plot(1:winSize, mean(resps), 'color', line_Color(color_Ind, :)); hold on
    h.LineWidth = 1;
    color_Ind = color_Ind  + inc_Ind;
    
end
legend(legend_Labels)
legend boxoff
h.Parent.Box = 'off';
h.Parent.TickDir = 'out';
h.Parent.XLabel.String = 'Time (ms)';
h.Parent.YLabel.String = 'Firing Rate (spk/s)';

