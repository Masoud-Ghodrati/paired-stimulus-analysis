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
dat = cbmex_Parse_data(NEV);
clear   NEV;
tRes = dat.MetaTags.TimeRes;  % sampling resolution
cStruct = dat.Data.Comments;  % comments
comments1 = double([cStruct.TimeStamp])/tRes*1000;  % comment times (ms)

spikes = double(dat.Data.Spikes.TimeStamp)/tRes*1000;  % spike times (ms)

% txt = reshape(NEV.Data.Comments.Text,[],92);
comment_txt        = reshape(cStruct.Comments,[],92);
[match, noMatch]   = regexp(cellstr(comment_txt(:, 1:22)),'\d','match','forceCellOutput');
trial_NumCellArray = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_NumArray     = cellfun(@str2num, trial_NumCellArray, 'UniformOutput', false);
% find(diff(cell2mat(trial_NumArray))>1)+1;
comments           = comments1(~[1 0 diff(cell2mat(trial_NumArray))'>1]);


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
select_Electrodes = electrodes;  % which electrode(s) you want to analyze

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
close all
SDF_binSize       = 10;  % ms
time_Window       = 1000*[1 3000];  % time window ro plot spikes (sec)
line_Color        = [1 0 0; 0 1 0; 0 0 1; 0 0 0]; % line colors for tag
% Raster plot properties
YAXIS_LIM          = [0, 192];
RASTER_AXIS_INCRIM = 2;
RASTER_START_POINT = 50;
RASTER_LINE        = [YAXIS_LIM(1) + RASTER_START_POINT, YAXIS_LIM(1) + RASTER_START_POINT + RASTER_AXIS_INCRIM];
RASTER_COLOR       = 0.5*[1 1 1];
RASTER_LINE_WIDTH  = 0.1;
PDF_RESOLUTION = '-r300';
FILE_NAME = 'PSTH_AND_RASTER';
raster_Increment = 0;

figure
sdf = 0;
for iElectrode = 1 : length(select_Electrodes)
    
    this_Electrode_Spikes = spikes(dat.Data.Spikes.Electrode == select_Electrodes(iElectrode));
    this_Electrode_Spikes = this_Electrode_Spikes(this_Electrode_Spikes >= time_Window(1) &  this_Electrode_Spikes <= time_Window(2));
    this_Electrode_Spikes = this_Electrode_Spikes-time_Window(1);
    if sum(this_Electrode_Spikes <= 0.5) > 0
        this_Electrode_Spikes = this_Electrode_Spikes(this_Electrode_Spikes >= 0.5);
    end
    sTrain = zeros(1, diff(time_Window));
    sTrain(round(this_Electrode_Spikes)) = 1;
    sdf = sdf + conv(ones(1, SDF_binSize), sTrain)*(1/(SDF_binSize/1000));
    
    this_Raster = repmat(this_Electrode_Spikes, [2 1]);
%     line(this_Raster, RASTER_LINE + raster_Increment, 'Color', RASTER_COLOR, 'LineWidth', RASTER_LINE_WIDTH); hold on
    
    raster_Increment = raster_Increment + RASTER_AXIS_INCRIM;
    
end

plot(1:diff(time_Window), sdf(1:diff(time_Window))./length(select_Electrodes))

% stim 1 onset
this_Tag = stim_OnTime1(stim_OnTime1 >= time_Window(1) &  stim_OnTime1 <= time_Window(2));
this_Tag = this_Tag - time_Window(1);
this_Tag = repmat(this_Tag, [2 1]);
line(this_Tag', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(1, :), 'LineWidth', RASTER_LINE_WIDTH)

% stim 1 offset
this_Tag = stim_OffTime1(stim_OffTime1 >= time_Window(1) &  stim_OffTime1 <= time_Window(2));
this_Tag = repmat(this_Tag, [2 1]);
this_Tag = this_Tag - time_Window(1);
line(this_Tag', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(2, :), 'LineWidth', RASTER_LINE_WIDTH)

% stim 1 offset
this_Tag = stim_OffTime2(stim_OffTime2 >= time_Window(1) &  stim_OffTime2 <= time_Window(2));
this_Tag = repmat(this_Tag, [2 1]);
this_Tag = this_Tag - time_Window(1);
line(this_Tag', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(3, :), 'LineWidth', RASTER_LINE_WIDTH)

% comment
this_Tag = comments(comments >= time_Window(1) &  comments <= time_Window(2));
this_Tag = repmat(this_Tag, [2 1]);
this_Tag = this_Tag - time_Window(1);
line(this_Tag', [YAXIS_LIM(1) raster_Increment+RASTER_START_POINT], 'linestyle', ':', 'color', line_Color(4, :), 'LineWidth', RASTER_LINE_WIDTH);

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
