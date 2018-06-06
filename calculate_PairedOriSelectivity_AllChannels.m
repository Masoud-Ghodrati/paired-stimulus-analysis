clear
close all
clc

% load the NEV file and do some pre-processing.
data_Path = 'F:\CJ194\Data\';
data_FileName = 'CJ194_datafile025.nev';

stimulus_Path = 'F:\CJ194\Stimulus\';
stimulus_FileName = 'Paired_Stimulus_File_CJ194_0001.mat';

load([stimulus_Path stimulus_FileName])

if exist([data_Path data_FileName(1:end-3) 'mat'], 'file')
    load([data_Path data_FileName(1:end-3) 'mat'])
else
    openNEV([data_Path data_FileName], 'read');
end

if exist('NEVdata', 'var')
    NEV   = NEVdata;
    clear   NEVdata;
end

%% Extract some event information and timing
dat = cbmex_Parse_data(NEV);
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
PDTimes = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode == 129))/tRes*1000;
PDTimes = PDTimes(PDTimes > comments(1));

% Channels information
electrodes = unique(dat.Data.Spikes.Electrode);  % electrode numbers
select_Electrodes = electrodes;  % which electrode(s) you want to analyze

%% extract stimulus information
if strcmpi(stim.textureType, 'texture')
    stim_LeadStim  = stim.TextFamilies(1:length(stim.TextFamilies)/2);  % leading stimulus names/indexes
    stim_TrailStim = stim.TextFamilies(1+(length(stim.TextFamilies)/2):end);  % trailing stimulus names/indexes
    stim_Family = stim.TextFamilies;
else
    stim_LeadStim  = 1:length(stim.oriList)/2;  % leading stimulus names/indexes
    stim_TrailStim = (1+(length(stim.oriList)/2)):length(stim.oriList);  % trailing stimulus names/indexes
    stim_Family    = stim.oriList;
end
stim_Train     = stim.allStimTrain;  % stimulus train. This should be a matrix of 3*n. 1st row: leading stim name/ind, 2nd trailing stim name/ind, last sample number
stim_Images    = stim.allStimFile;  % presented image file

%% make a spike train for each selected channel

sTrain = zeros(length(select_Electrodes), ceil(max(spikes)));
for iElectrode = 1 : length(select_Electrodes)
    sTrain(iElectrode, round(spikes(NEV.Data.Spikes.Electrode == select_Electrodes(iElectrode)))) = 1;
end

%% PSHT for a n*n pairing matrix for each selected channel individually
close all

SDF_binSize       = 20;  % ms
leadStimDuration  = 1000*stim.durationLeadStim;  % presentation time of leading stimulus
trailStimDuration = 1000*stim.durationTrailStim;  % presentation time of trailing stimulus
ISIDurartin       = 1000*stim.durationISI;  % ISI duration

winSize           = leadStimDuration + trailStimDuration + ISIDurartin +400;  % ms (PSTH length)

% there can be 4 Alignments: start of 1st stim, end of 1st stim, start of
% ISI, comments
other_Alignments = [0 winSize;
    -leadStimDuration winSize-leadStimDuration;
    -(leadStimDuration + trailStimDuration + ISIDurartin) winSize-(leadStimDuration + trailStimDuration + ISIDurartin)];

group_Trials = 100;  % group every "group_Trails" trials to see the effect of learning
line_Color   = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
FigureTab   = false;  % if ture it plots a figure with mutiple tabs, otheriwse, multiple figure
temp_Window = [120 250];

%%
close all
if FigureTab
    figure(1);
    tab_group = uitabgroup; % tabgroup
end
select_Alignments = 2;
line_Color2 = colormap('parula');
line_width  = 1;
for iElectrode = 1 : length(select_Electrodes)
    
    if FigureTab
        thistab = uitab(tab_group);  % build a tab
        axes('Parent', thistab); % somewhere to plot
    else
        figure(iElectrode);
    end
    sdf = conv(ones(1, SDF_binSize), sTrain(select_Electrodes(iElectrode),:))*(1/(SDF_binSize/1000));
    
    iPanel = 1;
    
    for iTrailStim = 1 : length(stim_TrailStim)
        subplot(6, 1, iPanel)
        tuning_Resp = [];
        for iLeadStim = 1 : length(stim_LeadStim)
            
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
                    resps = sdf(this_Epochs);
                case 2
                    
                    % PSTH aligned to the end of first event
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime1 + other_Alignments(2, 1)-1]';
                    resps = sdf(this_Epochs);
                case 3
                    
                    % PSTH aligned to the start of ISI
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_stim_OffTime2 + other_Alignments(3, 1)-1]';
                    resps = sdf(this_Epochs);
                case 4
                    % PSTH aligned to the start of comments
                    this_Epochs = repmat(1:winSize, [length(this_Pair) 1]) + [this_comments       + other_Alignments(3, 1)-1]';
                    resps = sdf(this_Epochs);
            end
            %             mean_resp = mean(resps(:, temp_Window(1): temp_Window(2)), 2);
            %             tuning_Resp(:, iTrailStim) = [mean(mean_resp) std(mean_resp)./sqrt(length(mean_resp)) stim_Family(stim_TrailStim(iTrailStim))];
            tuning_Resp = [tuning_Resp; mean(resps(:, temp_Window(1): temp_Window(2)))];
        end
        %         h = errorbar(tuning_Resp(3, :), tuning_Resp(1, :), tuning_Resp(2, :));
        %         h.Parent.Box = 'off';
        %         h.Parent.TickDir = 'out';
        imagesc(tuning_Resp)
        iPanel = iPanel + 1;
    end
    
    
    if FigureTab
        thistab.Title = ['Chn ' num2str(iElectrode)];
    else
        suptitle = ['Chn ' num2str(iElectrode)];
    end
    
end