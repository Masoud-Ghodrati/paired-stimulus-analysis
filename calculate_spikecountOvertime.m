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
bad_Electrode = [13 15 16 17 18 20 23 24 25 28 30 33 34 35 36 38 39 43 45 47 48 49 59 60 61 64 68 69 70 71 72 74 77 78 79 80 82 88 89 90 91 92 93 96];
select_Electrodes = electrodes(~ismember(electrodes, bad_Electrode));  % which electrode(s) you want to analyze

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

group_Trials = 100;  % group every "group_Trails" trials to see the effect of learning
line_Color   = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
baseline_Window = round(1000*[stim.durationITI(1)-0.1 stim.durationITI(1)]);
FigureTab   = true;  % if ture it plots a figure with mutiple tabs, otheriwse, multiple figure

%% find the time window of analysis

if FigureTab
    %     figure('units','normalized','outerposition',[0 0 1 1]);
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
    mx = mean(all_Resp);
    [mymx1, mymxind1] = max(mx(1:100));
    [mymx2, mymxind2] = max(mx(150:end));
    plot(mymxind1, mymx1, 'ro')
    h = plot(mymxind2+150,     mymx2, 'ro');
    h.Parent.Box = 'off';
    h.Parent.TickDir = 'out';
    h.Parent.XLabel.String = 'Time (ms)';
    h.Parent.YLabel.String = 'Firing Rate (spk/s)';
    if FigureTab
        thistab.Title = ['Chn ' num2str(iElectrode)];
    else
        suptitle = ['Chn ' num2str(iElectrode)];
    end
    
    cell_Ind_Max(:, iElectrode) = [mymxind1; mymxind2+150];
end


%%

for iElectrode = 1 : length(select_Electrodes)
    
    sdf = sTrain(iElectrode,:);
    
    for iTrailStim = 1 : length(stim_TrailStim)
        
        for iLeadStim = 1 : length(stim_LeadStim)
            
            if iTrailStim == iLeadStim
                
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
                
                peak_1 = sum(resps(:, cell_Ind_Max(1, iElectrode)-30 : cell_Ind_Max(1, iElectrode)+30 ), 2);
                peak_2 = sum(resps(:, cell_Ind_Max(2, iElectrode)-30 : cell_Ind_Max(2, iElectrode)+30 ), 2);
                
                data_Resp{1, iTrailStim}(:,:, iElectrode) =  [peak_1'; peak_2'];
            end
            
        end
        
    end
    
end
figure
for iTrailStim = 1 : length(stim_TrailStim)
    subplot(1,6,iTrailStim)
    
    h = plot(smooth(mean(squeeze(data_Resp{1, iTrailStim}(1,:, :)),2), 1), 'r'); hold on
    h = plot(smooth(mean(squeeze(data_Resp{1, iTrailStim}(2,:, :)),2), 1), 'b');
    h.Parent.Box = 'off';
    h.Parent.TickDir = 'out';
    if iTrailStim == 1
        h.Parent.YLabel.String = 'spike count';
        h.Parent.XLabel.String = 'Trial';
    end
end
