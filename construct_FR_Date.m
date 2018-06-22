function Param = construct_FR_Date(Param)


fprintf('\n Loading the data ...\n ');

load([Param.StimulusName])

if exist([Param.FileName(1:end-3) 'mat'], 'file')
    load([Param.FileName(1:end-3) 'mat'])
else
    openNEV([Param.FileName], 'read', 'nosave');
    NEV.Data.Spikes.Waveform = [];
    save([Param.FileName(1:end-3) 'mat'], 'NEV')
end

if exist('NEVdata', 'var')
    NEV   = NEVdata;
    clear   NEVdata;
end

%% Extract some event information and timing
dat = cbmex_Parse_data(NEV);
clear   NEV;
tRes      = dat.MetaTags.TimeRes;  % sampling resolution
cStruct   = dat.Data.Comments;  % comments
comments1 = double([cStruct.TimeStamp])/tRes*1000;  % comment times (ms)

spikes    = double(dat.Data.Spikes.TimeStamp)/tRes*1000;  % spike times (ms)

% txt = reshape(NEV.Data.Comments.Text,[],92);
% comment_txt        = reshape(cStruct.Comments,[],92);
% [match, noMatch]   = regexp(cellstr(comment_txt(:, 1:22)),'\d','match','forceCellOutput');
% trial_NumCellArray = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
% trial_NumArray     = cellfun(@str2num, trial_NumCellArray, 'UniformOutput', false);
% % find(diff(cell2mat(trial_NumArray))>1)+1;
% comments           = comments1(~[1 0 diff(cell2mat(trial_NumArray))'>1]);

comments      = comments1;
% Digital Timings
RawDIO        = dat.Data.SerialDigitalIO.UnparsedData;  % DIO tags
RawTimes      = double(dat.Data.SerialDigitalIO.TimeStamp)/tRes*1000;  % DIO digital time(ms)
DIO           = mod(RawDIO, 128);   % digital line without photodiode
stim_OnTime1  = RawTimes(DIO == 2);  % stim 1 onset
stim_OffTime1 = RawTimes(DIO == 3);  % stim 1 offset
stim_OnTime2  = RawTimes(DIO == 4);  % stim 2 onset
stim_OffTime2 = RawTimes(DIO == 5);  % stim 1 offset


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

sTrain = zeros(length(Param.SelectedChannels), ceil(max(spikes)));
for iElectrode = 1 : length(Param.SelectedChannels)
    sTrain(iElectrode, round(spikes(dat.Data.Spikes.Electrode == Param.SelectedChannels(iElectrode)))) = 1;
end

leadStimDuration  = 1000*stim.durationLeadStim;  % presentation time of leading stimulus
trailStimDuration = 1000*stim.durationTrailStim;  % presentation time of trailing stimulus
ISIDurartin       = 1000*stim.durationISI;  % ISI duration

% there can be 4 Alignments: start of 1st stim, end of 1st stim, start of
% ISI, comments
other_Alignments = [0 -leadStimDuration -(leadStimDuration + trailStimDuration + ISIDurartin)];

clc
for iElectrode = 1 : length(Param.SelectedChannels)
    fprintf(['Cell #:  ' num2str(Param.SelectedChannels(iElectrode))]);
    if Param.SDF  == true
        sdf = conv(ones(1, Param.SDF_binSize), sTrain(iElectrode,:))*(1/(Param.SDF_binSize/1000));
    else
        sdf = sTrain(iElectrode,:);
    end
    pair_Cnt = 1;
    for iTrailStim = 1 : length(stim_TrailStim)
        
        for iLeadStim = 1 : length(stim_LeadStim)
            
            if iTrailStim == iLeadStim
                this_Pair = find(stim_Train(1, :) == stim_LeadStim(iLeadStim) &  stim_Train(2, :) == stim_TrailStim(iTrailStim));
                
                switch Param.select_Alignments
                    case 1
                        % PSTH aligned to the start of first event
                        this_Stim_OnTime1  = round(stim_OnTime1(this_Pair));
                        this_Epochs = repmat(1:Param.PostStim, [length(this_Pair) 1]) + [this_Stim_OnTime1  + other_Alignments(1)-1]';
                        resps = sdf(this_Epochs);
                    case 2
                        % PSTH aligned to the end of first event
                        this_stim_OffTime1 = round(stim_OffTime1(this_Pair));
                        this_Epochs = repmat(1:Param.PostStim, [length(this_Pair) 1]) + [this_stim_OffTime1 + other_Alignments(2)-1]';
                        resps = sdf(this_Epochs);
                    case 3
                        % PSTH aligned to the start of ISI
                        this_stim_OffTime2 = round(stim_OffTime2(this_Pair));
                        this_Epochs = repmat(1:Param.PostStim, [length(this_Pair) 1]) + [this_stim_OffTime2 + other_Alignments(3)-1]';
                        resps = sdf(this_Epochs);
                    case 4
                        % PSTH aligned to the start of comments
                        this_comments      = round(comments(this_Pair));
                        this_Epochs = repmat(1:Param.PostStim, [length(this_Pair) 1]) + [this_comments      + other_Alignments(3)-1]';
                        resps = sdf(this_Epochs);
                end
                
                Param.all_Resp{pair_Cnt}(:,:, iElectrode) = resps;
                pair_Cnt = pair_Cnt +1 ;
            end
        end
        
        
    end
    
    fprintf(', ')
end

