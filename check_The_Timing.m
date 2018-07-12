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
cStruct   = dat.Data.Comments;  % comments
comments1 = double([cStruct.TimeStamp])/tRes*1000;  % comment times (ms)

spikes    = double(dat.Data.Spikes.TimeStamp)/tRes*1000;  % spike times (ms)

% txt = reshape(NEV.Data.Comments.Text,[],92);
comment_txt        = reshape(cStruct.Comments,[],92);
[match, noMatch]   = regexp(cellstr(comment_txt(:, 1:22)),'\d','match','forceCellOutput');
trial_NumCellArray = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
% trial_NumArray     = cellfun(@str2num, trial_NumCellArray, 'UniformOutput', false);
% trial_NumArray     = trial_NumArray(~ismember(1:length(trial_NumArray), [3829]));
trial_NumArray     = cellfun(@str2num, trial_NumCellArray(~ismember(1:length(trial_NumCellArray), [2329  4119])), 'UniformOutput', false);
find(diff(cell2mat(trial_NumArray))'>1)+1

[match, noMatch]   = regexp(cellstr(comment_txt(:, 23:38)),'\d','match','forceCellOutput');
trial_LeadCellArray= cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_LeadCellArray = trial_LeadCellArray(~cellfun(@isempty,trial_LeadCellArray));
trial_LeadArray    = cellfun(@str2num, trial_LeadCellArray, 'UniformOutput', false);

[match, noMatch]   = regexp(cellstr(comment_txt(:, 45:52)),'\d','match','forceCellOutput');
trial_TrailCellArray= cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_TrailCellArray = trial_TrailCellArray(~cellfun(@isempty,trial_TrailCellArray));
trial_TrialArray    = cellfun(@str2num, trial_TrailCellArray, 'UniformOutput', false);

[match, noMatch]   = regexp(cellstr(comment_txt(:, 53:end)),'\d','match','forceCellOutput');
trial_SampleCellArray= cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_SampleCellArray = trial_SampleCellArray(~cellfun(@isempty,trial_SampleCellArray));
trial_SampleArray  = cellfun(@str2num, trial_SampleCellArray, 'UniformOutput', false);

comment_IDs        = [cell2mat(trial_LeadArray)'; cell2mat(trial_TrialArray)'; cell2mat(trial_SampleArray)']; 

if any(any(comment_IDs - stim.allStimTrain(:, :))) == true
    error('number of stim in stim file and comments doesnt match')
end
comments           = comments1(~ismember(1:length(comments1), [1 2329+1  4119+1]));

% Digital Timings
RawDIO        = dat.Data.SerialDigitalIO.UnparsedData;  % DIO tags
RawTimes      = double(dat.Data.SerialDigitalIO.TimeStamp)/tRes*1000;  % DIO digital time(ms)
DIO           = mod(RawDIO, 128);   % digital line without photodiode
stim_OnTime1  = RawTimes(DIO == 2);  % stim 1 onset
stim_OffTime1 = RawTimes(DIO == 3);  % stim 1 offset
stim_OffTime1 = stim_OffTime1([1:2097 2099:end]);

stim_OnTime2  = RawTimes(DIO == 4);  % stim 2 onset
stim_OffTime2 = RawTimes(DIO == 5);  % stim 1 offset

% Photodiode
PDTimes = double(dat.Data.Spikes.TimeStamp(dat.Data.Spikes.Electrode == 129))/tRes*1000;
PDTimes = PDTimes(PDTimes > comments(1));

% Channels information
electrodes        = unique(dat.Data.Spikes.Electrode);  % electrode numbers
% select_Electrodes = [1,2,3,4,5,6,7,8,9,10,11,12,14,19,21,22,26,27,29,31,32,37,40,41,42,44,46,50,51,52,53,54,55,56,57,58,62,63,65,66,67,73,75,76,81,83,84,85,86,87,94,95]; % 25
% select_Electrodes = [1:14 16:19 21 22 26 27 29 31 32 37 40:42 44 46 47 50:58 61:70 73 75 76 81 83:88 91 93:96 ]; % 26
select_Electrodes = [1:12 17 19 21 23 26 27 29 32 32 37 40 41 42 44 46 50 51:57 66 73 75 76 83 85 86 87]; % 28
%% extract stimulus information
if strcmpi(stim.textureType, 'texture')
    stim_LeadStim  = stim.TextFamilies(1:length(stim.TextFamilies)/2);  % leading stimulus names/indexes
    stim_TrailStim = stim.TextFamilies(1+(length(stim.TextFamilies)/2):end);  % trailing stimulus names/indexes
else
    stim_LeadStim  = 1:length(stim.oriList)/2;  % leading stimulus names/indexes
    stim_TrailStim = (1+(length(stim.oriList)/2)):length(stim.oriList);  % trailing stimulus names/indexes
end
stim_Train     = stim.allStimTrain;  % stimulus train. This should be a matrix of 3*n. 1st row: leading stim name/ind, 2nd trailing stim name/ind, last sample number

%% 
figure(1)
subplot(221)
hist(diff(comments),0:1500)
title('comments')

subplot(222)
hist(diff(stim_OnTime1),0:1500)
title('s1 onset')

subplot(223)
hist(diff(stim_OffTime2),0:1500)
title('s2 offset')

subplot(224)
hist(diff(stim_OffTime1),0:1500)
title('s1 offset')

