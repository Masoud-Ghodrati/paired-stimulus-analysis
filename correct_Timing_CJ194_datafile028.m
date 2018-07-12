function [stim_OnTime1, stim_OffTime1, stim_OffTime2, comments] = correct_Timing_CJ194_datafile028(stim_OnTime1, stim_OffTime1, stim_OffTime2, dat, stim)

tRes                  = dat.MetaTags.TimeRes;  % sampling resolution
cStruct               = dat.Data.Comments;  % comments
comments_Uncorrected  = double([cStruct.TimeStamp])/tRes*1000;  % comment times (ms)

comment_txt           = reshape(cStruct.Comments,[],92);
[match, ~]            = regexp(cellstr(comment_txt(:, 1:22)),'\d','match','forceCellOutput');
trial_NumCellArray    = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_NumArray        = cellfun(@str2num, trial_NumCellArray, 'UniformOutput', false);
find(diff(cell2mat(trial_NumArray))'>1)

[match, ~]            = regexp(cellstr(comment_txt(:, 23:38)),'\d','match','forceCellOutput');
trial_LeadCellArray   = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_LeadCellArray   = trial_LeadCellArray(~cellfun(@isempty,trial_LeadCellArray));
trial_LeadArray       = cellfun(@str2num, trial_LeadCellArray, 'UniformOutput', false);

[match, ~]            = regexp(cellstr(comment_txt(:, 45:52)),'\d','match','forceCellOutput');
trial_TrailCellArray  = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_TrailCellArray  = trial_TrailCellArray(~cellfun(@isempty,trial_TrailCellArray));
trial_TrialArray      = cellfun(@str2num, trial_TrailCellArray, 'UniformOutput', false);

[match, ~]            = regexp(cellstr(comment_txt(:, 53:end)),'\d','match','forceCellOutput');
trial_SampleCellArray = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_SampleCellArray = trial_SampleCellArray(~cellfun(@isempty,trial_SampleCellArray));
trial_SampleArray     = cellfun(@str2num, trial_SampleCellArray, 'UniformOutput', false);

comment_IDs           = [cell2mat(trial_LeadArray)'; cell2mat(trial_TrialArray)'; cell2mat(trial_SampleArray)']; 

if any(any(comment_IDs - stim.allStimTrain)) == true
    error('number of stim in stim file and comments doesnt match')
end
comments              = comments_Uncorrected([false ~ismember(1:length(trial_NumArray), [1396  4611])]);

% Digital Timings
RawDIO                = dat.Data.SerialDigitalIO.UnparsedData;  % DIO tags
RawTimes              = double(dat.Data.SerialDigitalIO.TimeStamp)/tRes*1000;  % DIO digital time(ms)
DIO                   = mod(RawDIO, 128);   % digital line without photodiode
stim_OnTime1          = RawTimes(DIO == 2);  % stim 1 onset

stim_OffTime1         = RawTimes(DIO == 3);  % stim 1 offset
stim_OffTime1         = stim_OffTime1([1:4609 4611:end]);

stim_OnTime2          = RawTimes(DIO == 4);  % stim 2 onset
stim_OffTime2         = RawTimes(DIO == 5);  % stim 1 offset
