function [stim_OnTime1, stim_OffTime1, stim_OffTime2, comments] = correct_Timing_CJ194_datafile025(stim_OnTime1, stim_OffTime1, stim_OffTime2, dat, stim)

tRes                  = dat.MetaTags.TimeRes;  % sampling resolution
cStruct               = dat.Data.Comments;  % comments
comments_Uncorrected  = double([cStruct.TimeStamp])/tRes*1000;  % comment times (ms)

comment_txt           = reshape(cStruct.Comments,[],92);
[match, ~]            = regexp(cellstr(comment_txt(:, 1:22)),'\d','match','forceCellOutput');
trial_NumCellArray    = cellfun(@cell2mat,match(2:end), 'UniformOutput', false);
trial_NumArray        = cellfun(@str2num, trial_NumCellArray, 'UniformOutput', false);
trial_NumArray        = trial_NumArray(~ismember(1:length(trial_NumArray), [832 1689 2459]));
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
comments              = comments_Uncorrected([false true diff(cell2mat(trial_NumArray))'==1]);
comments              = comments(~ismember(1:length(comments), [832, 833, 1689, 2460]));
% Digital Timings
RawDIO                = dat.Data.SerialDigitalIO.UnparsedData;  % DIO tags
RawTimes              = double(dat.Data.SerialDigitalIO.TimeStamp)/tRes*1000;  % DIO digital time(ms)
DIO                   = mod(RawDIO, 128);   % digital line without photodiode
stim_OnTime1          = RawTimes(DIO == 2);  % stim 1 onset
stim_OnTime1          = stim_OnTime1([1:2179 2181:end]);

stim_OffTime1         = RawTimes(DIO == 3);  % stim 1 offset

stim_OffTime2         = RawTimes(DIO == 5);  % stim 1 offset

temp_On               = stim_OnTime1([2431:2434]);
temp_Off1             = stim_OffTime1([2431:2434]);
temp_Off2             = stim_OffTime2([2431:2434]);

corrected_OnTime      = [temp_On(1)       temp_Off2(2)-200 temp_On(3)      temp_Off2(4)-200];
corrected_OffTime1    = [temp_Off1(1)     temp_Off2(2)-100 temp_On(3)+100  temp_Off1(4)    ];
corrected_OffTime2    = [temp_Off1(1)+100 temp_Off2(2)     temp_On(3)+200  temp_Off2(4)    ];

stim_OnTime1          = [stim_OnTime1(1:2430)  corrected_OnTime   stim_OnTime1(2435:end) ];
stim_OffTime1         = [stim_OffTime1(1:2430) corrected_OffTime1 stim_OffTime1(2435:end)];
stim_OffTime2         = [stim_OffTime2(1:2430) corrected_OffTime2 stim_OffTime2(2435:end)];

comments = [comments(1:831) stim_OffTime2(832:833) comments(833:2429) stim_OffTime2(2431:2435)  comments(2435:2455) stim_OffTime2(2457:2460) comments(2460:end) stim_OffTime2(3598:end)];
