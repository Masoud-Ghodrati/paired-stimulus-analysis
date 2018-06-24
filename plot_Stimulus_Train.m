clear
close all
clc

stimulus_Path = 'F:\CJ194\Stimulus\';
stimulus_FileName = 'Paired_Stimulus_File_CJ194_0001.mat';

load([stimulus_Path stimulus_FileName])

unique_LeadStim  = unique(stim.allStimTrain(1,:));
unique_TrailStim = unique(stim.allStimTrain(2,:));
num_Stim         = 500;
stim_Train_Ind   = 1;
stim_Counter     = 0;
for iLead = 1 : length(unique_LeadStim)
    for iTrail = 1 : length(unique_TrailStim)
        
        this_Pair = [];
        for iStim = 1 : num_Stim
            
            if stim.allStimTrain(1, iStim) == unique_LeadStim(iLead) && stim.allStimTrain(2, iStim) == unique_TrailStim(iTrail)
                this_Pair = [this_Pair iStim];
            end
            
        end
        stim_Counter = stim_Counter + length(this_Pair);
        
        if ~isempty(this_Pair)
            
            h                 = plot(this_Pair, stim_Train_Ind * ones(size(this_Pair)), 's'); hold on
            h.MarkerSize      = 5;
            h.Marker          = 's';
            h.MarkerEdgeColor = 'none';
            h.MarkerFaceColor = rand(1,3);
            
        end
        stim_Train_Ind = stim_Train_Ind + 1;
    end
end

aX         = gca;
aX.Box     = 'off';
aX.TickDir = 'out';


