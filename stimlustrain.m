clear
close all
clc
load(['Z:\Shared\Marmoset data\CJ194 Stimulus Files\Paired_Stimulus_File_CJ194_0001']);

Lead  = unique(stim.allStimTrain(1,:));
Trail = unique(stim.allStimTrain(2,:));
numofStim = 500;
t = 1;
c = colormap;
u=0;
for iLead = 1 : length(Lead)
    for iTrail = 1 : length(Trail)
        a = [];
        
        for iStim = 1 : numofStim
            
            if stim.allStimTrain(1,iStim)==Lead(iLead) && stim.allStimTrain(2,iStim)==Trail(iTrail)
                a = [a iStim];
            end
            
        end
        u = u+length(a);
        if ~isempty(a)
            h = plot(a, t*ones(size(a)), 's'); hold on
            h.MarkerSize = 5;
            h.Marker = 's';
            h.MarkerEdgeColor = 'none';
            h.MarkerFaceColor = rand(1,3);
            
        end
        t = t + 1;
    end
end

aX = gca;
aX.Box = 'off';
aX.TickDir = 'out';


