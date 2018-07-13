% SCRIPT
% LZ - July 18, 2014
clear all;
% close all
clc

% load the NEV file and do some pre-processing.
dPath = 'F:\CJ194\Data\';
dFile = 'CJ194_datafile013.nev';


load('F:\CJ194\Stimulus\CJ194_FullTex0005.mat')

if exist([dPath dFile(1:end-3) 'mat'], 'file')
    load([dPath dFile(1:end-3) 'mat'])
else
    openNEV([dPath dFile], 'read');
end

if exist('NEVdata', 'var')
    NEV = NEVdata; clear NEVdata;
end

dat = cbmex_Parse_data(NEV);
tRes = dat.MetaTags.TimeRes; % sampling resolution
cStruct = dat.Data.cbmexComment;

comments = double([cStruct.TimeStamp])/tRes*1000; %comment times (ms)
% comments = double([cStruct.TimeStamp])/tRes*1000; %comment times (ms)
spikes = double(dat.Data.Spikes.TimeStamp)/tRes*1000; %spike times (ms)

% Digital Timings
RawDIO = dat.Data.SerialDigitalIO.UnparsedData; %digital value
% RawTimes= double(dat.Data.SerialDigitalIO.TimeStamp);
RawTimes= double(dat.Data.SerialDigitalIO.TimeStamp)/tRes*1000; %digital time(ms)
DIO = mod(RawDIO, 128); %digital line without photodiode
T2 = RawTimes(DIO == 2); % stim onset
T1 = RawTimes(DIO == 1); % blank onset

% Photodiode
PDTimes = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode == 129))/tRes*1000;
PDTimes = PDTimes(PDTimes > comments(1));

%% ------- HERE IS WHERE YOU CREATE THE PARAMETERS YOU'RE ANALYSING BY ------- %%

p1List = stim.allTrials(3,:);
p2List = stim.allTrials(4,:);
scrmList = stim.allTrials(2,:);
scrmList(~sum([p1List;p2List])==0)= 0;
p1List2 = p1List(p1List>0);
p2List2 = p2List(p2List>0);
% set up parameters 1 and 2
p1 = unique(p1List2);
p2 = unique(p2List2);


%% ------- MAKE A SPIKE TRAIN ------- %%

sTrain = zeros(32, ceil(max(spikes)));
for ch = 1:96
    sTrain(ch, round(spikes(NEV.Data.Spikes.Electrode == ch))) = 1;
end
 
%% ------- PSTH FOR EVERY p1 and p2 PAIRING ------- %%
% p1 is different plots, p2 is different lines


figure('units','normalized','outerposition',[0 0 1 1]);
SDF_binSize = 10; % ms
winSize = 250; %ms (PSTH length)
for ch = 1 : 96,
    
    subplot(10,10,ch)
    sdf = conv(ones(1, SDF_binSize), sTrain(ch,:))*(1/(SDF_binSize/1000));
    
    resps1 = zeros(sum(p1List>0),winSize);
    resps2 = zeros(sum(p2List>0),winSize);
    resps3 = zeros(sum(scrmList>0),winSize);
    tList = round(T2); % for times, use Digital 1 (could also use comments)
    % t2 is the stim onset, t1 is blank onset
    % tList = comments(3:end);  % if using comments, don't use the first init
    % comments
    
    q1 =1;
    for t = 1:length(p1List)
        
        if p1List(t)>0
            
            resps1(q1, :) = sdf(tList(t):tList(t)+winSize-1);
            q1=q1+1;
        end
    end
    q1 =1;
    for t = 1:length(p2List)
        
        if p2List(t)>0
            
            resps2(q1, :) = sdf(tList(t):tList(t)+winSize-1);
            q1=q1+1;
        end
    end
    q1 =1;
    for t = 1:length(scrmList)
        
        if scrmList(t)>0
            
            resps3(q1, :) = sdf(tList(t):tList(t)+winSize-1);
            q1=q1+1;
        end
    end
    
    %
    % figure(1);
    % clf;
    % for con = 1:length(p1)
    %     subplot(length(p1),1, con);
    %     for d = 1:length(p2)
    %
    %         plot(resps{con, d}/nCond(con,d), 'Color', colMap(d,:), 'LineWidth', 2); hold on;
    %     end
    % end
    
    m=mean(resps1);
    s=std(resps1)/sqrt(length(resps1));
    ind1=1:length(m);
    ind2=ind1(end:-1:1);
    hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],[0.4 0.4 0.9]);
    alpha(0.8)
    set(h,'edgecolor','none')
    plot(mean(resps1))
    
    m=mean(resps2);
    s=std(resps2)/sqrt(length(resps2));
    ind1=1:length(m);
    ind2=ind1(end:-1:1);
    hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],[0.9 0.4 0.4]);
    alpha(0.8)
    set(h,'edgecolor','none')
    plot(mean(resps2),'r')
    
    m=mean(resps3);
    s=std(resps3)/sqrt(length(resps3));
    ind1=1:length(m);
    ind2=ind1(end:-1:1);
    hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],0.4*[1 1 1]);
    alpha(0.8)
    set(h,'edgecolor','none')
    plot(mean(resps3),'k')
    set(gcf,'color','w')
    xlim([0 winSize])
    %     xlabel('Time (ms)'); ylabel('Firing Rate sp/s');
    
end
