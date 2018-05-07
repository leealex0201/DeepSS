function NEVAnalysis(NEV)

% electrode that I am interested in
ElectrodeNumb = 1;

%% Get waveforms and its indecies from this channel
[Wf_All,ChannelInd,NumbOfUnits,TimestampsAll] = GetAllWf(ElectrodeNumb,NEV);
ThisElectrodeInd = find(ChannelInd);
% ThisElectrodeInd is indecies of total event (all spiking timing) that
% falls in the "ElectrodeNumb" we assigned above.

% we will be only including those WFs that make up 90% of wf range
[WfAllUpBnd,WfAllLowBnd] = GetBnds(Wf_All,30,0.1);
TakeThese = GetTakeThese(Wf_All,WfAllUpBnd,WfAllLowBnd);
Wf_All = Wf_All(TakeThese,:);
TimestampsAll = TimestampsAll(TakeThese);
ThisElectrodeInd = ThisElectrodeInd(TakeThese);

%% Now start with subunits
TotalSpikingInd = [];
for i = 1:NumbOfUnits
    [ThisWf,SpikeInd] = GetWf(ElectrodeNumb,i,NEV);
    ThisTempElectrodInd = find(SpikeInd);
    
    % take only these
    TakeThese = GetTakeThese(ThisWf,WfAllUpBnd,WfAllLowBnd);
    
    ThisElectrodeSubUnitInd{i} = ThisTempElectrodInd(TakeThese);
    TotalSpikingInd = [TotalSpikingInd;ThisTempElectrodInd(TakeThese)];
end
TotalSpikingInd = sort(TotalSpikingInd);
IndexOfNoise = setdiff(ThisElectrodeInd,TotalSpikingInd);
% It's the same way. But now we found a spiking time information based on
% the electrode AND subunit number. Therefore, in sets perspective,
% ThisElectrodeSubUnitInd{1} + ThisElectrodeSubUnitInd{2} + IndexOfNoise =
% ThisElectrodeInd.

% spiking activity index in terms of "ThisElectrodeInd"
if NumbOfUnits > 0
    for i = 1:length(ThisElectrodeSubUnitInd)
        ThisSubUnit = ThisElectrodeSubUnitInd{i};
        [~,ThisSubInds,~] = intersect(ThisElectrodeInd,ThisSubUnit);
        SubInds{i} = ThisSubInds;
    end
end

% noise in terms of "ThisElectrodeInd"
[~,NoiseInds,~] = intersect(ThisElectrodeInd,IndexOfNoise);
% for i = 1:length(IndexOfNoise)
%     NoiseInds(i) = find(ThisElectrodeInd==IndexOfNoise(i));
% end

% PCA
[~,score] = pca(Wf_All);

FC = score(:,1); % First principal component
SC = score(:,2); % Second principal component

% peak-trough (max-min)
PT = GetPT(Wf_All);

%% Upto here, data is saved in ProgressSaved.mat
% Noise is thick black star
% First unit will be thick red
% Second unit will be thick blue
% load('ProgressSaved.mat')

NumbOfNoiseToPlot = 25; % 20 noise?
NumbOfUnitToPlot = 50; % 20 units

cols = hsv(NumbOfNoiseToPlot);

% noise first
SelectedNoiseInds = randperm(length(NoiseInds),NumbOfNoiseToPlot);

% % units next
% for i = 1:NumbOfUnits
%     ThisUnitInds{i} = randperm(length(SubInds{i}),NumbOfUnitToPlot);
% end

figure(1)
subplot(2,2,1)
ndhist(FC,SC,'bins',1);
hold on
% plot noise first
for i = 1:NumbOfNoiseToPlot
    subplot(2,2,1)
    plot(FC(NoiseInds(SelectedNoiseInds(i))),...
        SC(NoiseInds(SelectedNoiseInds(i))),'o','MarkerSize',10,...
        'MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:))
    subplot(2,2,2)
    plot(PT(NoiseInds(SelectedNoiseInds(i))),...
        TimestampsAll(NoiseInds(SelectedNoiseInds(i))),...
        'o','MarkerEdgeColor',cols(i,:),'MarkerFaceColor',cols(i,:),...
        'MarkerSize',10)
    hold on
    subplot(2,2,3:4)
    plot(linspace(1,16,48),Wf_All(NoiseInds(SelectedNoiseInds(i)),:),...
        'color',cols(i,:))
    hold on
end

% % plot units
% for i = 1:NumbOfUnitToPlot
%     subplot(2,2,1)
%     for ii = 1:NumbOfUnits
%         plot(FC(SubInds{ii}(ThisUnitInds{ii}(i))),...
%             SC(SubInds{ii}(ThisUnitInds{ii}(i))),'h','MarkerSize',7,...
%             'MarkerFaceColor',cols(ii,:),'MarkerEdgeColor',cols(ii,:))
%     end
%     subplot(2,2,2)
%     for ii = 1:NumbOfUnits
%         plot(PT(SubInds{ii}(ThisUnitInds{ii}(i))),...
%             TimestampsAll(SubInds{ii}(ThisUnitInds{ii}(i))),...
%             'o','MarkerEdgeColor',cols(ii,:),'MarkerFaceColor',cols(ii,:))
%         hold on
%     end
%     subplot(2,2,3:4)
%     for ii = 1:NumbOfUnits
%         plot(linspace(1,16,48),Wf_All(SubInds{ii}(ThisUnitInds{ii}(i)),:),...
%             'color',cols(ii,:))
%     end
% end
subplot(2,2,1)
set(gca,'XTickLabel',[],'YTickLabel',[])
subplot(2,2,2)
line([0 0],[0 max(TimestampsAll)],'color',[0 0 0],'LineWidth',4)
line([-100 max(PT)],[0 0],'color',[0 0 0],'LineWidth',4)
xlim([-100 max(PT)])
ylim([-(max(TimestampsAll)*.1) max(TimestampsAll)])
grid on
subplot(2,2,3:4)
xlim([1 16])
set(gca,'XTickLabel',[],'YTickLabel',[])
saveas(gcf,'ExampleImage','png')
end

function PT = GetPT(Wf)
for i = 1:size(Wf,1)
    ThisWf = Wf(i,:);
    PT(i) = max(ThisWf)-min(ThisWf);
end
end

function TakeThese = GetTakeThese(Wf,WfUpBnd,WfLowBnd)
TakeThese = [];
for i = 1:size(Wf,1)
    ThisWf = Wf(i,:);
    if max(ThisWf)<WfUpBnd && min(ThisWf)>WfLowBnd
        TakeThese = [TakeThese;i];
    end
end

end

function [UpBnd,LowBnd] = GetBnds(X,nbins,Thresh)

[HistVal,HistDom] = hist(X(:),nbins);
HistVal = HistVal./max(HistVal);
ThisInd = find(HistVal<Thresh);
for i = 1:length(ThisInd)-1
    if ThisInd(i+1)-ThisInd(i) ~= 1
        break
    end
end
LowBnd = HistDom(ThisInd(i));

for i = length(ThisInd):-1:2
    if ThisInd(i)-ThisInd(i-1) ~= 1
        break
    end
end
UpBnd = HistDom(ThisInd(i));

end

function [Wf,ThisInds,NumbOfUnits,Timestamps] = GetAllWf(ElectrodeNumb,NEV)
NEV.Data.Spikes.TimeStamp = NEV.Data.Spikes.TimeStamp';
NEV.Data.Spikes.Electrode = NEV.Data.Spikes.Electrode';
NEV.Data.Spikes.Unit = NEV.Data.Spikes.Unit';
NEV.Data.Spikes.Waveform = NEV.Data.Spikes.Waveform';

ThisInds = NEV.Data.Spikes.Electrode==ElectrodeNumb;
Wf = double(NEV.Data.Spikes.Waveform(ThisInds,:)); % waveforms

Timestamps = double(NEV.Data.Spikes.TimeStamp(ThisInds)); % timestamps 

are_real_units    = NEV.Data.Spikes.Unit > 0 & NEV.Data.Spikes.Unit < 255;
Electrodes       = unique([NEV.Data.Spikes.Electrode(are_real_units) ...
    NEV.Data.Spikes.Unit(are_real_units)], 'rows'); % electrodes
TempElectrodesArray = [];
for j = 1:size(Electrodes,1)
    if ismember(Electrodes(j,1),ElectrodeNumb)
        TempElectrodesArray = [TempElectrodesArray;Electrodes(j,:)];
    end
end
Electrodes = TempElectrodesArray;
NumbOfUnits = size(Electrodes,1);
end

function [Wf,spikes_of_this_unit] = GetWf(ElectrodeNumb,SubElectrodeNumb,NEV)
NEV.Data.Spikes.TimeStamp = NEV.Data.Spikes.TimeStamp';
NEV.Data.Spikes.Electrode = NEV.Data.Spikes.Electrode';
NEV.Data.Spikes.Unit = NEV.Data.Spikes.Unit';
NEV.Data.Spikes.Waveform = NEV.Data.Spikes.Waveform';
are_real_units    = NEV.Data.Spikes.Unit > 0 & NEV.Data.Spikes.Unit < 255;
Electrodes       = unique([NEV.Data.Spikes.Electrode(are_real_units) ...
    NEV.Data.Spikes.Unit(are_real_units)], 'rows'); % electrodes
TempElectrodesArray = [];
for j = 1:size(Electrodes,1)
    if ismember(Electrodes(j,1),ElectrodeNumb)
        TempElectrodesArray = [TempElectrodesArray;Electrodes(j,:)];
    end
end
Electrodes = TempElectrodesArray;
NEV.Electrodes = Electrodes; % insert electrodes in the structure
Eind = Electrodes(find(Electrodes(:,2)==SubElectrodeNumb),1);
Uind = Electrodes(find(Electrodes(:,2)==SubElectrodeNumb),2);
spikes_of_this_unit = NEV.Data.Spikes.Electrode == Eind & NEV.Data.Spikes.Unit == Uind;
Wf = double(NEV.Data.Spikes.Waveform(spikes_of_this_unit,:));

end