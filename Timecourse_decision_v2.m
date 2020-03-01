%I am trying to write a script that can generate the roc along the time.
%% Step1: I am trying to extract the id for Right/Left planning neuron
% clearvars -except Sum % load Sum2.mat (there is one neuron in Sum.mat which does not fire during delay, but captured)
% for i = 1:length(Sum)
%     Idx(i) = Sum(i).auROC.Rcorr_vsLcorr.stats.respFlag;
%     Idx_dir(i) = Sum(i).auROC.Rcorr_vsLcorr.p;
% end

function Timecourse_decision_v2(data,Idx,Idx_dir)
%         Input:
%               Sum: the data saved in Sum2.mat
%               Idx: the Index of selective neuron
%               Idx: the direction of the selectivity
Idx_rp = find(Idx==1 & Idx_dir>0); % id for right  neuron
Idx_lp = find(Idx==1 & Idx_dir<0); % id for left  neuron
%% Step2: Calculate the right direction preference across time ((auROC-0.5)*2)
type = {'Right','Left'};
clear PSTH_auROC_R PSTH_auROC_L

for i = 1:length(Idx_rp)
    PSTH_auROC_R(i,:) = psth_auROC_ke(data(Idx_rp(i)).(type{2}).Correct.scmatrix, data(Idx_rp(i)).(type{1}).Correct.scmatrix);
end
m_auROC= mean(PSTH_auROC_R,1);
sem_auROC = std(PSTH_auROC_R,1)/sqrt(size(PSTH_auROC_R,1));
time = data(1).Right.Correct.timepoint;
figure
h1= boundedline(time,m_auROC,sem_auROC,'r');
ylim([-0.4,0.4])
%% Step2: Calculate the Left direction preference across time ((auROC-0.5)*2)
type = {'Right','Left'};
for i = 1:length(Idx_lp)
    PSTH_auROC_L(i,:) = psth_auROC_ke(data(Idx_lp(i)).(type{2}).Correct.scmatrix, data(Idx_lp(i)).(type{1}).Correct.scmatrix);
end
m_auROC= mean((PSTH_auROC_L),1);
sem_auROC = std(PSTH_auROC_L,1)/sqrt(size(PSTH_auROC_L,1));
time = data(1).Right.Correct.timepoint;
hold on
h2 = boundedline(time,m_auROC,sem_auROC);
ylim([-0.3,0.3])
ylabel('Direction Planning Preference')
xlabel('Time (s)')
%% Step3:  generate the colormap
auROC_colormap_decision(time,PSTH_auROC_R, PSTH_auROC_L)