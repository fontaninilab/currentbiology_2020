%% load data
load('data_PSTH_taste.mat')
%% Let's construct the firing rate matrix for each neuron
f = fieldnames(data);
for i = 1:length(data)
    for j = 1:length(f)
        firing_rate.(f{j})(i,:) = data(i).(f{j}).FR_avg; 
    end    
end
%% PCA analysis
[coeff,score,latent,tsquared,explained,mu] = pca((firing_rate.unit_all(:,11:end))'); % taste sampling start from 11th bins
%% apply the PCA transformation to each case
S_traj = ((firing_rate.unit_S(:,11:end))'-mu)*coeff(:,1:3);   % taste sampling start from 11th bins
Q_traj = ((firing_rate.unit_Q(:,11:end))'-mu)*coeff(:,1:3);   % taste sampling start from 11th bins
M_traj = ((firing_rate.unit_M(:,11:end))'-mu)*coeff(:,1:3);   % taste sampling start from 11th bins
SO_traj = ((firing_rate.unit_SO(:,11:end))'-mu)*coeff(:,1:3); % taste sampling start from 11th bins

% Plot the population trajectory
figure;
h1 = plot3(S_traj(:,1),S_traj(:,2),S_traj(:,3),'-ro')
hold on
h2 = plot3(Q_traj(:,1),Q_traj(:,2),Q_traj(:,3),'-bo')
h3 = plot3(M_traj(:,1),M_traj(:,2),M_traj(:,3),'-mo')
h4 =plot3(SO_traj(:,1),SO_traj(:,2),SO_traj(:,3),'-ko')
box on
grid on
rotate3d on
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
legend([h1,h2,h3,h4],{'S','Q','M','SO'})

%% calculate eucleadian distance in PCA space
S_Qdist = sqrt(sum((S_traj - Q_traj).^2,2));    % Eucleadian distance between S and Q
M_SOdist = sqrt(sum((M_traj - SO_traj).^2,2));  % Eucleadian distance between M and SO


figure;
hold on
plot(data(1).unit_S.timepoint(11:end),S_Qdist,'-ro')
plot(data(1).unit_S.timepoint(11:end),M_SOdist,'-bo')
ylim([0,25])
ylabel('Eucleadian distance in PCA space')
xlabel('Time (s)')
legend({'S-Q distance','M-SO distance'})

%% Calculated the auROC normalized distance between tastatns
for i = 1:length(data)
    data(i).SM = abs(psth_auROC_ke(data(i).unit_S.scmatrix,data(i).unit_M.scmatrix));    % distance between S and M
    data(i).SQ = abs(psth_auROC_ke(data(i).unit_S.scmatrix,data(i).unit_Q.scmatrix));    % distance between S and Q
    data(i).MSO = abs(psth_auROC_ke(data(i).unit_M.scmatrix,data(i).unit_SO.scmatrix));  % distance between M and SO
    data(i).QSO = abs(psth_auROC_ke(data(i).unit_Q.scmatrix,data(i).unit_SO.scmatrix));  % distance between Q and SO
    fprintf('Finish processing neuron # %0.1f\n', i)
end

for i = 1:length(data)
    SQ(i,:) = data(i).SQ;
    MSO(i,:) = data(i).MSO;
    SM (i,:) = data(i).SM;
    QSO(i,:) = data(i).QSO;
end


similarAction = 1/2 * (SQ + MSO);            % distance between tastatns associated with the same actions
m_similarAction = mean(similarAction);
sem_similarAction = std(similarAction)/sqrt(size(similarAction,1));

similarQuality = 1/2 * (SM + QSO);           % distance between tastatns associated with the same quality
m_similarQuality = mean(similarQuality);
sem_similarQuality = std(similarQuality)/sqrt(size(similarQuality,1));

figure; 
h1 = boundedline(data(1).unit_all.timepoint, m_similarAction, sem_similarAction,'r')
hold on
h2 = boundedline(data(1).unit_all.timepoint, m_similarQuality, sem_similarQuality,'k')
xlim([0,2.5])
ylim([0.06,0.12])
xlabel('Time (s)')
ylabel('Normalized auROC')
legend([h1,h2],{'Taste pairs for tastants with same action','Taste pairs for stastants with same quality'})

