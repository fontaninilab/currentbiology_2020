%% load
load('data_PSTH_direction.mat')
%% Calculate direction index and infer statistical significance.
for j = 1:length(data)
    fprintf('Process Neuron # %d\n', j)
    indx = find(data(1).Right.Correct.timepoint>-1 & data(1).Right.Correct.timepoint<0); % analyze window [-1,0]
    data(j).auROC.PlanRcorr_vsLcorr.Right.FRmatrix = mean(data(j).Right.Correct.FRmatrix(:,indx),2);
    data(j).auROC.PlanRcorr_vsLcorr.Left.FRmatrix = mean(data(j).Left.Correct.FRmatrix(:,indx),2);
    [data(j).auROC.PlanRcorr_vsLcorr.p]     = Test_auROC_dISCRIMINATION  (data(j).auROC.PlanRcorr_vsLcorr.Right.FRmatrix,...
        data(j).auROC.PlanRcorr_vsLcorr.Left.FRmatrix); % calculate the direction index
    
    %start permutation test
    [data(j).auROC.PlanRcorr_vsLcorr.shuffl.allp] = [data(j).auROC.PlanRcorr_vsLcorr.Right.FRmatrix;data(j).auROC.PlanRcorr_vsLcorr.Left.FRmatrix];
    
    % prepare for shuffle to compute statistic
    for jj =1:1000
        [data(j).auROC.PlanRcorr_vsLcorr.shuffl.Right(:,jj), idx]= datasample(data(j).auROC.PlanRcorr_vsLcorr.shuffl.allp,...
            size(data(j).auROC.PlanRcorr_vsLcorr.Right.FRmatrix,1),'Replace',false);
        idx1 = 1:size(data(j).auROC.PlanRcorr_vsLcorr.shuffl.allp,1); % shuffling
        data(j).auROC.PlanRcorr_vsLcorr.shuffl.Left(:,jj) = data(j).auROC.PlanRcorr_vsLcorr.shuffl.allp(setdiff(idx1,idx));
        data(j).auROC.PlanRcorr_vsLcorr.shuffl.p (:,jj)   = Test_auROC_dISCRIMINATION(data(j).auROC.PlanRcorr_vsLcorr.shuffl.Right(:,jj),...
            data(j).auROC.PlanRcorr_vsLcorr.shuffl.Left(:,jj)); % calculate the index for shuffled data
    end
    if data(j).auROC.PlanRcorr_vsLcorr.p<0 % calculate p value for stats.
        [data(j).auROC.PlanRcorr_vsLcorr.stats.pvalue] = length(find(data(j).auROC.PlanRcorr_vsLcorr.shuffl.p<data(j).auROC.PlanRcorr_vsLcorr.p))/1000;
    else
        [data(j).auROC.PlanRcorr_vsLcorr.stats.pvalue] = length(find(data(j).auROC.PlanRcorr_vsLcorr.shuffl.p>data(j).auROC.PlanRcorr_vsLcorr.p))/1000;
    end
    
    if data(j).auROC.PlanRcorr_vsLcorr.stats.pvalue<0.01 % criteria for significant index
        data(j).auROC.PlanRcorr_vsLcorr.stats.respFlag=1;
    else
        data(j).auROC.PlanRcorr_vsLcorr.stats.respFlag=0;
    end
    
    if data(j).auROC.PlanRcorr_vsLcorr.p==0
        data(j).auROC.PlanRcorr_vsLcorr.stats.respFlag=0;
    end
    
end
%% Plot the histogram of the direction index
a=1;
b=1;
c=1;
d=1;
e=1;
figure(400)
colors=distinguishable_colors(5);
for kk = 1:size(data,2)
    if data(kk).auROC.PlanRcorr_vsLcorr.stats.respFlag==0
        her2=plot(data(kk).auROC.PlanRcorr_vsLcorr.p,kk); hold on;
        set(her2,'Color',colors(3,:),'Marker','o','MarkerSize',6,...
            'MarkerEdgeColor',[.7 .7 .7],'MarkerFaceColor' ,colors(1,:) );
        %         pp(b)=data(kk).auROC.Rcorr_vsLcorr.p; % this is a error
        p(a)=data(kk).auROC.PlanRcorr_vsLcorr.p;
        a=a+1;
    else
        her3=plot(data(kk).auROC.PlanRcorr_vsLcorr.p,kk);hold on;
        set(her3,'Color',colors(3,:),'Marker','o','MarkerSize',6,...
            'MarkerEdgeColor',[.7 .7 .7],'MarkerFaceColor' ,colors(3,:) );
        %         pp(b)=data(kk).auROC.Rcorr_vsLcorr.p; % this is a error
        % Ke modified the followin part on 2/13/2019
        if data(kk).auROC.PlanRcorr_vsLcorr.stats.pvalue ==0 % pay attention to the case that stats p ==0
            shu = data(kk).auROC.PlanRcorr_vsLcorr.shuffl.p;
            p_ok = length(find(data(kk).auROC.PlanRcorr_vsLcorr.p>=shu))/1000; % recalculate the stats p
            p_oko = length(find(data(kk).auROC.PlanRcorr_vsLcorr.p<=shu))/1000; % recalculate the stats p
            if p_ok>=0.01 & p_oko >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
                fprintf('Neuron %4.2f is not good\n',kk)
                data(kk).auROC.PlanRcorr_vsLcorr.stats.respFlag =0;
                p(a)=data(kk).auROC.PlanRcorr_vsLcorr.p;
                a = a+1;
            else
                pp(b) =data(kk).auROC.PlanRcorr_vsLcorr.p;
                idPlan(b)=kk;
                if pp(b)<0
                    idPlanL(d)=kk;
                    d=d+1;
                elseif pp(b)>0
                    idPlanR(e)=kk;
                    e=e+1;
                end
                b=b+1;
            end
        else
            pp(b) =data(kk).auROC.PlanRcorr_vsLcorr.p;
            idPlan(b)=kk;
            if pp(b)<0
                idPlanL(d)=kk;
                d=d+1;
            elseif pp(b)>0
                idPlanR(e)=kk;
                e=e+1;
            end
            b=b+1;
        end
    end
    P(kk)=data(kk).auROC.PlanRcorr_vsLcorr.p;
end
hold on;
%figure settings
plot([0 0],[0 kk],'--k');
xlim([-1 1]);
ylim([0 kk]);
yticks([0 kk]);
yticklabels({'0',num2str(kk)});
ylabel('Neuron ID');
xlabel('Direction Index');
xticks([-1 0 1]);
xticklabels({'-1','0','1'});
box('off');
title('Direction-selective neurons')
set(gca,'fontsize',8);
%
figure(401);
h = histogram(p,30);hold on;h1=histogram(pp,30);
bin      = 0.06;
%histogram color settings
set(h,'FaceColor',colors(1,:),'EdgeColor','w','BinWidth',bin);
set(h1,'FaceColor',colors(3,:),'EdgeColor','w','BinWidth',bin);
% x axes settings
xlabel('Planning Preference');
xlim([-1 1]);
xticks([-1 -0.5 0 0.5 1]);
% y axes settings
yticks([0 max(h.Values)]);
ylabel('Neurons');
yticklabels({'0',num2str(max(h.Values))});
ylabel('Neuron ID');

title ('Population histogram of planning preference');
text(-0.95,15,{'Planning preference is calculated';...
    'using an ROC analysis comparing';...
    'right and left correct trials.'},'FontSize',5);
legend('NS','p<0.01');
box('off');
set(gca,'fontsize',8);
%% Time Course of direction index
for i = 1:length(data)
    Idx(i) = data(i).auROC.PlanRcorr_vsLcorr.stats.respFlag;
    Idx_dir(i) = data(i).auROC.PlanRcorr_vsLcorr.p;
end
%Plot the population auROC for planning
Timecourse_decision_v2(data,Idx,Idx_dir)
title('Planning')

%% try to calculate taste selectivity between M and SO for neurons showing direction selectivity
clear resp
for i = 1:length(data)
    resp(i) = data(i).auROC.PlanRcorr_vsLcorr.stats.respFlag; % get index of neurons showing direction selectivity
end
idx_p = find(resp ==1);
taste =[];
for j = 1:length(idx_p)
    fprintf('Process neuron %4.2f \n',j)
    indx = find(data(1).Right.Mplanning.timepoint>-1 & data(1).Right.Mplanning.timepoint<0); % analyze window [-1,0]
    taste(j).Mplanning.FRmatrix = mean(data(idx_p(j)).Right.Mplanning.FRmatrix(:,indx),2);
    taste(j).Oplanning.FRmatrix =mean(data(idx_p(j)).Right.Oplanning.FRmatrix(:,indx),2);
    [taste(j).MO.p]    = Test_auROC_dISCRIMINATION(taste(j).Mplanning.FRmatrix, taste(j).Oplanning.FRmatrix); % calculate taste selectivity
    
    % prepare for shuffle to compute statistic
    [taste(j).MO.allp] = [taste(j).Mplanning.FRmatrix;taste(j).Oplanning.FRmatrix];
    for jj =1:1000
        [taste(j).MO.shuffl.M(:,jj), idx]= datasample(taste(j).MO.allp,size(taste(j).Mplanning.FRmatrix,1),'Replace',false);
        idx1 = 1:size(taste(j).MO.allp,1);
        taste(j).MO.shuffl.O(:,jj)       =  taste(j).MO.allp(setdiff(idx1,idx));
        taste(j).MO.shuffl.p (:,jj)       = Test_auROC_dISCRIMINATION(taste(j).MO.shuffl.M(:,jj),taste(j).MO.shuffl.O(:,jj));
    end
    if taste(j).MO.p<0
        [taste(j).MO.stats.pvalue] = length(find(taste(j).MO.shuffl.p<taste(j).MO.p))/1000;
    else
        [taste(j).MO.stats.pvalue] = length(find(taste(j).MO.shuffl.p>taste(j).MO.p))/1000;
    end
    
    if taste(j).MO.stats.pvalue<0.01
        taste(j).MO.stats.respFlag=1;
    else
        taste(j).MO.stats.respFlag=0;
    end
    if taste(j).MO.p==0
        taste(j).MO.stats.respFlag=0;
    end
end
%%  try to calculate taste selectivity between S and Q
for j = 1:length(idx_p)
    fprintf('Process neuron %4.2f \n',j)
    indx = find(data(1).Left.Splanning.timepoint>-1 & data(1).Left.Splanning.timepoint<0); % analyze window [-1,0]
    taste(j).Splanning.FRmatrix = mean(data(idx_p(j)).Left.Splanning.FRmatrix(:,indx),2);
    taste(j).Qplanning.FRmatrix =mean(data(idx_p(j)).Left.Qplanning.FRmatrix(:,indx),2);
    
    [taste(j).SQ.p]    = Test_auROC_dISCRIMINATION  (taste(j).Splanning.FRmatrix, taste(j).Qplanning.FRmatrix);
    
    % prepare for shuffle to compute statistic
    % put in a single vector the p values from Right and Left
    [taste(j).SQ.allp] = [taste(j).Splanning.FRmatrix;taste(j).Qplanning.FRmatrix];
    for jj =1:1000
        [taste(j).SQ.shuffl.S(:,jj), idx]= datasample(taste(j).SQ.allp,size(taste(j).Splanning.FRmatrix,1),'Replace',false);
        idx1 = 1:size(taste(j).SQ.allp,1);
        taste(j).SQ.shuffl.Q(:,jj)       =  taste(j).SQ.allp(setdiff(idx1,idx));
        taste(j).SQ.shuffl.p (:,jj)       = Test_auROC_dISCRIMINATION(taste(j).SQ.shuffl.S(:,jj),taste(j).SQ.shuffl.Q(:,jj));
    end
    if taste(j).SQ.p<0
        [taste(j).SQ.stats.pvalue] = length(find(taste(j).SQ.shuffl.p<taste(j).SQ.p))/1000;
    else
        [taste(j).SQ.stats.pvalue] = length(find(taste(j).SQ.shuffl.p>taste(j).SQ.p))/1000;
    end
    
    if taste(j).SQ.stats.pvalue<0.01
        taste(j).SQ.stats.respFlag=1;
    else
        taste(j).SQ.stats.respFlag=0;
    end
    if taste(j).SQ.p==0
        taste(j).SQ.stats.respFlag=0;
    end
end
%% try to summarize to see which neurons show taste selectivity
for i = 1:length(taste)
    if taste(i).MO.stats.pvalue ==0
        shu = taste(i).MO.shuffl.p;
        ok_p = length(find(taste(i).MO.p >= shu))/1000;
        oko_p= length(find(taste(i).MO.p <= shu))/1000;
        if ok_p >=0.01 & oko_p >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
            fprintf('Neuron %4.2f is not good Left\n',i)
            taste(i).MO.stats.respFlag = 0;
        end
        
    elseif taste(i).SQ.stats.pvalue ==0
        shu2 = taste(i).SQ.shuffl.p;
        ok_p = length(find(taste(i).SQ.p >= shu2))/1000;
        oko_p= length(find(taste(i).SQ.p <= shu2))/1000;
        if ok_p >=0.01 & oko_p >=0.01 % if this exist, it means shuffling preference and real preference has lots in common
            fprintf('Neuron %4.2f is not good Right\n',i)
            taste(i).SQ.stats.respFlag = 0;
        end
        
    else
    end
end
%%
clear tasteSelective
x=1;
for i = 1:length(taste)
    if taste(i).MO.stats.respFlag || taste(i).SQ.stats.respFlag
        fprintf('Neuron %4.2f is taste selective \n',i)
        tasteSelective(x) = i;  % index for taste selective neurons
        x = x+1;
    end
end
%% plot the left/right preference vs taste preference
figure;
for i = 1:length(taste)
    taste_pall(i) = max([abs(taste(i).MO.p),abs(taste(i).SQ.p)]); % max taste selectivity
    lr_pall(i)    = data(idx_p(i)).auROC.PlanRcorr_vsLcorr.p;     % direction index
end
h1 = plot(abs(taste_pall),abs(lr_pall),'bo');

for i = 1:length(tasteSelective)
    taste_p(i) = max([abs(taste(tasteSelective(i)).MO.p),abs(taste(tasteSelective(i)).SQ.p)]);
    lr_p(i)    = Sum(idx_p(tasteSelective(i))).auROC.PlanRcorr_vsLcorr.p;
end

hold on
h2 = plot(abs(taste_p),abs(lr_p),'ro');
xlim([0,1])
ylim([0,1])
hold on
plot([0,1],[0,1])
axis square
xlabel('Max taste selectivity')
ylabel('Direction index')
legend([h1,h2],{'p_y < 0.01', 'p_x_,_y < 0.01'})
%% try to show neurons
for i = 1:length(data)
    data(i).Right.Correct = [];
    data(i).Left.Correct = [];
end
