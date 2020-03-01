%% This function is used for bar plots for mutiple classes
% data should be in matrix
function bar_plot_multi(X,stats)
if nargin < 2
    stats = false;
end
val_m = mean(X,1);
val_sem = std(X,1)./sqrt(size(X,1));
figure;
color = {'c','m','r','y','b','k','g'};
hold on
for i = 1:length(val_m)
    bar(i,val_m(i),'FaceColor',color{i},'EdgeColor',[0,0,0],'LineWidth',1)
end
e=errorbar(val_m, val_sem,'LineWidth',1.5);
ylim([0,2*max(val_m)])
e.Color ='k';
e.LineStyle = 'none';
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.01 0.01]);
if stats
    p = anova1(X);
    fprintf('p value is %f \n',p)
end
end

