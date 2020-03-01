%% sort the preparatory activity with the [-1,0] range
function auROC_colormap_decision(time,PSTH_auROC_R, PSTH_auROC_L)
iddx = find(time> -2 & time <0);
reduce_R = PSTH_auROC_R(:,iddx);
for i = 1:size(reduce_R,1)
    len_th(i) = length(find(reduce_R(i,:) > 0.1));
end

[B, I] = sort(len_th,'descend');
clear len_th
for i = 1:size(reduce_R,1)
    PSTH_auROC_R_reorg(i,:) = PSTH_auROC_R(I(i),:);
end
% figure;
% imagesc(time,[],PSTH_auROC_R_reorg)

reduce_L = PSTH_auROC_L(:,iddx);
for i = 1:size(reduce_L,1)
    len_th(i) = length(find(reduce_L(i,:) <- 0.1));
end

[B, I] = sort(len_th,'descend');

for i = 1:size(reduce_L,1)
    PSTH_auROC_L_reorg(i,:) = PSTH_auROC_L(I(i),:);
end
% figure;
% imagesc(time,[],PSTH_auROC_L_reorg)
figure;
imagesc(time,[],[PSTH_auROC_R_reorg;PSTH_auROC_L_reorg])
xlim([-2,1])
caxis([-0.4,0.4])
colormap(jet)