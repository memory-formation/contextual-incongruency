%% Behavioral Analysis
save_root = '';
load([save_root '\Info_matrix_sujes_bons.mat'])

%% Exp 1

%% Subjective Congruence rating
ind_subcong = squeeze(Info_matrix_sujes_bons(:,5,:));

for s= 1:29
mean_rate_score(s,:) = mean(ind_subcong(s,ind_subcong(s,:)~=0));
end 
mean_rate=mean(mean_rate_score);
SD_rate=std(mean_rate_score);

% n_recall
sub_coh = zeros(size(ind_subcong,1),4);
sub_coh(:,1)=sum(ind_subcong==1,2)./size(ind_subcong,2)*100;
sub_coh(:,2)=sum(ind_subcong==2,2)./size(ind_subcong,2)*100;
sub_coh(:,3)=sum(ind_subcong==3,2)./size(ind_subcong,2)*100;
sub_coh(:,4)=sum(ind_subcong==4,2)./size(ind_subcong,2)*100;

mean_per = mean(sub_coh/100,1);
SD_all = std(sub_coh/100,0,1);

%% Memory
% Gist
ind_mem = squeeze(Info_matrix_sujes_bons(:,2,:));

mem_gist = zeros(size(ind_mem,1),2);
mem_gist(:,1)=sum(ind_mem==2,2)./size(ind_mem,2)*100;
mem_gist(:,2)=sum(ind_mem==1,2)./size(ind_mem,2)*100;

mean_per = mean(mem_gist/100,1);
SD_all = std(mem_gist/100,0,1);

% recollection
recle = squeeze(Info_matrix_sujes_bons(:,3,:));

recle_per = zeros(29,2);
for s=1:29
    recle_per(s,1)=sum(recle(s,ind_mem(s,:)==2)==2)/sum(ind_mem(s,:)==2)*100;
    recle_per(s,2)=sum(recle(s,ind_mem(s,:)==2)==1)/sum(ind_mem(s,:)==2)*100;
end 

mean_recle_per = mean(recle_per/100,1);
SD_all = std(recle_per/100,0,1);


% Recle vs non-recle
all_recle_per = zeros(29,2);
for s=1:29
    all_recle_per(s,1)=sum(recle(s,:)==2)/size(recle,2)*100;
    all_recle_per(s,2)=sum(recle(s,:)~=2)/size(recle,2)*100;
end 

mean_all_recle_per = mean(all_recle_per/100,1);
SD_all_recle = std(all_recle_per/100,0,1);

[p,h,stats] = signrank(all_recle_per(:,1),all_recle_per(:,2));

%% Memory separated by Congruence rating
% Gist
acc_gist= zeros(29,2);
for s=1:29
    acc_gist(s,1)=sum(ind_mem(s,ind_subcong(s,:)>2)==2)/sum(ind_subcong(s,:)>2)*100;
    acc_gist(s,2)=sum(ind_mem(s,ind_subcong(s,:)<3&ind_subcong(s,:)>0)==2)/sum(ind_subcong(s,:)<3&ind_subcong(s,:)>0)*100;
end 

[~,p_acc_gist,~,t_acc_gist] = ttest(acc_gist(:,1),acc_gist(:,2));

for fig=2
figure(21);clf;
err = std(acc_gist);
b=bar(mean(acc_gist,1),0.8);hold on
er = errorbar(mean(acc_gist,1),err); 
scatter(ones(size(acc_gist(:,1))).*(1+(rand(size(acc_gist(:,1)))-0.5)/5),acc_gist(:,1),'k','filled');
scatter(ones(size(acc_gist(:,2))).*(2+(rand(size(acc_gist(:,2)))-0.5)/5),acc_gist(:,2),'k','filled');hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;
% b.FaceColor = 'flat';
% b.FaceAlpha = 0.6;
%     b.CData(1,:) = [0,0,1];
%     b.CData(2,:) = [1,0,0];
xticks([1 2 ])
xticklabels({'High','Low'})
ylim([0,100])
ylabel('Percentage of trials')
title('Gist')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end 

% Object item
acc_pic= zeros(29,2);
for s=1:29
    acc_pic(s,1)=sum(recle(s,ind_subcong(s,:)>2 & ind_mem(s,:)==2)==2)/sum(ind_mem(s,:)==2 & ind_subcong(s,:)>2 )*100;
    acc_pic(s,2)=sum(recle(s,ind_subcong(s,:)<3 & ind_subcong(s,:)>0 & ind_mem(s,:)==2)==2)/sum(ind_mem(s,:)==2 & ind_subcong(s,:)<3 & ind_subcong(s,:)>0)*100;
end 

[~,p_acc_pic,~,t_acc_pic] = ttest(acc_pic(:,1),acc_pic(:,2));
for fig=2
figure(22);clf;
err = std(acc_pic);
b=bar(mean(acc_pic,1),0.8);hold on
er = errorbar(mean(acc_pic,1),err); 
scatter(ones(size(acc_pic(:,1))).*(1+(rand(size(acc_pic(:,1)))-0.5)/5),acc_pic(:,1),'k','filled');
scatter(ones(size(acc_pic(:,2))).*(2+(rand(size(acc_pic(:,2)))-0.5)/5),acc_pic(:,2),'k','filled');hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;
% b.FaceColor = 'flat';
% b.FaceAlpha = 0.6;
%     b.CData(1,:) = [0,0,1];
%     b.CData(2,:) = [1,0,0];
xticks([1 2 ])
xticklabels({'High','Low'})
ylim([0,100])
ylabel('Percentage of trials')
title('Item')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end 


%% Exp 2
load([save_root '\Results_Experiment2.mat'])

%% Coh rating 
acc_per_coh = Res_exp2(:,1:2);

mean_per = mean(acc_per_coh,1);
SD_all = std(acc_per_coh,0,1);

acc_coh = acc_per_coh*100;

figure(1);clf;
for fig=1
err = std(acc_coh);
b=bar(mean(acc_coh,1),0.8);hold on
er = errorbar(mean(acc_coh,1),err); 
scatter(ones(size(acc_coh(:,1))).*(1+(rand(size(acc_coh(:,1)))-0.5)/5),acc_coh(:,1),'k','filled');
scatter(ones(size(acc_coh(:,2))).*(2+(rand(size(acc_coh(:,2)))-0.5)/5),acc_coh(:,2),'k','filled');hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;
% b.FaceColor = 'flat';
% b.FaceAlpha = 0.6;
%     b.CData(1,:) = [0,0,1];
%     b.CData(2,:) = [1,0,0];
xticks([1 2 3 4])
xticklabels({'High congruence','Low congruence'})
ylim([0,100])
ylabel('Percentage of trials')
%title('Subjective Coh_rating')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end 
