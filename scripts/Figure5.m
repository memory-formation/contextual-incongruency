%% ERPs analysis
save_root = '';
load([save_root '\Erps_0.1_20hz.mat'])
load([save_root '\reject_index.mat'])


%% Obejct encoding ERPs
for s=1:29
    ERP_data{s}=Data_offset{s};
    
end 
clearvars -except ERP_data reject* save_root

%% Smoothing
t_select = 52:1301;

N=50; % N of points to be smoothed
for isubj=1:29
    display (['Smothing Subject Number: ' int2str(isubj) '/' int2str(29)]);
    for ichan=1:29
       
        for itrial=1:80
            Data_general(isubj).offset_2_long(ichan,:,itrial)=smoothdata(ERP_data{isubj}(ichan,t_select,itrial),'movmean',N);
        end   
    end
end 
%% Downsample for Similarity analysis
Data_general_down.offset_2_long=zeros(0);
for s=1:29
    
   
    %%%%%%%% offset2 long %%%%%%%%
    for i=1:80
        for j=1:29
            offset_2_long(j,:,i)=downsample(squeeze(Data_general(s).offset_2_long(j,:,i)),5); % downsample and removal of the baseline
            
        end 
    end 
    ERP_ds{s}=offset_2_long;clearvars offset_2_long
end 

clearvars -except ERP_ds reject* save_root

%% Prepare Behavioral results structure
load([save_root '\Info_matrix_sujes_bons.mat'])

ind_subcong = squeeze(Info_matrix_sujes_bons(:,5,:));
ind_mem = squeeze(Info_matrix_sujes_bons(:,2,:));
recle = squeeze(Info_matrix_sujes_bons(:,3,:));

% beh results
Mem = zeros(0);Sub_ID = zeros(0); Cog = zeros(0);Recl=zeros(0);
for s=1:29
    ind=reject_index_offset2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
    
    n_trial = sum(ind);
    %% sub ID
    Sub_ID = cat(2,Sub_ID,rep_num(s,1,n_trial));
    %% Mem
    Mem = cat(2,Mem,ind_mem(s,ind));
    %% Recollection
    Recl = cat(2,Recl,recle(s,ind)); 
    %% cog rating
    Cog = cat(2,Cog,ind_subcong(s,ind)); 
end 


%% ERPs average within eahc region
Y_all = zeros(0);
for ichan=1:29
    disp(ichan)
    Data_ERP_chan=zeros(0);
    for s=1:29
        ind=reject_index_offset2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
        
        Data_ERP_chan = cat(2,Data_ERP_chan,squeeze(ERP_ds{s}(ichan,:,ind)));
    end
    Y_all(ichan,:,:)=Data_ERP_chan; clearvars Data_ERP_chan
end

channel_list = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','FCz','Oz','F7','F8',...
    'T3','T4','T5','T6','Fz','Cz','Pz','Fc1','Fc2','Fc5','Fc6','Cp1','Cp2',...
    'Cp5','Cp6','Po1','Po2'};
Zone_1 = [1 3 9 11 17 20 22];
Zone_2 = [2 4 9 12 17 21 23];
Zone_3 = [5 9 13 18 20 22 24 26];
Zone_4 = [6 9 14 18 21 23 25 27];
Zone_5 = [7 10 15 19 24 26 28];
Zone_6 = [8 10 16 19 25 27 29];


Y_all_reginal=zeros([6,size(Y_all,[2,3])]);
Y_all_reginal(1,:,:)= mean(Y_all(Zone_1,:,:),1);
Y_all_reginal(2,:,:)= mean(Y_all(Zone_2,:,:),1);
Y_all_reginal(3,:,:)= mean(Y_all(Zone_3,:,:),1);
Y_all_reginal(4,:,:)= mean(Y_all(Zone_4,:,:),1);
Y_all_reginal(5,:,:)= mean(Y_all(Zone_5,:,:),1);
Y_all_reginal(6,:,:)= mean(Y_all(Zone_6,:,:),1);

%% Normalize the data before next step

for i_reg = 1:6
    for i_time = 1:250
      Y_all_reginal_z(i_reg,i_time,:) = zscore(Y_all_reginal(i_reg,i_time,:));
    end 
end 
Cog(:)=zscore(Cog(:));    
%% LMM on Offset ERPs for each region   
t_map_real = zeros(3,250,6); %1,intercept; 
beta_map_real = zeros(3,250,6); %1,intercept; 
p_map_real = zeros(3,250,6); %1,intercept; 
CI_low_map_real = zeros(3,250,6); 
CI_high_map_real = zeros(3,250,6);

for i_reg=1:6
    disp([ 'Region',int2str(i_reg)])
    
    Ydata_alltp_rs = squeeze(Y_all_reginal_z(i_reg,:,:));
    for i_time =1:size(Ydata_alltp_rs,1)

        tbl_time = table(Ydata_alltp_rs(i_time,:)',Recl(:)==2,nominal(Sub_ID(:)),Cog(:),'VariableNames',{'ERP','Recl','Subject_ID','Coh_rating'});
        lme_time = fitlme(tbl_time,'ERP~Recl+Coh_rating+(1|Subject_ID)');
            
        t_map_real(1,i_time,i_reg) = lme_time.Coefficients{1,4};
        t_map_real(2,i_time,i_reg) = lme_time.Coefficients{2,4};
        t_map_real(3,i_time,i_reg) = lme_time.Coefficients{3,4};
        
        beta_map_real(1,i_time,i_reg) = lme_time.Coefficients{1,2};
        beta_map_real(2,i_time,i_reg) = lme_time.Coefficients{2,2};
        beta_map_real(3,i_time,i_reg) = lme_time.Coefficients{3,2};

        p_map_real(1,i_time,i_reg) = lme_time.Coefficients{1,6};
        p_map_real(2,i_time,i_reg) = lme_time.Coefficients{2,6};
        p_map_real(3,i_time,i_reg) = lme_time.Coefficients{3,6};
        
        CI_low_map_real(1,i_time,i_reg) = lme_time.Coefficients{1,7};
        CI_low_map_real(2,i_time,i_reg) = lme_time.Coefficients{2,7};
        CI_low_map_real(3,i_time,i_reg) = lme_time.Coefficients{3,7};
        
        CI_high_map_real(1,i_time,i_reg) = lme_time.Coefficients{1,8};
        CI_high_map_real(2,i_time,i_reg) = lme_time.Coefficients{2,8};
        CI_high_map_real(3,i_time,i_reg) = lme_time.Coefficients{3,8};
 
        
        clearvars tbl_time lme_time
    end
    clearvars Ydata_alltp_rs
end

%% Plot
for fig=5
    figure(501);clf;
    imagesc(squeeze(t_map_real(2,:,:))');
    colormap('hot')
    caxis([-5,5])
    xticks(([0.1 0.6 1.1 1.6 2.1 2.6 ]-0.1)*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Region')
    yticks([1,2,3,4,5,6])
    yticklabels({'1','2','3','4','5','6'})
    title('Memory')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(502);clf;
    imagesc(squeeze(p_map_real(2,:,:))'<0.05/(250*6));hold on
    B = bwboundaries(squeeze(p_map_real(2,:,:))'<0.05/(250*6));
    for k=1:size(B,1)
        b = B{k};    plot(b(:,2),b(:,1),'k','LineWidth',3);
    end
    hold off
    xticks(([0.1 0.6 1.1 1.6 2.1 2.6 ]-0.1)*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Region')
    yticks([1,2,3,4,5,6])
    yticklabels({'1','2','3','4','5','6'})
    title('Memory (Sig)');
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(503);clf;
    imagesc(squeeze(beta_map_real(2,:,:))');
    caxis([-0.1,0.2])
    xticks(([0.1 0.6 1.1 1.6 2.1 2.6 ]-0.1)*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Region')
    yticks([1,2,3,4,5,6])
    yticklabels({'1','2','3','4','5','6'})
    title('Memory beta estimate')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    
    %
    figure(601);clf;
    imagesc(squeeze(t_map_real(3,:,:))');
    colormap('hot')
    caxis([-5,5])
    xticks(([0.1 0.6 1.1 1.6 2.1 2.6 ]-0.1)*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Region')
    yticks([1,2,3,4,5,6])
    yticklabels({'1','2','3','4','5','6'})
    title('Coherence')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(602);clf;
    imagesc(squeeze(p_map_real(3,:,:))'<0.05/(250*6));hold on
    B = bwboundaries(squeeze(p_map_real(3,:,:))'<0.05/(250*6));
    for k=1:size(B,1)
        b = B{k};    plot(b(:,2),b(:,1),'k','LineWidth',3);
    end
    hold off
    xticks(([0.1 0.6 1.1 1.6 2.1 2.6 ]-0.1)*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Region')
    yticks([1,2,3,4,5,6])
    yticklabels({'1','2','3','4','5','6'})
    title('Coherence (Sig)');
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(603);clf;
    imagesc(squeeze(beta_map_real(3,:,:))');
    caxis([-0.1,0.1])
    xticks(([0.1 0.6 1.1 1.6 2.1 2.6 ]-0.1)*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Region')
    yticks([1,2,3,4,5,6])
    yticklabels({'1','2','3','4','5','6'})
    title('Congruence beta estimate')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    %% raw EPRs
    ind_subcong = squeeze(Info_matrix_sujes_bons(:,5,:));
    ind_mem = squeeze(Info_matrix_sujes_bons(:,2,:));
    recle = squeeze(Info_matrix_sujes_bons(:,3,:));
    
    % beh results
    Mem = zeros(0);Sub_ID = zeros(0); Cog = zeros(0);Recl=zeros(0);
    for s=1:29
        ind=reject_index_offset2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
        
        n_trial = sum(ind);
        %% sub ID
        Sub_ID = cat(2,Sub_ID,rep_num(s,1,n_trial));
        %% Mem
        Mem = cat(2,Mem,ind_mem(s,ind));
        %% Recollection
        Recl = cat(2,Recl,recle(s,ind));
        %% cog rating
        Cog = cat(2,Cog,ind_subcong(s,ind));
    end

    for s=1:29
        ERP_coghigh_zone4(s,:) = squeeze(mean(mean(Y_all(Zone_4,:,Sub_ID==s & Cog>2),3),1));
        ERP_coglow_zone4(s,:) = squeeze(mean(mean(Y_all(Zone_4,:,Sub_ID==s & Cog<3 & Cog>0),3),1));
        
        
        ERP_recle_zone1(s,:) = squeeze(mean(mean(Y_all(Zone_1,:,Sub_ID==s & Recl==2),3),1));
        ERP_nonrecle_zone1(s,:) = squeeze(mean(mean(Y_all(Zone_1,:,Sub_ID==s & Recl~=2),3),1));
        
    end
    
    % Coh zone4
    for subfig=93
        options.handle     = figure(subfig);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        options.alpha      = 0.7;
        options.line_width = 2;
        options.error      = 'sem';
        
        plot_areaerrorbar(ERP_coghigh_zone4,options);hold on
        
        options.color_area = [219 193 128]./255;    %
        options.color_line = [ 186 148 52]./255;
        
        plot_areaerrorbar(ERP_coglow_zone4,options); hold off
        
        xticks([0 0.5 1 1.5 2 2.5]*100)
        xticklabels({'0','500','1000','1500','2000','2500'})
        xlabel('Offset [ms]')
        ylabel('mV')
        title('Zone 4')
        legend({'High','','Low Coh',' '})
        ylim([-10,10])
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    end
    
    % Recle zone1
    for subfig=94
        options.handle     = figure(subfig);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        options.alpha      = 0.7;
        options.line_width = 2;
        options.error      = 'sem';
        
        plot_areaerrorbar(ERP_recle_zone1,options);hold on
        
        options.color_area = [219 193 128]./255;    %
        options.color_line = [ 186 148 52]./255;
        
        plot_areaerrorbar(ERP_nonrecle_zone1,options); hold off
        
        xticks([0 0.5 1 1.5 2 2.5]*100)
        xticklabels({'0','500','1000','1500','2000','2500'})
        xlabel('Offset [ms]')
        ylabel('mV')
        title('Zone 1')
        legend({'Recle','',' Non Recle',' '})
        ylim([-10,10])
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    end


end
