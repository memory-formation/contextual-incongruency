%% Representational Similarity Analysis (RSA)
save_root = '';
load([save_root '\Erps_0.1_20hz.mat'])
load([save_root '\reject_index.mat'])
%% Smoothing
N=50; % N of points to be smoothed
for isubj=1:29
    display (['Smothing Subject Number: ' int2str(isubj) '/' int2str(29)]);
    SIMI_general(isubj).S11=zeros(size(DATA_ERPs_img(isubj).S11));
    SIMI_general(isubj).S22=zeros(size(DATA_ERPs_img(isubj).S22));
    SIMI_general(isubj).S33=zeros(size(DATA_ERPs_img(isubj).S33));
    SIMI_general(isubj).S44=zeros(size(DATA_ERPs_img(isubj).S44));
    SIMI_general(isubj).offset_2_long=zeros(size(Data_offset{isubj}));
    
    for ichan=1:29
        for itrial=1:size(DATA_ERPs_img(isubj).S11,3)
            SIMI_general(isubj).S11(ichan,:,itrial)=smoothdata(DATA_ERPs_img(isubj).S11(ichan,:,itrial),'movmean',N);
        end
        
        for itrial=1:size(DATA_ERPs_img(isubj).S22,3)
            SIMI_general(isubj).S22(ichan,:,itrial)=smoothdata(DATA_ERPs_img(isubj).S22(ichan,:,itrial),'movmean',N);
        end
        
        for itrial=1:size(DATA_ERPs_img(isubj).S33,3)
            SIMI_general(isubj).S33(ichan,:,itrial)=smoothdata(DATA_ERPs_img(isubj).S33(ichan,:,itrial),'movmean',N);
        end
        
        for itrial=1:size(DATA_ERPs_img(isubj).S44,3)
            SIMI_general(isubj).S44(ichan,:,itrial)=smoothdata(DATA_ERPs_img(isubj).S44(ichan,:,itrial),'movmean',N);
        end
        
        for itrial=1:size(Data_offset{isubj},3)
            SIMI_general(isubj).offset_2_long(ichan,:,itrial)=smoothdata(Data_offset{isubj}(ichan,:,itrial),'movmean',N);
        end
       
    end
end 
%% Downsample for Similarity analysis

SIMI_general_down.S11=zeros(0);

for s=1:29
    
    %%%%%%%% S11 %%%%%%%%%%%
    for i=1:size(SIMI_general(s).S11,3)
        for j=1:29
            SIMI_general_down(s).S11(j,:,i)=downsample(squeeze(SIMI_general(s).S11(j,52:end,i)),5); % downsample and removal of the baseline
            
        end 
    end 
    clearvars i j 
    
    %%%%%%%% S22 %%%%%%%%%%%
    for i=1:size(SIMI_general(s).S22,3)
        for j=1:29
            SIMI_general_down(s).S22(j,:,i)=downsample(squeeze(SIMI_general(s).S22(j,52:end,i)),5); % downsample and removal of the baseline
            
        end 
    end 
    clearvars i j 
    %%%%%%%% S33 %%%%%%%%%%%
    for i=1:size(SIMI_general(s).S33,3)
        for j=1:29
            SIMI_general_down(s).S33(j,:,i)=downsample(squeeze(SIMI_general(s).S33(j,52:end,i)),5); % downsample and removal of the baseline
            
        end 
    end 
    clearvars i j 
    %%%%%%%% S44 %%%%%%%%%%%
    for i=1:size(SIMI_general(s).S44,3)
        for j=1:29
            SIMI_general_down(s).S44(j,:,i)=downsample(squeeze(SIMI_general(s).S44(j,52:end,i)),5); % downsample and removal of the baseline
            
        end 
    end 
    clearvars i j 
    

    %%%%%%%% offset %%%%%%%%
    for i=1:size(SIMI_general(s).offset_2_long,3)
        for j=1:29
            SIMI_general_down(s).offset_2_long(j,:,i)=downsample(squeeze(SIMI_general(s).offset_2_long(j,52:end,i)),5); % downsample and removal of the baseline
            
        end 
    end 
end 

clearvars -except SIMI_general_down reject* save_root

%% RSA

for s=1:29
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S11 with offset %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['S11 Subject Number: ' int2str(s) '/' int2str(29)]);  
    for i=1:80
        
        A1=corr_matrix(SIMI_general_down(s).offset_2_long(:,1:200,i),SIMI_general_down(s).S11(:,:,i));

        A2=corr_matrix(SIMI_general_down(s).offset_2_long(:,51:250,i),SIMI_general_down(s).S11(:,:,i));
        %  combine the matrix
        SIMI_result_general(s).S11_offset_2_long(:,:,i) = [A1 A2(:,151:200)];

    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S22 with offset %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['S22 Subject Number: ' int2str(s) '/' int2str(29)]);  
    for i=1:80
        
        A1=corr_matrix(SIMI_general_down(s).offset_2_long(:,1:200,i),SIMI_general_down(s).S22(:,:,i));

        A2=corr_matrix(SIMI_general_down(s).offset_2_long(:,51:250,i),SIMI_general_down(s).S22(:,:,i));
        %  combine the matrix
        SIMI_result_general(s).S22_offset_2_long(:,:,i) = [A1 A2(:,151:200)];

    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S33 with offset %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['S33 Subject Number: ' int2str(s) '/' int2str(29)]);  
    for i=1:80
        
        A1=corr_matrix(SIMI_general_down(s).offset_2_long(:,1:200,i),SIMI_general_down(s).S33(:,:,i));

        A2=corr_matrix(SIMI_general_down(s).offset_2_long(:,51:250,i),SIMI_general_down(s).S33(:,:,i));
        %  combine the matrix
        SIMI_result_general(s).S33_offset_2_long(:,:,i) = [A1 A2(:,151:200)];

    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S44 with offset %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['S44 Subject Number: ' int2str(s) '/' int2str(29)]);  
    for i=1:80
        
        A1=corr_matrix(SIMI_general_down(s).offset_2_long(:,1:200,i),SIMI_general_down(s).S44(:,:,i));

        A2=corr_matrix(SIMI_general_down(s).offset_2_long(:,51:250,i),SIMI_general_down(s).S44(:,:,i));
        %  combine the matrix
        SIMI_result_general(s).S44_offset_2_long(:,:,i) = [A1 A2(:,151:200)];

    end 
end 



%% Linear Mixed-effect Model on Representational Similarity Analysis Results
% load Similarity results
save_root = '';
load([save_root '\reject_index.mat'])
load([save_root '\Info_matrix_sujes_bons.mat'])

%% Prepare behavioral results structure
ind_subcong = squeeze(Info_matrix_sujes_bons(:,5,:));
ind_mem = squeeze(Info_matrix_sujes_bons(:,2,:));
ind_recl = squeeze(Info_matrix_sujes_bons(:,3,:));

ind_S11off2_long = reject_index_S11.*reject_index_offset2_long;
ind_S22off2_long = reject_index_S22.*reject_index_offset2_long;
ind_S33off2_long = reject_index_S33.*reject_index_offset2_long;
ind_S44off2_long = reject_index_S44.*reject_index_offset2_long;

Mem = zeros(0);Sub_ID = zeros(0);Seq = zeros(0); Cog = zeros(0);Y_all = zeros(0);Recl=zeros(0);
for s=1:29
    disp(s)
    ind_S11=ind_S11off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
    ind_S22=ind_S22off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
    ind_S33=ind_S33off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
    ind_S44=ind_S44off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
    
    n_S11 = sum(ind_S11);n_S22 = sum(ind_S22);n_S33 = sum(ind_S33);n_S44 = sum(ind_S44);
    %% sub ID
    Sub_ID = cat(2,Sub_ID,rep_num(s,1,n_S11+n_S22+n_S33+n_S44));
    %% Mem
    Mem_sub = [ind_mem(s,ind_S11),ind_mem(s,ind_S22),ind_mem(s,ind_S33),ind_mem(s,ind_S44)];
    Mem = cat(2,Mem,Mem_sub); clearvars Mem_sub
    %% Recollection
    Recl_sub = [ind_recl(s,ind_S11),ind_recl(s,ind_S22),ind_recl(s,ind_S33),ind_recl(s,ind_S44)];
    Recl = cat(2,Recl,Recl_sub); clearvars Recl_sub
    %% Seq
    Seq_sub = [rep_num(1,1,n_S11),rep_num(2,1,n_S22),rep_num(3,1,n_S33),rep_num(4,1,n_S44)];
    Seq = cat(2,Seq,Seq_sub); clearvars Seq_sub
    %% cog rating
    Cog_sub = [ind_subcong(s,ind_S11),ind_subcong(s,ind_S22),ind_subcong(s,ind_S33),ind_subcong(s,ind_S44)];
    Cog = cat(2,Cog,Cog_sub); clearvars Cog_sub
    %% SIMI data
    SIMI_S11 = SIMI_result_general(s).S11_offset_2_long(:,:,ind_S11);
    SIMI_S22 = SIMI_result_general(s).S22_offset_2_long(:,:,ind_S22);
    SIMI_S33 = SIMI_result_general(s).S33_offset_2_long(:,:,ind_S33);
    SIMI_S44 = SIMI_result_general(s).S44_offset_2_long(:,:,ind_S44);
    
    SIMI_sub = cat(3,SIMI_S11,SIMI_S22,SIMI_S33,SIMI_S44);
    Y_all = cat(2,Y_all,reshape(SIMI_sub,[size(SIMI_sub,1)*size(SIMI_sub,2),size(SIMI_sub,3)]));
    
end 


%% Prepare RSA result
Simi_all = reshape(Y_all,[200,250,size(Y_all,2)]);

% select the key time window
[M_offset,I_offset]=max(mean(Simi_all,3),[],1);

[M,I]=max(mean(Simi_all,3),[],2);

[ps,I_enc]=max(M); % timpepoint with the highest Simi value
[ps_off,I_off]=max(M_offset); % timpepoint with the highest Simi value

tw = I_enc-5:I_enc+5; % [-50,50] of peak

SIMI_1D = squeeze(mean(Simi_all(tw,:,:),1));

%% Normalize variables

for i_time = 1:250
    SIMI_1D_z(i_time,:) = zscore(SIMI_1D(i_time,:));
    
end 
Seq(:) = zscore(Seq(:));Cog(:)=zscore(Cog(:));

%% LMM on RSA (sequence with object encoding)
t_map_real = zeros(4,250); 
beta_map_real = zeros(4,250); 
p_map_real = zeros(4,250); 
CI_low_map_real = zeros(4,250); 
CI_high_map_real = zeros(4,250); 
Ydata_alltp_rs = SIMI_1D_z;

for i_time =1:size(Ydata_alltp_rs,1)
    if floor(i_time/40)==ceil(i_time/40)
        disp([int2str(i_time/(size(Ydata_alltp_rs,1))*100) '% complete'])
    end
    tbl_time = table(Ydata_alltp_rs(i_time,:)',Recl(:)==2,Sub_ID(:),Seq(:),Cog(:),'VariableNames',{'SIMI','Recl','Subject_ID','Seq','Coh_rating'});
    lme_time = fitlme(tbl_time,'SIMI~Recl+Seq+Coh_rating+(1|Subject_ID)');

    
    t_map_real(1,i_time) = lme_time.Coefficients{1,4};
    t_map_real(2,i_time) = lme_time.Coefficients{2,4};
    t_map_real(3,i_time) = lme_time.Coefficients{3,4};
    t_map_real(4,i_time) = lme_time.Coefficients{4,4};
    
    
    beta_map_real(1,i_time) = lme_time.Coefficients{1,2};
    beta_map_real(2,i_time) = lme_time.Coefficients{2,2};
    beta_map_real(3,i_time) = lme_time.Coefficients{3,2};
    beta_map_real(4,i_time) = lme_time.Coefficients{4,2};
    
    
    p_map_real(1,i_time) = lme_time.Coefficients{1,6};
    p_map_real(2,i_time) = lme_time.Coefficients{2,6};
    p_map_real(3,i_time) = lme_time.Coefficients{3,6};
    p_map_real(4,i_time) = lme_time.Coefficients{4,6};
    
    CI_low_map_real(1,i_time) = lme_time.Coefficients{1,7};
    CI_low_map_real(2,i_time) = lme_time.Coefficients{2,7};
    CI_low_map_real(3,i_time) = lme_time.Coefficients{3,7};
    CI_low_map_real(4,i_time) = lme_time.Coefficients{4,7};
    
    CI_high_map_real(1,i_time) = lme_time.Coefficients{1,8};
    CI_high_map_real(2,i_time) = lme_time.Coefficients{2,8};
    CI_high_map_real(3,i_time) = lme_time.Coefficients{3,8};
    CI_high_map_real(4,i_time) = lme_time.Coefficients{4,8};
    
    clearvars tbl_time lme_time
end


%% LMM on RSA (sequence with episodic sequence offset)

%% Ploting
% Figure 6 a&b
%% plot
for fig=6
    %% RSA (Figure 6a)
    figure(51);clf;
    imagesc(mean(Simi_all,3))
    caxis([-0.4,0.4])
    xticks([ 1 50 100 150 200 250])
    xticklabels({'0','500','1000','1500','2000','2500'})
    yticks([ 1 50 100 150 200])
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('Object encoding [ms]')
    ylabel('Sequence encoding [ms]')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
    
    %% LMM Memory (Figure 6b, upper)

    figure(64);clf;
    plot(beta_map_real(2,:),'linewidth',3);hold on 
    plot(CI_low_map_real(2,:),'linewidth',2,'color','black');
    plot(CI_high_map_real(2,:),'linewidth',2,'color','black');  
    hold off
    ylim([-0.25,0.1])
    xticks([0 0.5 1 1.5 2 2.5 ]*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Beta')
    title('Memory (LMM)')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(641);clf;
    imagesc(t_map_real(2,:))
    colormap('hot')
    caxis([-5,2])
    
    %% LMM Congruence (Figure 6b, middle)
    
    figure(65);clf;
    plot(beta_map_real(4,:),'linewidth',3);hold on 
    plot(CI_low_map_real(4,:),'linewidth',2,'color','black');
    plot(CI_high_map_real(4,:),'linewidth',2,'color','black');  
    hold off
    ylim([-0.1,0.06])
    xticks([0 0.5 1 1.5 2 2.5 ]*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Beta')
    title('Congruence (LMM)')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(651);clf;
    imagesc(t_map_real(4,:))
    colormap('hot')
    caxis([-5,2])
    
    %% LMM Sequence order (Figure 6b, lower)
    
    figure(65);clf;
    plot(beta_map_real(3,:),'linewidth',3);hold on 
    plot(CI_low_map_real(3,:),'linewidth',2,'color','black');
    plot(CI_high_map_real(3,:),'linewidth',2,'color','black');  
    hold off
    ylim([-0.1,0.06])
    xticks([0 0.5 1 1.5 2 2.5 ]*100)
    xticklabels({'0','500','1000','1500','2000','2500'})
    xlabel('Offset [ms]')
    ylabel('Beta')
    title('Sequence order (LMM)')
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(651);clf;
    imagesc(t_map_real(3,:))
    colormap('hot')
    caxis([-5,2])
end

