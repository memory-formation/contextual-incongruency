%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMM on stability for sequential pictures  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_root = '';
load([save_root '/Erps_0.1_20hz.mat'])


%%%%%%%%%%%%%%%
%% Smoothing %%
%%%%%%%%%%%%%%%

N=50; % N of points to be smoothed
for isubj=1:29
    display (['Smoothing Subject Number: ' int2str(isubj) '/' int2str(29)]);
    SIMI_general(isubj).S44 = zeros(size(DATA_ERPs_img(isubj).S44));
    SIMI_general(isubj).offset = zeros(size(Data_offset{isubj}(:,1:1051,:)));


    for ichan=1:29

        for itrial=1:size(DATA_ERPs_img(isubj).S44,3)
            SIMI_general(isubj).S44(ichan,:,itrial)=smoothdata(DATA_ERPs_img(isubj).S44(ichan,:,itrial),'movmean',N);
        end

        for itrial=1:size(Data_offset{isubj},3)
            SIMI_general(isubj).offset(ichan,:,itrial)=smoothdata(Data_offset{isubj}(ichan,1:1051,itrial),'movmean',N);
        end

    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Downsample for Stability analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SIMI_general_down.S44=zeros(0);
for s = 1:29
    
    %%%%%%%% S44 %%%%%%%%%%%
    for i=1:size(SIMI_general(s).S44,3)
        for j=1:29
            SIMI_general_down(s).S44(j,:,i)=downsample(squeeze(SIMI_general(s).S44(j,52:end,i)),5); % downsample and removal of the baseline
            
        end 
    end 
    clearvars i j 

    %%%%%%%% Offset %%%%%%%%%%%
    for i=1:size(SIMI_general(s).S44,3)
        for j=1:29
            SIMI_general_down(s).offset(j,:,i)=downsample(squeeze(SIMI_general(s).offset(j,52:end,i)),5); % downsample and removal of the baseline

        end
    end
    clearvars i j
end






%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:29

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S44 with Offset  %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['Stability analysis for S44-Offset Subject Number: ' int2str(s) '/' int2str(29)]);  
    for i=1:80 
        [R,P] = corr (SIMI_general_down(s).S44(:,:,i),SIMI_general_down(s).offset(:,:,i),'type','pearson');
        SIMI_result_general(s).S44_offset(:,i) = fisherz(diag(R));
    end 

end





%%%%%%%%%%%%%%%%%%%%%%%%
%% LMM %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%% Linear Mixed-effect Model on Stability data

save_root = '';
load([save_root '/reject_index.mat'])
load([save_root '/Info_matrix_sujes_bons.mat'])


%% Prepare behavioral results structure
ind_subcong = squeeze(Info_matrix_sujes_bons(:,5,:));
ind_mem = squeeze(Info_matrix_sujes_bons(:,2,:));
ind_recl = squeeze(Info_matrix_sujes_bons(:,3,:));

ind_S44off2_long = reject_index_S44.*reject_index_offset2_long;


%% Exclude possible NaNs EEG trials

NaN1 = ones (29,80);

for s = 1:29

    ind_S44=ind_S44off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;

    [row_1, col_1] = find(isnan(SIMI_result_general(s).S44_offset(:,:)));

    if length(col_1) > 0; NaN1 (s,col_1(1,:)) = 0; end

end

ind_S44off2_long = reject_index_S44.*reject_index_offset2_long.*NaN1;



Mem = zeros(0); Sub_ID = zeros(0); Cog = zeros(0); Y_all = zeros(0); 


for s = 1:29

    display (['Preparing LMM for Subject Number: ' int2str(s) '/' int2str(29)]);


    ind_S44=ind_S44off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & ind_subcong(s,:)~=0;
    n_S44 = sum(ind_S44);

    %% sub ID
    Sub_ID = cat(2,Sub_ID,rep_num(s,1,n_S44));

    %% Mem
    Mem_sub = [ind_mem(s,ind_S44)];
    Mem = cat(2,Mem,Mem_sub); clearvars Mem_sub

    %% cog rating
    Cog_sub = [ind_subcong(s,ind_S44)];
    Cog = cat(2,Cog,Cog_sub); clearvars Cog_sub

    %% SIMI data
    SIMI_S44 = SIMI_result_general(s).S44_offset(:,ind_S44);

    SIMI_sub = cat(2,SIMI_S44);
    Y_all = cat(2,Y_all,SIMI_sub);
end
 


Cog(:)=zscore(Cog(:));


%% LMM on scalp ERP  

Ydata_alltp_rs = Y_all;

t_map_real = zeros(3,size(Ydata_alltp_rs,1)); %1,intercept; 
beta_map_real = zeros(3,size(Ydata_alltp_rs,1)); %1,intercept; 
p_map_real = zeros(3,size(Ydata_alltp_rs,1)); %1,intercept; 


for i_time = 1:size(Ydata_alltp_rs,1)

    display (['LMM for time point: ' int2str(i_time) '/' int2str(size(Ydata_alltp_rs,1))]);


    tbl_time = table(Ydata_alltp_rs(i_time,:)',nominal(Sub_ID(:)),Mem(:),Cog(:),'VariableNames',{'SIMI','Subject_ID','Mem','Coh_rating'});
    lme_time = fitlme(tbl_time,'SIMI~Mem+Coh_rating+(1|Subject_ID)');

    t_map_real(1,i_time) = lme_time.Coefficients{1,4};
    t_map_real(2,i_time) = lme_time.Coefficients{2,4};
    t_map_real(3,i_time) = lme_time.Coefficients{3,4};

    beta_map_real(1,i_time) = lme_time.Coefficients{1,2};
    beta_map_real(2,i_time) = lme_time.Coefficients{2,2};
    beta_map_real(3,i_time) = lme_time.Coefficients{3,2};

    CILow_map_real(1,i_time) = lme_time.Coefficients{1,7};
    CILow_map_real(2,i_time) = lme_time.Coefficients{2,7};
    CILow_map_real(3,i_time) = lme_time.Coefficients{3,7};

    CIHigh_map_real(1,i_time) = lme_time.Coefficients{1,8};
    CIHigh_map_real(2,i_time) = lme_time.Coefficients{2,8};
    CIHigh_map_real(3,i_time) = lme_time.Coefficients{3,8};

    p_map_real(1,i_time) = lme_time.Coefficients{1,6};
    p_map_real(2,i_time) = lme_time.Coefficients{2,6};
    p_map_real(3,i_time) = lme_time.Coefficients{3,6};

    clearvars tbl_time lme_time
end


%% Figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Stability in Figure 3 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SIMI_S44

for s = 1:29

    ind_Cong = ind_S44off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & (ind_subcong(s,:)== 3 | ind_subcong(s,:)== 4);
    ind_Incong = ind_S44off2_long(s,:)==1 & ind_mem(s,:) ~= 0 & (ind_subcong(s,:)== 1 | ind_subcong(s,:)== 2);

    %% Stability data
    SIMI_Cong (s,:) = mean(SIMI_result_general(s).S44_offset(:,ind_Cong),2);
    SIMI_Incong (s,:) = mean(SIMI_result_general(s).S44_offset(:,ind_Incong),2);

    %% Number of trials
    NTr_total(s,1) = sum(ind_Cong);
    NTr_total(s,2) = sum(ind_Incong);
end




