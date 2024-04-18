


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
    SIMI_general(isubj).S11=zeros(size(DATA_ERPs_img(isubj).S11));
    SIMI_general(isubj).S22=zeros(size(DATA_ERPs_img(isubj).S22));
    SIMI_general(isubj).S33=zeros(size(DATA_ERPs_img(isubj).S33));
    SIMI_general(isubj).S44=zeros(size(DATA_ERPs_img(isubj).S44));



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

    end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Downsample for Stability analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
end 









%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:29

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S11 with S22   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['Stability analysis for S11-S22 Subject Number: ' int2str(s) '/' int2str(29)]);  
    for i=1:80 
        [R,P] = corr (SIMI_general_down(s).S11(:,:,i),SIMI_general_down(s).S22(:,:,i),'type','pearson');
        SIMI_result_general(s).S11_S22(:,i) = fisherz(diag(R));
    end 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S22 with S33   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['Stability analysis for S22-S33 Subject Number: ' int2str(s) '/' int2str(29)]);  
    for i=1:80 
        [R,P] = corr (SIMI_general_down(s).S22(:,:,i),SIMI_general_down(s).S33(:,:,i),'type','pearson');
        SIMI_result_general(s).S22_S33(:,i) = fisherz(diag(R));
    end 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S33 with S44   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display (['Stability analysis for S33-S44 Subject Number: ' int2str(s) '/' int2str(29)]);
    for i=1:80
        [R,P] = corr (SIMI_general_down(s).S33(:,:,i),SIMI_general_down(s).S44(:,:,i),'type','pearson');
        SIMI_result_general(s).S33_S44(:,i) = fisherz(diag(R));
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

ind_S1S2 = reject_index_S11.*reject_index_S22;
ind_S2S3 = reject_index_S22.*reject_index_S33;
ind_S3S4 = reject_index_S33.*reject_index_S44;

%% Exclude possible NaNs EEG trials

NaN1 = ones (29,80);
NaN2 = ones (29,80);
NaN3 = ones (29,80);

for s = 1:29

    ind_S11=ind_S1S2(s,:)==1 & ind_subcong(s,:)~=0;
    ind_S22=ind_S2S3(s,:)==1 & ind_subcong(s,:)~=0;
    ind_S33=ind_S3S4(s,:)==1 & ind_subcong(s,:)~=0;

    [row_1, col_1] = find(isnan(SIMI_result_general(s).S11_S22(:,:)));
    [row_2, col_2] = find(isnan(SIMI_result_general(s).S22_S33(:,:)));
    [row_3, col_3] = find(isnan(SIMI_result_general(s).S33_S44(:,:)));


    if length(col_1) > 0; NaN1 (s,col_1(1,:)) = 0; end
    if length(col_2) > 0; NaN2 (s,col_2(1,:)) = 0; end
    if length(col_3) > 0; NaN3 (s,col_3(1,:)) = 0; end
end

ind_S1S2 = reject_index_S11.*reject_index_S22.*NaN1;
ind_S2S3 = reject_index_S22.*reject_index_S33.*NaN2;
ind_S3S4 = reject_index_S33.*reject_index_S44.*NaN3;


Sub_ID = zeros(0); Seq = zeros(0); Y_all = zeros(0); 


for s = 1:29

    display (['Preparing LMM for Subject Number: ' int2str(s) '/' int2str(29)]);

    ind_S11=ind_S1S2(s,:)==1 & ind_subcong(s,:)~=0;
    ind_S22=ind_S2S3(s,:)==1 & ind_subcong(s,:)~=0;
    ind_S33=ind_S3S4(s,:)==1 & ind_subcong(s,:)~=0;


    n_S11 = sum(ind_S11);n_S22 = sum(ind_S22);n_S33 = sum(ind_S33);

    %% sub ID
    Sub_ID = cat(2,Sub_ID,rep_num(s,1,n_S11+n_S22+n_S33));


    %% Seq
    Seq_sub = [rep_num(1,1,n_S11),rep_num(2,1,n_S22),rep_num(3,1,n_S33)];
    Seq = cat(2,Seq,Seq_sub); clearvars Seq_sub


    %% Stability data
    SIMI_S11_S22 = SIMI_result_general(s).S11_S22(:,ind_S11);
    SIMI_S22_S33 = SIMI_result_general(s).S22_S33(:,ind_S22);
    SIMI_S33_S44 = SIMI_result_general(s).S33_S44(:,ind_S33);


    SIMI_sub = cat(2,SIMI_S11_S22,SIMI_S22_S33,SIMI_S33_S44);
    Y_all = cat(2,Y_all,SIMI_sub);
end


Seq(:) = zscore(Seq(:));

%% LMM on scalp ERP  

Ydata_alltp_rs = Y_all;

t_map_real = zeros(2,size(Ydata_alltp_rs,1)); %1,intercept; 
beta_map_real = zeros(2,size(Ydata_alltp_rs,1)); %1,intercept; 
p_map_real = zeros(1,size(Ydata_alltp_rs,1)); %1,intercept; 


for i_time = 1:size(Ydata_alltp_rs,1)

    display (['LMM for time point: ' int2str(i_time) '/' int2str(size(Ydata_alltp_rs,1))]);

    tbl_time = table(Ydata_alltp_rs(i_time,:)',Seq(:),nominal(Sub_ID(:)),'VariableNames',{'ERP','Seq','Subject_ID'});
    lme_time = fitlme(tbl_time,'ERP~Seq+(1|Subject_ID)');

    t_map_real(1,i_time) = lme_time.Coefficients{1,4};
    t_map_real(2,i_time) = lme_time.Coefficients{2,4};

    beta_map_real(1,i_time) = lme_time.Coefficients{1,2};
    beta_map_real(2,i_time) = lme_time.Coefficients{2,2};

    CILow_map_real(1,i_time) = lme_time.Coefficients{1,7};
    CILow_map_real(2,i_time) = lme_time.Coefficients{2,7};

    CIHigh_map_real(1,i_time) = lme_time.Coefficients{1,8};
    CIHigh_map_real(2,i_time) = lme_time.Coefficients{2,8};

    p_map_real(1,i_time) = lme_time.Coefficients{1,6};
    p_map_real(2,i_time) = lme_time.Coefficients{2,6};

    clearvars tbl_time lme_time
end

clearvars Ydata_alltp_rs





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Stability in Figure 3 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SIMI_S11_S22 SIMI_S22_S33 SIMI_S33_S44
for s = 1:29

    ind_S11=ind_S1S2(s,:)==1 & ind_subcong(s,:)~=0;
    ind_S22=ind_S2S3(s,:)==1 & ind_subcong(s,:)~=0;
    ind_S33=ind_S3S4(s,:)==1 & ind_subcong(s,:)~=0;

    %% Stability data
    SIMI_S11_S22 (s,:) = mean(SIMI_result_general(s).S11_S22(:,ind_S11),2);
    SIMI_S22_S33 (s,:) = mean(SIMI_result_general(s).S22_S33(:,ind_S22),2);
    SIMI_S33_S44 (s,:) = mean(SIMI_result_general(s).S33_S44(:,ind_S33),2);
end










