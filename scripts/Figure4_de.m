
save_root = '';
load([save_root '\ERPs_smoth_down_enc_Exp2.mat'])
load ([save_root '\reject_index_Exp2.mat'])


%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%

time_window = 1:200;
m = 1;

for s = 2:28

    display (['Stability for Subject Number: ' int2str(s) '/' int2str(28)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S11 with S22   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    II = 1:60;
    for i=1:length(II)

        [R,P] = corr (DATA_study2_SIMIdown(s).S11(:,time_window,(II(i))),DATA_study2_SIMIdown(s).S22(:,time_window,(II(i))),'type','pearson');
        SIMI_result_general_Exp2(m).S11_S22(:,i) = fisherz(diag(R));
    end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S22 with S33   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     [I II] = find (reject_index_S22(s,:).*reject_index_S33(s,:) > 0);
    for i=1:length(II)

        [R,P] = corr (DATA_study2_SIMIdown(s).S22(:,time_window,(II(i))),DATA_study2_SIMIdown(s).S33(:,time_window,(II(i))),'type','pearson');
        SIMI_result_general_Exp2(m).S22_S33(:,i) = fisherz(diag(R));
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S33 with S44   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     [I II] = find (reject_index_S33(s,:).*reject_index_S44(s,:) > 0);
    for i=1:length(II)

        [R,P] = corr (DATA_study2_SIMIdown(s).S33(:,time_window,(II(i))),DATA_study2_SIMIdown(s).S44(:,time_window,(II(i))),'type','pearson');
        SIMI_result_general_Exp2(m).S33_S44(:,i) = fisherz(diag(R));
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S44 with S55   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     [I II] = find (reject_index_S44(s,:).*reject_index_S55(s,:) > 0);
    for i=1:length(II)

        [R,P] = corr (DATA_study2_SIMIdown(s).S44(:,time_window,(II(i))),DATA_study2_SIMIdown(s).S55(:,time_window,(II(i))),'type','pearson');
        SIMI_result_general_Exp2(m).S44_S55(:,i) = fisherz(diag(R));
    end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S55 with S66   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     [I II] = find (reject_index_S55(s,:).*reject_index_S66(s,:) > 0);
    for i=1:length(II)

        [R,P] = corr (DATA_study2_SIMIdown(s).S55(:,time_window,(II(i))),DATA_study2_SIMIdown(s).S66(:,time_window,(II(i))),'type','pearson');
        SIMI_result_general_Exp2(m).S55_S66(:,i) = fisherz(diag(R));
    end

    m = m+1;
end




%%%%%%%%%%%%%%%%%%%%%%%%
%% Averaged Stability %%
%%%%%%%%%%%%%%%%%%%%%%%%

reject_index_S11 = reject_index_S11;
reject_index_S22 = reject_index_S22;
reject_index_S33 = reject_index_S33;
reject_index_S44 = reject_index_S44;
reject_index_S55 = reject_index_S55;
reject_index_S66 = reject_index_S66;

%% Prepare structure
ind_S1S2 = reject_index_S11.*reject_index_S22;
ind_S2S3 = reject_index_S22.*reject_index_S33;
ind_S3S4 = reject_index_S33.*reject_index_S44;
ind_S4S5 = reject_index_S44.*reject_index_S55;
ind_S5S6 = reject_index_S55.*reject_index_S66;

%% Exclude possible NaNs EEG trials

NaN1 = ones (27,60);
NaN2 = ones (27,60);
NaN3 = ones (27,60);
NaN4 = ones (27,60);
NaN5 = ones (27,60);
NaN6 = ones (27,60);

for s = 1:27

    ind_S11=ind_S1S2(s,:)==1;
    ind_S22=ind_S2S3(s,:)==1;
    ind_S33=ind_S3S4(s,:)==1;
    ind_S44=ind_S4S5(s,:)==1;
    ind_S55=ind_S5S6(s,:)==1;

    [row_1, col_1] = find(isnan(sum(SIMI_result_general_Exp2(s).S11_S22(:,:),1)));
    [row_2, col_2] = find(isnan(sum(SIMI_result_general_Exp2(s).S22_S33(:,:),1)));
    [row_3, col_3] = find(isnan(sum(SIMI_result_general_Exp2(s).S33_S44(:,:),1)));
    [row_4, col_4] = find(isnan(sum(SIMI_result_general_Exp2(s).S44_S55(:,:),1)));
    [row_5, col_5] = find(isnan(sum(SIMI_result_general_Exp2(s).S55_S66(:,:),1)));

    if length(col_1) > 0; NaN1 (s,col_1(1,:)) = 0; end
    if length(col_2) > 0; NaN2 (s,col_2(1,:)) = 0; end
    if length(col_3) > 0; NaN3 (s,col_3(1,:)) = 0; end
    if length(col_4) > 0; NaN4 (s,col_4(1,:)) = 0; end
    if length(col_5) > 0; NaN5 (s,col_5(1,:)) = 0; end
end

ind_S1S2 = reject_index_S11.*reject_index_S22.*NaN1;
ind_S2S3 = reject_index_S22.*reject_index_S33.*NaN2;
ind_S3S4 = reject_index_S33.*reject_index_S44.*NaN3;
ind_S4S5 = reject_index_S44.*reject_index_S55.*NaN4;
ind_S5S6 = reject_index_S55.*reject_index_S66.*NaN5;


for s = 1:27

    a = ind_S1S2(s,:);
    [U UU] = find(a > 0);
    SIMI_exp2(s,1,:) = mean(SIMI_result_general_Exp2(s).S11_S22(:,UU),2)';

    a = ind_S2S3(s,:);
    [U UU] = find(a > 0);
    SIMI_exp2(s,2,:) = mean(SIMI_result_general_Exp2(s).S22_S33(:,UU),2)';

    a = ind_S3S4(s,:);
    [U UU] = find(a > 0);
    SIMI_exp2(s,3,:) = mean(SIMI_result_general_Exp2(s).S33_S44(:,UU),2)';

    a = ind_S4S5(s,:);
    [U UU] = find(a > 0);
    SIMI_exp2(s,4,:) = mean(SIMI_result_general_Exp2(s).S44_S55(:,UU),2)';

    a = ind_S5S6(s,:);
    [U UU] = find(a > 0);
    SIMI_exp2(s,5,:) = mean(SIMI_result_general_Exp2(s).S55_S66(:,UU),2)';
end

clearvars -except SIMI*






%%%%%%%%%%%%%%%%%%%%%%%%
%% LMM %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%% Linear Mixed-effect Model on Stability data

save_root = '';

%% Prepare structure
ind_S1S2 = reject_index_S11.*reject_index_S22;
ind_S2S3 = reject_index_S22.*reject_index_S33;
ind_S3S4 = reject_index_S33.*reject_index_S44;
ind_S4S5 = reject_index_S44.*reject_index_S55;
ind_S5S6 = reject_index_S55.*reject_index_S66;

%% Exclude possible NaNs EEG trials

NaN1 = ones (27,60);
NaN2 = ones (27,60);
NaN3 = ones (27,60);
NaN4 = ones (27,60);
NaN5 = ones (27,60);
NaN6 = ones (27,60);

m = 1;
for s = 1:27

    ind_S11=ind_S1S2(m,:)==1;
    ind_S22=ind_S2S3(m,:)==1;
    ind_S33=ind_S3S4(m,:)==1;
    ind_S44=ind_S4S5(m,:)==1;
    ind_S55=ind_S5S6(m,:)==1;

    [row_1, col_1] = find(isnan(SIMI_result_general_Exp2(s).S11_S22(:,:)));
    [row_2, col_2] = find(isnan(SIMI_result_general_Exp2(s).S22_S33(:,:)));
    [row_3, col_3] = find(isnan(SIMI_result_general_Exp2(s).S33_S44(:,:)));
    [row_4, col_4] = find(isnan(SIMI_result_general_Exp2(s).S44_S55(:,:)));
    [row_5, col_5] = find(isnan(SIMI_result_general_Exp2(s).S55_S66(:,:)));

    if length(col_1) > 0; NaN1 (m,col_1(1,:)) = 0; end
    if length(col_2) > 0; NaN2 (m,col_2(1,:)) = 0; end
    if length(col_3) > 0; NaN3 (m,col_3(1,:)) = 0; end
    if length(col_4) > 0; NaN4 (m,col_4(1,:)) = 0; end
    if length(col_5) > 0; NaN5 (m,col_5(1,:)) = 0; end
    m = m+1;
end

ind_S1S2 = reject_index_S11.*reject_index_S22.*NaN1;
ind_S2S3 = reject_index_S22.*reject_index_S33.*NaN2;
ind_S3S4 = reject_index_S33.*reject_index_S44.*NaN3;
ind_S4S5 = reject_index_S44.*reject_index_S55.*NaN4;
ind_S5S6 = reject_index_S55.*reject_index_S66.*NaN5;


Sub_ID = zeros(0); Seq = zeros(0); Y_all = zeros(0); 


for s = 1:27

    display (['Preparing LMM for Subject Number: ' int2str(s) '/' int2str(29)]);

    ind_S11=ind_S1S2(s,:)==1;
    ind_S22=ind_S2S3(s,:)==1;
    ind_S33=ind_S3S4(s,:)==1;
    ind_S44=ind_S4S5(s,:)==1;
    ind_S55=ind_S5S6(s,:)==1;


    n_S11 = sum(ind_S11);n_S22 = sum(ind_S22);n_S33 = sum(ind_S33);
    n_S44 = sum(ind_S44); n_S55 = sum(ind_S55);

    %% sub ID
    Sub_ID = cat(2,Sub_ID,rep_num(s,1,n_S11+n_S22+n_S33+n_S44+n_S55));


    %% Seq
    Seq_sub = [rep_num(1,1,n_S11),rep_num(2,1,n_S22),rep_num(3,1,n_S33),rep_num(4,1,n_S44),rep_num(5,1,n_S55)];
    Seq = cat(2,Seq,Seq_sub); clearvars Seq_sub


    %% Stability data
    SIMI_S11_S22 = SIMI_result_general_Exp2(s).S11_S22(:,ind_S11);
    SIMI_S22_S33 = SIMI_result_general_Exp2(s).S22_S33(:,ind_S22);
    SIMI_S33_S44 = SIMI_result_general_Exp2(s).S33_S44(:,ind_S33);
    SIMI_S44_S55 = SIMI_result_general_Exp2(s).S44_S55(:,ind_S44);
    SIMI_S55_S66 = SIMI_result_general_Exp2(s).S55_S66(:,ind_S55);


    SIMI_sub = cat(2,SIMI_S11_S22,SIMI_S22_S33,SIMI_S33_S44,SIMI_S44_S55,SIMI_S55_S66);
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

