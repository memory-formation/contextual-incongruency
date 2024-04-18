function m = corr_matrix(X,Y)
% m = corr_matrix(X,Y)
% 
% Fast computation of correlation matrix between two channels*timepoints data sets
% with function fast_corr that enables the quick vectorized computation of 
% pair-wise correlations between corresponding columns of two large matrices.
% 
% 
% 
% Inputs:
% X & Y,    size [channels * timepoints] matrices/vectors; both must
%           be of equal size to enable computation of pair-wise
%           correlations, column-by-column
% 
% Output:
% m,        a correlation matrix of correlation coefficients between X & Y
%           size [timepoints_Y x timepoints_X]
%           -------------------------------------------------------------
%           e.g.       ######################################
%                      ###       _________________        ###
%                      ###    t6|__|__|__|__|__|__|       ###
%                      ###    t5|__|__|__|__|__|__|       ###
%                      ### Y  t4|__|__|__|__|__|__|       ###
%                      ###    t3|__|__|__|__|__|__|       ###
%                      ###    t2|__|__|__|__|__|__|       ###
%                      ###    t1|__|__|__|__|__|__|       ###
%                      ###     0 t1 t2 t3 t4 t5 t6        ###
%                      ###               X                ###
%                      ######################################
% Author:  Xiongbo Wu, May, 2019

    [r1,c1] = size(X); [r2,c2] = size(Y); 
    if r1~=r2 
        error('''X'' and ''Y'' must contain the same number of rows/observations.')
    end
    if c1~=c2 
        error('''X'' and ''Y'' must contain the same number of columns/variables.')
    end
    
    m=zeros(size(X,2),size(Y,2));
    
    
    % Step 1: the first correlation (with exact matching timepoints between two matrices)
     cor=fast_corr(X,Y);
     ind_b=size(Y,2):-1:1;
     ind_a=1:size(X,2);
     idx = sub2ind(size(m), ind_b, ind_a);
     m(idx)=cor;  
     clearvars ind_b  ind_a inx cor
     
     
     % Step 2: For the rest of correlation of matrix with sliding
     % timepoints 
    for t=2:size(Y,2)
        
        % we need to flip the matrix for fast_corr function to work on
        % columnized way
        a=X;
        
        b=Y(:,[t:size(Y,2) 1:t-1]);
        
        
        cor=fast_corr(a,b); % correlation coefficient
        
        % ind_a: column index vector / ind_b: row index vector
        ind_a=1:size(X,2);
        ind_b= [size(Y,2)-t+1:-1:1, size(Y,2):-1:size(Y,2)-t+2 ];
        
        
        
        % idx is the single index number in the correlation matrix 
        % (e.g. 4 stands for the 4th element in the matrix)
        idx = sub2ind(size(m), ind_b, ind_a);
        
        m(idx)=cor;
        clearvars ind_b ind_a idx cor
    end
    
    
end
    
    
