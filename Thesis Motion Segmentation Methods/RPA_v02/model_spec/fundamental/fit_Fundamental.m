function [ M ] = fit_Fundamental( X, C )
%RANSAC_FIT Summary of this function goes here
%   Detailed explanation goes here

if  ~exist('C','var') || isempty(C)

    x = X(1:3,:);
    y = X(4:6,:);
    
    
    f_start = fund(x,y);
    f_start=reshape(f_start,3,3);
    f_nonlin = F_from_x_nonlin(f_start,x,y);
    M = f_nonlin(:);
    
else
    
   label = sort(unique(C));
    N = length(label); % number of cluster
    M = nan(9,N);      % preallocate N models
    
    for i=1:N
        
        L  = label(i);
        %if(sum(Cs==L)>8)
        points2fit = X(:,C==L);
        x = points2fit(1:3,:);
        y = points2fit(4:6,:);
         
        f_start = fund(x,y);
        f_start = reshape(f_start ,3,[]); 
        % non linear refinement
        
        f_nonlin = F_from_x_nonlin(f_start,x,y);
        
%         m = [x(1:2,:); y(1:2,:)];
%         [m,T1,T2]=normalonetoone(m);
%         f= funmatMestTorr(m,10,10^(-5));
%         f=T1'*f*T2;
%         f=f./norm(f);
        
        M(:,i) = f_nonlin(:);
        %end
        
    end
end


