
function [ M ] =fit_trifocal( X, C)
% fit trifocal models according to the clustering C

if nargin==1 || isempty(C)
    C=ones(size(X,2),1);
end

label=unique(C);
N=length(label); % number of clusters
M=nan(27,N); % preallocate N models

for i=1:N
    
    L  = label(i);
    points2fit = X(:,C==L);
    
    T= trifocal_nonlin( points2fit );
    %T= trifocal( points2fit );
   	
    M(:,i)=T;
    
end


end