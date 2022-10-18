function [ M ] = fit_Planes( X,C)
%FIT_LINES  fit line-models according to the clustering C
%   Detailed explanation goes here

if nargin==1 || isempty(C)
    C=ones(size(X,2),1);
end

label=sort(unique(C));
N=length(label); %numero di cluster;
M=nan(4,N);



for i=1:N
    
    L  = label(i);
    points2fit = X(:,C==L);
    U = fitplanes(points2fit);
    M(:,i)=U'./norm(U);
    
end




