function [ M ] =fit_Homography( X, C)
%FIT HOMOGRAPHY fit homography-models according to the clustering C
%   Detailed explanation goes here

if nargin==1 || isempty(C)
    C=ones(size(X,2),1);
end

label=unique(C);
N=length(label); %numero di cluster;
f=size(X,1);
cardmss=4;
M=zeros(9,N);

for i=1:N
    L  = label(i);
    points2fit = X(:,C==L);
    card = sum(C==L);
    x = points2fit(1:2,:);
    y = points2fit(4:5,:);
   	Omega = HomographyDLT(x,y);
   	
%	Omega=vgg_H_from_x_nonlin(Omega,x,y);
%	Omega=vgg_H_from_x_lin(x,y);
%   if(card>8)
%     Omega=vgg_H_from_x_nonlin(Omega,x,y);
%    end
% 	Omega=h2m(x,y);
    M(:,i)=Omega(:);
end
    
    
