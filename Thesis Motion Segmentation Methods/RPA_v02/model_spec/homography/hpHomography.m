function [ H ] = hpHomography( X, S )
%HPFUNDAMENTAL Summary of this function goes here

m=size(S,1);
H=zeros(9,m); % preallocating for speed
for i=1:m

    x = X(1:3,S(i,:));
    y = X(4:6,S(i,:));

    
    %Omega=HomographyDLT(x,y);
    Omega= vgg_H_from_x_lin(x,y);
    %[Omega,~] = vgg_H_from_x_nonlin(Omega,x,y);

    H(:,i)=Omega(:);

    
end



