function [ d ] = distPointAffspace( X, L )
%DISTPOINTAFFSPACE compute the distance between a point X and an affine
%subspace L


f=size(X,1);
L=reshape(L,f,4);
% proietto  X su U ottenendo Y
U=L(:,2:end); %giacitura di L
p=L(:,1); % punto di L
k=(X-p)'*U; %coefficienti di Y in U
Y=p;
for i=1:length(k)
Y = Y + k(i).*U(:,i); %proiezione di X di L
end
d= norm(X-Y);


