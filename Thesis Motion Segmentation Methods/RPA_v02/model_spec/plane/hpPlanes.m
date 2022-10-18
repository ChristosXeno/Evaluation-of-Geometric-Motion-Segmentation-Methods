function [ H ] = hpPlanes( X,S)
%HPPLANES
m=size(S,1);
H=zeros(4,m); 
for i=1:m
   H(:,i)= fitplanes(X(:,S(i,:)));   
    
end
