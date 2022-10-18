function [ H] = hpSeg( X,S )
%HPSEG Summary of this function goes here
%   Detailed explanation goes here
m=size(S,1);
H=zeros(4,m); 
for i=1:m
    a=[X(:,S(i,1))];
    b=[X(:,S(i,2))];
    H(:,i)=[a;b];
    
end

end

