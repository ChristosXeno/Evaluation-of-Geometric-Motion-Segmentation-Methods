function [ d ] = distPointPlane( x,H )
%DISTPOINTLINE compute the euclidean distance between a point x and an
%hyperplane H
%   Detailed explanation goes here

d=abs(H'*[x;1])/norm(H);

end

