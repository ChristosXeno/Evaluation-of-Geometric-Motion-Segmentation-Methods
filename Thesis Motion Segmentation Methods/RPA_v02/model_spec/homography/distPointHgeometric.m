function [ d ] = distPointHgeometric( x,H )
%DISTPOINHOMOGRAPHY compute the sampson distance between a couple of point
%arranged in the vector x and an homography H
%   Detailed explanation goes here
H=reshape(H,3,3);
m1=x(1:3,:);
m2=x(4:6,:);

% forward mapping
mf = homo2cart(H*m1);
% backward mapping
mb = homo2cart(H\m2);


d1 = sum((m1(1:2, :)-mb).^2, 1); % distance in the first image
d2 = sum((m2(1:2, :)-mf).^2, 1); % distance in the second image

d = d1 + d2;

end

