function [ d ] = distPointH( x,H )
%DISTPOINHOMOGRAPHY compute the sampson distance between a couple of point
%arranged in the vector x and an homography H
%   Detailed explanation goes here
H = reshape(H,3,3);
m1 = x(1:3,:);
m2 = x(4:6,:);

d = sqrt(abs(vgg_H_sampson_distance_sqr(H,m1,m2)));
if(isnan(d))
    keyboard
end
end

