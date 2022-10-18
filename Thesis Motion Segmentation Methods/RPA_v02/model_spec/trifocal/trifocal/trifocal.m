
function [ T ] = trifocal( X )
% Compute the trifocal tensor from N correspondences.
%   INPUT:  
%           X a 9xN matrix of homogeneous points
%   OUTPUT: 
%           T a 27x1 vector containing the trifocal tensor associated to 
%             this triplet of cameras (estimated with linear algorithm)

% non-homogeneous coordinates
Corresp=[X(1:2,:);X(4:5,:);X(7:8,:)];

% Normalization of the data
[x1,Normal1]=Normalize2Ddata(Corresp(1:2,:));
[x2,Normal2]=Normalize2Ddata(Corresp(3:4,:));
[x3,Normal3]=Normalize2Ddata(Corresp(5:6,:));

% Model to estimate T: linear equations
T=linearTFT(x1,x2,x3);

% tensor denormalization
T=transform_TFT(T,Normal1,Normal2,Normal3,1);

T=T(:);

end