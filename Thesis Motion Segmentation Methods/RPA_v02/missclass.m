
% Segmentation = labels found by motion segmentation algorithm = group
% npoints = points in each group = N
% ngroups = number of groups = 3

function [miss] = missclass(Segmentation,npoints,ngroups,gtSeg)

% [miss,index] = missclass(Segmentation,npoints,ngroups)
%
% Computes the number of missclassified points in the vector Segmentation. 
%
% Segmentation: 1 by sum(npoints) or sum(ngroups) by 1 vector containing 
% the label for each group, ranging from 1 to n

% npoints: 1 by ngroups or ngroups by 1 vector containing the number of 
% points in each group.

% ngroups: number of groups

Permutations = perms(1:ngroups);
if(size(Segmentation,2)==1)
    Segmentation=Segmentation';
end

miss = zeros(size(Permutations,1),size(Segmentation,1));

tempSeg = Segmentation;
inds = {};
err_rates = [];
    
for j=1:size(Permutations,1)
        
    for k = 1:ngroups
            
        inds{k} = find(tempSeg==k);
        Segmentation(inds{k}) = Permutations(j,k);
    end
        
    err_rates(j) = errorCalc(gtSeg,Segmentation);
        
end

miss = min(err_rates);
    
end