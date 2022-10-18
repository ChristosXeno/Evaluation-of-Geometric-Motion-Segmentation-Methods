
% Segmentation = labels found by motion segmentation algorithm = group
% npoints = points in each group = N
% ngroups = number of groups = 3

function [miss,index] = missclass(Segmentation,npoints,ngroups)

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
for k=1:size(Segmentation,1)
    for j=1:size(Permutations,1)
        
        missingIndexes1 = find(Segmentation(k,1:npoints(1))~=Permutations(j,1));
        
        
        miss(j,k) = length(missingIndexes1);
        
        
        for i=2:ngroups
            
            missingIndexes2 = find(Segmentation(k,sum(npoints(1:i-1))+1:sum(npoints(1:i)))~=Permutations(j,i));
            
            miss(j,k) = miss(j,k) + length(missingIndexes2);
        end
    end
end

[miss,temp] = min(miss,[],1);
miss=mean(miss);
index = Permutations(temp,:);