%Generate data for the tests from the ground truth
% - ngroup      number of groups
% - PPG         points in each group
% - xord        normalized coordinates ordered per group
% - yord        unnormalized coordinates ordered per group
% - xordindex   vector of indeces to pass from unordered coordinates to
%               coordinates ordered per group
% - x1ord       coordinates normalized with Hartley's normalization
%

%order the true segmentation in groups
ngroups=max(G);                                             %total number of groups
xordindex=[];
PPG=[];
x1ord = [];
for (i=1:ngroups)
    xordindex=[xordindex; find(G==i)];                      %index which pass from non-ordered to ordered features
    PPG=[PPG length(find(G==i))];                               %array with the number of features for each group
end
xord=x(:,xordindex,:);
yord=y(:,xordindex,:);
for(i=1:size(x,3))
    x1ord(:,:,i)=normalise2dpts(yord(:,:,i));
end