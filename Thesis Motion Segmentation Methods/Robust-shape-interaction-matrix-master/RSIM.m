function [missrate, grp, bestRank, minNcutValue,W] = RSIM(X, s, UpperD, LowerD)
% Input: X --  data matrix, s -- groundtruth label, UpperD -- largest rank
% Pan Ji, pan.ji@anu.edu.au
if(nargin<4)
	LowerD = 1;
end
if(nargin<3)
	UpperD = 4;
end
K = max(s);
r = LowerD*K:UpperD*K; % rank from lower bound K to upper bound 4K
[~,~,VR] = svd(X,'econ'); % take the right singular vector of X
clusterLabel = {};
approxBound = [];
Aff = {};
eigenValues = [];

for ii = 1:length(r)
	rnk = r(ii);
	V = VR(:,1:rnk);
	
	V = normr(V); % normalize each row
	
	Z = V*V'; % new shape interaction matrix	
		
	W = real(Z.^3.5); % enhance block-diagonal structure; 
	                  %	On hopkins155, average err = 0.79% with powering value gamma = 3.8.
					  % You can also try other powering values in [3,5].
					 	
	Aff{ii} = W;
	[clusterLabel{ii},~,~] = ncutW(W,K); % install your ncutW function 
	                                     % from https://www.cis.upenn.edu/~jshi/software/ 
		
	D = diag(1./sum(W));
	L = D*W;	
	eigenValues = eigs(L,K+1);	% you can also easily modify the ncutW function and 
	                            % let it output the eignvalues to save the above three steps
	approxBound(ii) = ComputeNcutValue(W,clusterLabel{ii})/(eigenValues(K)-eigenValues(K+1));
end

[minNcutValue, idx] = min(approxBound);
W = Aff{idx};
grp = clusterLabel{idx};
bestRank = r(idx);

tmp = zeros(length(s),1);
for i = 1:K
    tmp = tmp+grp(:,i)*i;
end
grp = bestMap(s,tmp);
missrate = sum(s(:) ~= grp(:)) / length(s);
%missrate = ErrorRate(grp, s); % calculate the error rate

% figure;
% graphSeg(grp,s);

figure;
subplot(1,2,1); gscatter(X(1,:),X(2,:),s); axis equal; title('GroundTruth'); legend off; set(gca,'Color','k');
subplot(1,2,2); gscatter(X(1,:),X(2,:),grp); axis equal; title('RSIM'); legend off; set(gca,'Color','k');

end
