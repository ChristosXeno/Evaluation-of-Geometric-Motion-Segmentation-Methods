function [ Z ] = non_degenerate_mss( X,S )
%NON_DEGENERATE_MSS Discard mss that contains collinear three points
%   
[r]=size(S,1);
flg=true(r,1);
for i=1:r
    T=X(:,S(i,:));
    if(isdegenrate(T))
        flg(i)=false;
        disp('degenerate mss')
    end
end
Z=S(flg,:);
end

