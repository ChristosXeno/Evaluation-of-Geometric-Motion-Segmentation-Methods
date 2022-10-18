function [ r ] = islinedegenerate(X, theta, cardmss)
%ISDEGENRATE check if a minimal sample set contains three collinear points
%   r=1 if X is degenrate r=0 otherwise
if nargin==1 
    theta=[];
    cardmss=4;
end



[d,n]=size(X); %number of points
if(n~=cardmss)
    error('minimal sample set of 4 points are needed')
elseif(d>2)
    X=X([1,2],:);
end
R=nan(1,n);
I=[1,2,3,4];
for i=1:4
    T=X(:,I~=i);
    R(i)=iscolinear(T(:,1),T(:,2),T(:,3));

end
r=any(R);
end

