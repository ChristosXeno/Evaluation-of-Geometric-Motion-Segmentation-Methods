function [ H ] = hpAffspace( X,S )
%HPAFFSSPACE generate affine subspace. Each subspace is rapresented as
%L(p,U) the set of point q= p+v s.t. v\in U: for this reason a space is
%represented as a matrix where the first column collect the components of p
%and the others the component of an orthonormal basis for U
%

[f,n] = size(X);% lunghezza della traiettoria: delle colonne di X
m = size(S,1); %numero di ipotesi
cardmss=size(S,2);
d= cardmss-1; %dimensioni del sottospazio
H = zeros(f*cardmss,m);



for i=1:m
    points2fit=X(:,S(i,:));
  
    centroid=sum(points2fit,2)./n;
    momentum_matrix=zeros(f);
    for j=1:cardmss
        momentum_matrix=momentum_matrix+(points2fit(:,j)-centroid)*(points2fit(:,j)-centroid)';
    end
    [U,~,~]=svd(momentum_matrix);
    U=U(:,1:d);
    U=U(:);
    
    H(:,i)=vertcat(centroid, U);

    
end

%alternativa 2


