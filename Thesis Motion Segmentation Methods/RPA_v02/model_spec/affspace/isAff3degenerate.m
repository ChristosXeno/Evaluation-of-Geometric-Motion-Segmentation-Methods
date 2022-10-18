function  y  = isAff3degenerate( X, theta , cardmss )
%ISHDEGENERATE Check if a sample of 4 points X is degenerate i.e <X> <4
%   INPUT:
%          X: fx4 matrix-sample of four points in homogeneous coordinate
%          in two views
%          y: boolean- true if X is H-degenerate.
%          theta: parameter vectors of subspace matrices in columns

%%
tol=1e-5;
[~,lambda,~]=svd(X(:,2:end));
L=diag(lambda);
y=any(L<tol);
%disp('degenerate')
