
function d = distPointTrifocal(X,T)

T=reshape(T,[3 3 3]);

P=cameras_from_T(T);

d=ReprError(P,X);

%  The error is computed by projecting the space points onto the M images
%  using their projective matrices and computing the RMS of the distances
%  from the repojected point to the original image point.
%
%  Input arguments:
%  ProjM      - 1xM-cell of 3x4 projection matrices.
%  Corresp    - 2MxN matrix with the image points in each image or
%               3MxN matrix if homogeneous coordinates are used.

end

function P=cameras_from_T(T)

% Algorithm 15.1 pag 375 Hartley&Zisserman
P=cell(1,3);

% epipoles
[~,~,V]=svd(T(:,:,1)); v1=V(:,end);
[~,~,V]=svd(T(:,:,2)); v2=V(:,end);
[~,~,V]=svd(T(:,:,3)); v3=V(:,end);
[~,~,V]=svd([v1 v2 v3].'); epi31=V(:,end);

[~,~,V]=svd(T(:,:,1).'); v1=V(:,end);
[~,~,V]=svd(T(:,:,2).'); v2=V(:,end);
[~,~,V]=svd(T(:,:,3).'); v3=V(:,end);
[~,~,V]=svd([v1 v2 v3].'); epi21=V(:,end);

epi31=epi31/norm(epi31);
epi21=epi21/norm(epi21);

% Camera projection matrices
P{1}=eye(3,4);
P{2}=[T(:,:,1)*epi31, T(:,:,2)*epi31, T(:,:,3)*epi31, epi21];
P{3}=[(epi31*epi31'-eye(3))*[T(:,:,1)'*epi21 T(:,:,2)'*epi21 T(:,:,3)'*epi21], epi31];

end
