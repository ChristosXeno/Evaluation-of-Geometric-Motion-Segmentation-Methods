function flag = validateTheta_homography( Theta)

% flag = validateMSS_homography(X, Theta, s)
%
% DESC:
% Validates the parameter vector
%
% INPUT:
% X                 = samples on the manifold
% Theta             = parameter vector
% s                 = indices of the MSS
%
% OUTPUT:
% flag              = true if the MSS is valid

% condition number threshold
T_K = 1e12;

% perform here the check on the parameter vector Theta
H = reshape(Theta, 3, 3);
K = cond(H);
R=rcond(H);
if (K > T_K || R< 1e-20)
    flag = false;
else
    flag = true;
end;

return
