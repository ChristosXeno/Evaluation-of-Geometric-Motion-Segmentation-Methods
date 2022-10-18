function [d ] = distPointSeg( P,T )
%DISTANCEPOINTSEG compute the distance between a point P and a segment T
%   Detailed explanation goes here
A=T(1:2);B=T(3:4);
v=B-A;
p=(P-A)'*v/norm(v,2); %proiezione di PA Su AB
if(p<0)% è vicino ad A
    d=norm(P-A);
elseif(norm(p)>norm(v))
    d=norm(P-B);
else
    d=distPointLine(P,cross([A;1],[B;1]));
end


end

