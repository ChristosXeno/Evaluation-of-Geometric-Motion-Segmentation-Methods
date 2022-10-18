function [ ] = drawLines( M ,c)
%DRAWLINES Summary of this function goes here
%   Detailed explanation goes here
  
if(nargin==1 || isempty(c))
    c=[0, 0 ,0];
end
    p=[-M(1)/M(2),-M(3)/M(2)];
        xx=linspace(-1,1);
        yy=polyval(p,xx);
        hold on
        xlim([-1 1]);
        ylim([-1 1]);
        plot(xx,yy,'-','Markersize',5, 'Color',c);
        
        axis equal
end

