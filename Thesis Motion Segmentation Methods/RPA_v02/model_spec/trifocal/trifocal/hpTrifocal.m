function [ H ] = hpTrifocal( X, S )

% TO DO: versione con soluzione multiple 

m=size(S,1);
H=zeros(27,m); % preallocating for speed

for i=1:m
    
    y=X(:,S(i,:));
    
    % T= trifocal_nonlin( y );
    T= trifocal( y );

    H(:,i)=T;

end

end

