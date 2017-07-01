function [Korrespondenzen_robust] = F_ransac(Korrespondenzen,varargin)
% Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
% robusten Korrespondenzpunktpaaren
epsilon=50;
p=90;
tolerance=100;
if(nargin>1)
    epsilon=varargin{1};
    if(nargin>2)
        p=varargin{2};
        if(nargin>3)
            tolerance=varargin{3};
        end
    end
end
k=8;
s = log(1-p)/log(1-(1-epsilon)^k);

end
end