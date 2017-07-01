function [Korrespondenzen_robust] = F_ransac(Korrespondenzen,varargin)
% Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
% robusten Korrespondenzpunktpaaren

% define epsilon,p and tolerance
epsilon=0.50;
p=0.9;
tolerance=1000;
if(nargin>1)
    epsilon=varargin{1};
    if(nargin>2)
        p=varargin{2};
        if(nargin>3)
            tolerance=varargin{3};
        end
    end
end
load('K.mat');
% compute how many models we compute
S = floor(log(1-p)/log(1-(1-epsilon)^8))

amount_best=0;
sumdist_best=0;
ConsensusSet_best=zeros(size(Korrespondenzen,2),1);
    for s=1:S
        Korr=choose_random(Korrespondenzen);
        F=achtpunktalgorithmus(Korr,K);
        
        ConsensusSet=zeros(size(Korrespondenzen,2),1);
        amount=0;
        sumdist=0;
        for k=1:size(Korrespondenzen,2)
            dist=sampson_dist(Korrespondenzen(:,k),F);
            if dist<tolerance
                amount = amount+1;      % amount of correspondences fitting the model s
                sumdist = sumdist+dist; % sum of sampson-distances for model s
                % (sum of all pairs or only of pairs that are part of consensus set?)
                ConsensusSet(k)=1;
            end
        end

        % Find the best model
        if amount>amount_best
            F_best=F;
            amount_best=amount;
            sumdist_best=sumdist;
            ConsensusSet_best=ConsensusSet;
        elseif amount==amount_best
            if sumdist>sumdist_best
                F_best=F;
                sumdist_best=sumdist;
                ConsensusSet_best=ConsensusSet;
            end
        end
    end
Korrespondenzen_robust=Korrespondenzen(logical(ConsensusSet_best));
%display(F_best);
end

function eight_correspondences = choose_random(Korrespondenzen)
    k=floor(rand(8,1)*size(Korrespondenzen,2))+1;
    eight_correspondences = Korrespondenzen(k);
end

function dist = sampson_dist(Ko,F)
    x1=[Ko(1:2,:);1];
    x2=[Ko(3:4);1];
    e3=[0,0,1]';
    dist = ( (x2'*F*x1)^2 )/( norm(cross(e3,F*x1)) + norm(cross(x2'*F,e3')) );
end
