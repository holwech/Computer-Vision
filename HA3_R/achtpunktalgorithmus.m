function [EF] = achtpunktalgorithmus(Korrespondenzen,varargin)
% Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
% mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
% vorliegt oder nicht
cali=0;
if(nargin==2)
    K=varargin{1};
    cali=1;
end

% Kronecker Product of correspondences
A = zeros(size(Korrespondenzen,2),9);
for k=1:size(Korrespondenzen,2)
    x1=[Korrespondenzen(1:2,k);1];
    x2=[Korrespondenzen(3:4,k);1];
    if(cali==1)
        x1=K\x1;
        x2=K\x2;
    end
    A(k,:)=kron(x1,x2)';
end

% minimize Epipolarequality
[U_A,S_A,V_A] = svd(A);
G = [V_A(1:3,9),V_A(4:6,9),V_A(7:9,9)];
[U_G,S_G,V_G] = svd(G);

if(cali==0)
    % estimate F
    S_G(3,3)=0;
    EF = U_G*S_G*V_G';
    
elseif(cali==1)
    % estimate E
    S=[1,0,0;0,1,0;0,0,0];
    EF = U_G*S*V_G';
    % 3D-Reconstruction?
end

end