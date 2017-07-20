function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden
[U,~,V] = svd(E);
 
%Make sure that E is normalized (singular values 1, 1, 0)
S = diag([1 1 0]);

%Make sure U and V are orthonormal rotationmatrices ==> det(U) = 1 det(V)=1
d=sign(det(U));
U = U*[1, 0, 0; 
    0, 1, 0; 
    0, 0, d];
d=sign(det(V));
V = [1, 0, 0;
    0, 1, 0;
    0, 0, d]*V;


R_z1 = [ 0, -1, 0;
    1, 0, 0;
    0, 0, 1];
R_z2 = R_z1';

T1_hat = U*R_z1*S*U';
T2_hat = U*R_z2*S*U';
R1 = U*R_z1'*V';
R2 = U*R_z2'*V';

% rekonstrukt T from T_hat
T1 = [T1_hat(3,2); T1_hat(1,3); T1_hat(2,1)];
T2 = [T2_hat(3,2); T2_hat(1,3); T2_hat(2,1)];

end