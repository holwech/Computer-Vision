function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden

Rz_pos = [0 -1 0; 1 0 0; 0 0 1];
Rz_neg = [0 1 0; -1 0 0; 0 0 1];

[U, S, V] = svd(E);
%Make sure U and V are orthonormal
d1 = int16(det(U));
d2 = int16(det(V));
if (d1 ~= 1 || d2 ~= 1)
    disp('U and V are not rotation matrices!');
end

T1_hat = U*Rz_pos*S*U';
T1 = [T1_hat(3,2); T1_hat(3,1)*-1; T1_hat(2,1)];
R1 = U*Rz_pos'*V';
T2_hat = U*Rz_neg*S*U';
T2 = [T2_hat(3,2); T2_hat(3,1)*-1; T2_hat(2,1)];
R2 = U*Rz_neg'*V';

end