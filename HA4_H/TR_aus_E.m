function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden

Rz_pos = [0 -1 0; 
        1 0 0; 
        0 0 1];
Rz_neg = [0 1 0; 
        -1 0 0; 
        0 0 1];

%Make sure that E is normalized (singular values 1, 1, 0)
[U, ~, V] = svd(E);
S = diag([1 1 0]);
E = U*S*V';

%Find T_hat and R by using svd
[U, ~, V] = svd(E);

%Make sure U and V are orthonormal rotationmatrices ==> det(U) = 1 det(V)=1

d1 = int16(det(U));
d2 = int16(det(V));
if (d1 ~= 1)
    U = U*[1 0 0; 0 1 0; 0 0 -1];
end
if (d2 ~= 1)
    V = V*[1 0 0; 0 1 0; 0 0 -1];
end
% d1 = int16(det(U))
% d2 = int16(det(V))
% Error = U*S*V' - E

%[U, S, V] = svd(E);

T1_hat = U*Rz_pos*S*U';
T1 = [T1_hat(3,2); T1_hat(3,1)*-1; T1_hat(2,1)];
R1 = U*Rz_pos'*V';

T2_hat = U*Rz_neg*S*U';
T2 = [T2_hat(3,2); T2_hat(3,1)*-1; T2_hat(2,1)];
R2 = U*Rz_neg'*V';

end