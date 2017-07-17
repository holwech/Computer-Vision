function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden
 [U,S,V] = svd(E);
 R_z1 = [ 0, -1, 0;
        1, 0, 0;
        0, 0, 1];
 R_z2 = [ 0, 1, 0;
        -1, 0, 0;
        0, 0, 1];   
 
 T1_hat=U*R_z1*S*U';
 T2_hat=U*R_z2*S*U';
 R1=U*R_z1'*V';
 R2=U*R_z2'*V';
 
 % rekonstrukt T from T_hat
 T1=
 T2=
end