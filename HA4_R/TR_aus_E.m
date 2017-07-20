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
 
 T1_hat = U*R_z1*S*U';
 T2_hat = U*R_z2*S*U';
 R1 = U*R_z1'*V';
 R2 = U*R_z2'*V';
 
 % rekonstrukt T from T_hat
 T1 = [T1_hat(3,2); T1_hat(3,1)*-1; T1_hat(2,1)];
 T2 = [T2_hat(3,2); T2_hat(3,1)*-1; T2_hat(2,1)];
% E=eye(3);
%  T1 = 0;
%  for i = 1:3
%      for j = 1:3
%          T1 = T1 + cross(T1_hat(i,j)*E(:,i),E(:,j));
%          T1 = -(1.0/2) * T1;
%      end
%  end
%  T2 = 0;
%  for i = 1:3
%      for j = 1:3
%          T2 = T2 + cross(T2_hat(i,j)*E(:,i),E(:,j));
%          T2 = -(1.0/2) * T2;
%      end
%  end
end