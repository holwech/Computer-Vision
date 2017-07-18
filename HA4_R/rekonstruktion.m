function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

n = size(Korrespondenzen,2);
X1 = [Korrespondenzen(1:2,:);ones(1,n)];
X2 = [Korrespondenzen(3:4,:);ones(1,n)];
X1 = K\X1;
X2 = K\X2;
    
% compute Matrices M
M = zeros(3*n,n+1,4);    
for i = 1:n
    x2=X2(:,i);
    x1=X1(:,i);
    
    M(i:i+2,i,1) = hat(x2)*R1*x1;
    M(i:i+2,n+1,1) = hat(x2)*T1;
    M(i:i+2,i,2) = hat(x2)*R1*x1;
    M(i:i+2,n+1,2) = hat(x2)*T2;
    M(i:i+2,i,3) = hat(x2)*R2*x1;
    M(i:i+2,n+1,3) = hat(x2)*T1;
    M(i:i+2,i,4) = hat(x2)*R2*x1;
    M(i:i+2,n+1,4) = hat(x2)*T2;
end

% Solve the optimization problem to obtain lambda-vector for each matrix M(:,:,i)
lambdas = zeros(n+1,4);
for i=1:4
    [U,S,V] = svd(M(:,:,i));
    lambdas(:,i) = V(:,size(V,2));
end

% return the euklidean movement for which most of the lambdas are positive
bool = lambdas>0; % returns bool matrix
count_positive = sum(bool); % sums the values in each column and returns row vector
[~,index] = max(count_positive);
switch index
    case 1
        R=R1;
        T=T1;
        lambdas = lambdas(:,1);
    case 2
        R=R1;
        T=T2;
        lambdas = lambdas(:,2);
    case 3
        R=R2;
        T=T1;
        lambdas = lambdas(:,3);
    case 4
        R=R2;
        T=T2;
        lambdas = lambdas(:,4);
    otherwise
        disp('Where am I?')
end

% reconstrucion of P
P1 = X1*diag(lambdas(1:n));
figure;
scatter3(P1(1,:),P1(2,:),P1(3,:),'bo');
hold on
O = -R'*T*lambdas(n+1); % coordinates of camera 2 in coordinate system 1
scatter3(O(1),O(2),O(3),'rd');
scatter3(0,0,0,'rd');
hold off
end

function x_hat = hat(x)
    x_hat = [0, -x(3), x(2);
            x(3), 0, -x(1);
            -x(2), x(1), 0];
end
