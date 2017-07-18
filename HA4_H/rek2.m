function [T,R, lambdas, P1] = rek2(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1



% compute Matrices M
M = zeros(3*size(Korrespondenzen,2),size(Korrespondenzen,2)+1,4);
for i = 1:size(Korrespondenzen,2)
    x2=K\[Korrespondenzen(3:4,i);1];
    x1=K\[Korrespondenzen(1:2,i);1];
    M(i:i+2,i,1) = hat(x2)*R1*x1;
    M(i:i+2,size(Korrespondenzen,2)+1,1) = hat(x2)*T1;
    M(i:i+2,i,2) = hat(x2)*R1*x1;
    M(i:i+2,size(Korrespondenzen,2)+1,2) = hat(x2)*T2;
    M(i:i+2,i,3) = hat(x2)*R2*x1;
    M(i:i+2,size(Korrespondenzen,2)+1,3) = hat(x2)*T1;
    M(i:i+2,i,4) = hat(x2)*R2*x1;
    M(i:i+2,size(Korrespondenzen,2)+1,4) = hat(x2)*T2;
end

% Solve the optimization problem to obtain lambda-vector for each matrix M(:,:,i)
lambdas = zeros(size(Korrespondenzen,2)+1,4);
for j=1:4
    [U,S,V] = svd(M(:,:,j));
    lambdas(:,j) = V(:,size(Korrespondenzen,2)+1);
end

% return the euklidean movement for which most of the lambdas are positive
bool = lambdas>0; % returns bool matrix
count_positive = sum(bool); % sums the values in each column and returns row vector
[~,index] = max(count_positive);

n = size(Korrespondenzen,2);
gamma = lambdas(n + 1, index);
switch index
    case 1
        R=R1;
        T=gamma*T1;
    case 2
        R=R1;
        T=gamma*T2;
    case 3
        R=R2;
        T=gamma*T1;
    case 4
        R=R2;
        T=gamma*T2;
    otherwise
        disp('Where am I?')
end

lambdas(:, index);
x1 = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
x1 = K\x1;
P1 = zeros(3, n);
for i = 1:n
    P1(:,i) = lambdas(i,1)*x1(:,i);
end

end

function x_hat = hat(x)
    x_hat = [0, -x(3), x(2);
            x(3), 0, -x(1);
            -x(2), x(1), 0];
end
