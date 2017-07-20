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
T_set = {T1,T2};
R_set = {R1,R2};
for i = 1:n
    set=1;
    for k=1:2
        for l=1:2
            M(3*i:3*i+2,i,set) = hat(X2(:,i))*R_set{k}*X1(:,i);
            M(3*i:3*i+2,n+1,set) = hat(X2(:,i))*T_set{l};
            set=set+1;
        end
    end
end

% Solve the optimization problem to obtain lambda-vector for each matrix M(:,:,i)
lambdas = zeros(n+1,4);
for set=1:4
    [~,~,V] = svd(M(:,:,set));
    lambdas(:,set) = V(:,n+1)/V(n+1,n+1);
end

% return the euklidean movement for which most of the lambdas are positive
bool = lambdas>0; % returns bool matrix
count_positive = sum(bool,1); % sums the values in each column and returns row vector
[~,index] = max(count_positive);
switch index
    case 1
        R=R1;
        T=T1;
    case 2
        R=R1;
        T=T2;
    case 3
        R=R2;
        T=T1;
    case 4
        R=R2;
        T=T2;
    otherwise
        disp('Where am I?')
end
lambdas=lambdas(:,index);
% reconstruction of P
P1 = X1*diag(lambdas(1:n));

% 3D Plot
figure;
scatter3(P1(1,:),P1(2,:),P1(3,:),'bo');
hold on
O = -R'*T; % coordinates of camera 2 in coordinate system 1
plotCamera('Location',[O(1),O(2),O(3)],'Color',[1,0,0], 'Size', 0.2)
plotCamera('Location',[0,0,0],'Color',[1,0,0], 'Size', 0.2)
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on 
axis equal
hold off
end

function x_hat = hat(x)
    x_hat = [0, -x(3), x(2);
            x(3), 0, -x(1);
            -x(2), x(1), 0];
end
