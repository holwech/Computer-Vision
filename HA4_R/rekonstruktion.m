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
M_inv = zeros(3*n,n+1,4);
T_cell = {T1,T2,T1,T2};
R_cell= {R1,R1,R2,R2};
for i = 1:n
    for set = 1:4
        M(3*(i-1)+1:3*i,i,set) = hat(X2(:,i))*R_cell{set}*X1(:,i);
        M(3*(i-1)+1:3*i,n+1,set) = hat(X2(:,i))*T_cell{set};
        M_inv(3*(i-1)+1:3*i,i,set) = hat(X1(:,i))*R_cell{set}'*X2(:,i);
        M_inv(3*(i-1)+1:3*i,n+1,set) = hat(X1(:,i))*T_cell{set};
    end
end

% Solve the optimization problem to obtain lambda-vector for each matrix M(:,:,i)
lambdas = zeros(n+1,4);
lambdas_inv = zeros(n+1,4);
for set=1:4
    [~,~,V] = svd(M(:,:,set));
    lambdas(:,set) = V(:,n+1)/V(n+1,n+1);
    [~,~,V] = svd(M_inv(:,:,set));
    lambdas_inv(:,set) = V(:,n+1)/V(n+1,n+1);
end

% return the euklidean movement for which most of the lambdas are positive
bool = uint8(lambdas>0) % returns matrix with zeros and ones
count_positive = sum(bool,1) % sums the values in each column and returns row vector
% bool = lambdas_inv>0; % returns bool matrix
% count_positive_inv = sum(bool,1); % sums the values in each column and returns row vector
% if max(count_positive_inv)>max(count_positive)
%     count_positive = count_positive_inv;
% end
[~,index] = max(count_positive)
R = R_cell{index};
T = T_cell{index};
% reconstruction of P
P1 = X1*diag(lambdas(1:n,index));
%lambdas=lambdas(:,index);

% 3D Plot
figure;
%subplot(1,3,2)
plot(P1,R,T)
title('chosen R,T')

for index = 1:4
    
    R = R_cell{index};
    T = T_cell{index};


    % reconstruction of P
    P1 = X1*diag(lambdas(1:n,index));

    % 3D Plot
    figure
    %subplot(index+2,3,2)
    plot(P1,R,T)
end

end

function x_hat = hat(x)
    x_hat = [0, -x(3), x(2);
            x(3), 0, -x(1);
            -x(2), x(1), 0];
end

function plot(P1,R,T)
    scatter3(P1(1,:),P1(2,:),P1(3,:),'bo');
    hold on
    O = -R'*T; % coordinates of camera 2 in coordinate system 1
    plotCamera('Location',[O(1),O(2),O(3)],'Color',[1,0,0],'Orientation',R, 'Size', 0.2)
    plotCamera('Location',[0,0,0],'Color',[1,0,0], 'Size', 0.2)
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on
    axis equal
    rotate3d
    hold off
end