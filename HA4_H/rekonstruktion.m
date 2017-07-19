function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

%Want to solve the equation [x2_hat*R*x1 x2_hat*T][lambda; gamma] = 0
% => Md = 0, by minimization of Md with respect to d and norm(d) = 1;

n = size(Korrespondenzen, 2);

%% Preprocess the KPs
%Put the pixel coordinates in homogenous form and use the calibration
%matrix K to find the homogenous coordinates.
x1 = [Korrespondenzen(1:2,:);ones(1,n)];
x2 = [Korrespondenzen(3:4,:);ones(1,n)];
x1 = K\x1;
x2 = K\x2;

%% Constuct M for all four possible (R,T):

%We have four combinations of (R1,R2,T1,T2)
%(R1,T1) ==> E
%(R2,T2)
%(R1,T2) ==> -E
%(R2,T1)

diagonal_elements = {zeros(3,n), zeros(3,n), zeros(3,n), zeros(3,n)};
column_elements = {zeros(3,n), zeros(3,n), zeros(3,n), zeros(3,n)};
R_cell = {R1, R1, R2, R2};
T_cell = {T1, T2, T1, T2};
M_diag_cell = {0, 0, 0, 0};
M_cell = {0, 0, 0, 0};

%Calculate x2_hat*R*x1 and x2_hat*T for all KP
for set = 1:4
    for i = 1:n
        x2_hat = skew(x2(:,i));
        diagonal_elements{set}(:,i) = x2_hat*R_cell{set}*x1(:,i);
        column_elements{set}(:,i) = x2_hat*T_cell{set};
    end

    %Gather in matrix form
    [r,c] = size(diagonal_elements{set});
    i = 1:numel(diagonal_elements{set});
    j = repmat(1:c,r,1);
    M_diag_cell{set} = full(sparse(i',j(:),diagonal_elements{set}(:)));
    M_cell{set} = [M_diag_cell{set} reshape(column_elements{set},3*n, 1)];
    %M_cell{set} = full(M_cell{set});
    
end
%M_cell{1}

%% Minimize M*d using svd
d = {zeros(n + 1,1), zeros(n + 1,1), zeros(n + 1,1), zeros(n + 1,1)};
for set = 1:4
    [~, ~, V] = svd(M_cell{set});
    sz = size(V,2);
    d{set}(:,1) = V(:,sz);
end

%% Find the correct (R,T) 
%by checking lambda2*x2 = lambda1*R*x1 + gamma*T for
%both lambda1, lambda2 > 0
counter = {0, 0, 0, 0};
for set = 1:4
    for i = 1:n
        if(d{set}(i) > 0)
            counter{set} = counter{set} + 1;
        end
    end
end

index = 1;
value = 0;
for set = 1:4
    if(counter{set} > value)
        index = set;
        value = counter{set};
    end
end

T = d{index}(n + 1)*T_cell{index}; %T = gamma*Ti
R = R_cell{index};
lambdas = d{index}(1:n);

P1 = zeros(3, n);
for i = 1:n
    P1(:,i) = lambdas(i,1)*x1(:,i);
end
close all
% 3D Plot
figure;
scatter3(P1(1,:),P1(2,:),P1(3,:),'bo');
hold on
O = -R'*T; % coordinates of camera 2 in coordinate system 1
scatter3(O(1),O(2),O(3),'rd');
scatter3(0,0,0,'rd');
hold off

end

function [sm] = skew(vec)
sm = [0 -vec(3) vec(2);
    vec(3) 0 -vec(1);
    -vec(2) vec(1) 0];
end