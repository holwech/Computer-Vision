function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

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

R_cell = {R1, R2, R2, R1};
T_cell = {T1, T2, T1, T2};
%For finding lambda_1s
M_diag = [zeros(3*n,n,4)];
M_column = [zeros(3*n, 1, 4)];
M = [zeros(3*n, n + 1, 4)];
%For finding lambda_2s
M2_diag = [zeros(3*n,n,4)];
M2_column = [zeros(3*n, 1, 4)];
M2 = [zeros(3*n, n + 1, 4)];
for set = 1:4
    j = 1;
    for i = 1:n
        x2_hat = skew(x2(:,i));
        x1_hat = skew(x1(:,i));
        
        M_diag(j:j+2, i, set) = x2_hat*R_cell{set}*x1(:,i);
        M_column(j:j+2, 1, set) = x2_hat*T_cell{set};
        
        M2_diag(j:j+2, i, set) = x1_hat*R_cell{set}'*x2(:,i);
        M2_column(j:j+2, 1, set) = -x1_hat*R_cell{set}'*T_cell{set};
        
        j = j + 3;
    end
    M(:,:,set) = [M_diag(:,:,set) M_column(:,1,set)];
    M2(:,:,set) = [M2_diag(:,:,set) M2_column(:,1,set)];
end
        
    



%% Minimize M*d using svd
d1 = [zeros(n + 1, 4)];
d2 = [zeros(n + 1, 4)];
for set = 1:4
    [~, ~, V] = svd(M(:,:,set));
    sz = size(V,2);
    d1(:,set) = V(:,sz);
    d1(:,set) = d1(:,set)/d1(n + 1, set);
    
    [~, ~, V2] = svd(M2(:,:,set));
    sz = size(V2,2);
    d2(:,set) = V2(:,sz);
    d2(:,set) = d2(:,set)/d2(n + 1, set);
end


%% Find the correct (R,T) 
%by checking lambda2*x2 = lambda1*R*x1 + gamma*T for
%both lambda1, lambda2 > 0
counter = [zeros(1,4)];
for set = 1:4
    for i = 1:n
        if(d1(i, set) > 0 && d2(i, set) > 0)
            counter(set) = counter(set) + 1;
        end
    end
end

index = 1;
value = 0;
for set = 1:4
    if(counter(set) > value)
        index = set;
        value = counter(set);
    end
end

T = T_cell{index}; %T = gamma*Ti
R = R_cell{index};
lambdas = [d1(1:n, index), d2(1:n, index)];

P1 = zeros(3, n);
for i = 1:n
    P1(:,i) = lambdas(i,1)*x1(:,i);
end

P1 = [P1(:,:); ones(1,n)];



close all
%% 3D Plot w.r.t cameraframe 1
figure;
%Change coordinates for plot so it looks nicer with respect to 
%the way matlab plots it
%Axis change: z_camerea = x_matlab, y_camera = -y_matlab
%To achieve this we rotate 90deg about x and 90deg about z
Rx = [1 0 0; 0 0 -1; 0 1 0]; %Rotate 90deg about x-axis
Rz = [0 -1 0; 1 0 0; 0 0 1]; %Rotate 90deg about z-axis
M = Rx*Rz; % x_camera_coord = M*x_plot_coord, M is the change of basis
%This means that x_plot_coord = M'*x_camera_coord
P_plot = M'*P1(1:3,:);
scatter3(P_plot(1,:),P_plot(2,:),P_plot(3,:),'bo');
hold on
O = -R'*T; % coordinates of camera 2 in coordinate system 1
O_plot = M'*O;
scatter3(O_plot(1),O_plot(2),O_plot(3),'rd');
scatter3(0,0,0,'rs');

xlabel('x') % x-axis label
ylabel('y') % y-axis label
zlabel('z') % z-axis label

hold off;

%% 3D plot w.r.t cameraframe 2
% figure;
% P2 = zeros(3, n);
% for i = 1:n
%     P2(:,i) = lambdas(i,2)*x2(:,i);
% end
% 
% P2 = [P2(:,:); ones(1,n)];
% 
% %Change coordinates for plot so it looks nicer with respect to 
% %the way matlab plots it
% %Axis change: z_camerea = x_matlab, y_camera = -y_matlab
% %To achieve this we rotate 90deg about x and 90deg about z
% Rx = [1 0 0; 0 0 -1; 0 1 0]; %Rotate 90deg about x-axis
% Rz = [0 -1 0; 1 0 0; 0 0 1]; %Rotate 90deg about z-axis
% M = Rx*Rz; % x_camera_coord = M*x_plot_coord, M is the change of basis
% %This means that x_plot_coord = M'*x_camera_coord
% P_plot = M'*P2(1:3,:);
% scatter3(P_plot(1,:),P_plot(2,:),P_plot(3,:),'bo');
% hold on
% O = T; % coordinates of camera 2 in coordinate system 1
% O_plot = M'*O;
% scatter3(O_plot(1),O_plot(2),O_plot(3),'rd');
% scatter3(0,0,0,'rs');
% 
% xlabel('x') % x-axis label
% ylabel('y') % y-axis label
% zlabel('z') % z-axis label
% 
% hold off;

end

function [sm] = skew(vec)
sm = [0 -vec(3) vec(2);
    vec(3) 0 -vec(1);
    -vec(2) vec(1) 0];
 end