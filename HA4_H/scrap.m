% function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% % Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% % Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1
% 
% %Want to solve the equation [x2_hat*R*x1 x2_hat*T][lambda; gamma] = 0
% % => Md = 0, by minimization of Md with respect to d and norm(d) = 1;
% 
% n = size(Korrespondenzen, 2);
% 
% %% Preprocess the KPs
% %Put the pixel coordinates in homogenous form and use the calibration
% %matrix K to find the homogenous coordinates.
% x1 = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
% x2 = [Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];
% x1 = K\x1;
% x2 = K\x2;
% 
% %% Constuct M for both sets of (R,T):
% diagonal_elements_1 = zeros(3,n);
% last_column_1 = zeros(3,n);
% diagonal_elements_2 = zeros(3,n);
% last_column_2 = zeros(3,n);
% 
% %Calculate x2_hat*R*x1 and x2_hat*T for all KP
% for i = 1:n
%     x2_hat = skew(x2(:,i));
%     diagonal_elements_1(:,i) = x2_hat*R1*x1(:,i);
%     diagonal_elements_2(:,i) = x2_hat*R2*x1(:,i);
%     last_column_1(:,i) = x2_hat*T1;
%     last_column_2(:,i) = x2_hat*T2;
% end
% 
% %Gather in matrix form
% [r,c] = size(diagonal_elements_1);
% i = 1:numel(diagonal_elements_1);
% j = repmat(1:c,r,1);
% M1_diag = sparse(i',j(:),diagonal_elements_1(:));
% M2_diag = sparse(i',j(:),diagonal_elements_2(:));
% M1 = [M1_diag reshape(last_column_1,3*n, 1)];
% M2 = [M2_diag reshape(last_column_2,3*n, 1)];
% M1 = full(M1);
% M2 = full(M2);
% 
% %% Minimize M*d using svd
% [~, ~, V1] = svd(M1);
% sz = size(V1,2);
% d1 = V1(:,sz);
% 
% [~, ~, V2] = svd(M2);
% sz = size(V2,2);
% d2 = V2(:,sz);
% 
% %% Find lambda2
% lambda_M1 = zeros(1,n);
% lambda_M2 = zeros(1,n);
% for i = 1:n
%     x1_M1_estimate = d1(i)*R1*x1(:,i) + d1(n + 1)*T1;
%     x1_M2_estimate = d2(i)*R2*x1(:,i) + d2(n + 1)*T2;
%     %We assume errors in our system leads to different lambda2 values for
%     %the different entries of x1, and thus take the mean
%     lambda_M1(i) = mean(x1_M1_estimate./x2(:,i));
%     lambda_M2(i) = mean(x1_M2_estimate./x2(:,i));
% end
% % disp(x1_KP(:,1:5))
% % disp(x2_KP(:,1:5))
% % disp(lambda_M1(:,1:5))
% 
% %% Find the correct (R,T) 
% %by checking lambda2*x2 = lambda1*R*x1 + gamma*T for
% %both lambda1, lambda2 > 0
% counter_M1 = 0;
% counter_M2 = 0;
% for i = 1:n
%     if(d1(i) > 0 && lambda_M1(i) > 0)
%         counter_M1 = counter_M1 + 1;
%     end
%     
%     if(d2(i) > 0 && lambda_M2(i) > 0)
%         counter_M2 = counter_M2 + 1;
%     end
% end
% 
% if(counter_M1 > counter_M2)
%     T = T1;
%     R = R1;
%     lambdas = [d1(:)'; lambda_M1 0];
%     P1 = d1(1,:)*x1(:,:);
% else
%     T = T2;
%     R = R2;
%     lambdas = [d2(:)'; lambda_M2 0];
%     P1 = d2(1,:)*x1(:,:);
% end
% 
% disp(counter_M1)
% disp(counter_M2)
% 
% end
% 
% function [sm] = skew(vec)
% sm = [0 -vec(3) vec(2);
%     vec(3) 0 -vec(1);
%     -vec(2) vec(1) 0];
% end



%% Find lambda2
lambda = {zeros(1,n), zeros(1,n), zeros(1,n), zeros(1,n)};
for set = 1:4
    for i = 1:n
        x1_estimate = d{set}(i)*R_cell{set}*x1(:,i) + d{set}(n + 1)*T_cell{set};
        %We assume errors in our system leads to different lambda2 values for
        %the different entries of x1, and thus take the mean
        lambda{set}(1,i) = mean(x1_estimate./x2(:,i));
    end
end
% disp(x1_KP(:,1:5))
% disp(x2_KP(:,1:5))
% disp(lambda_M1(:,1:5))

%% Find the correct (R,T) 
%by checking lambda2*x2 = lambda1*R*x1 + gamma*T for
%both lambda1, lambda2 > 0
counter = {0, 0, 0, 0};
for set = 1:4
    for i = 1:n
        if(d{set}(i) > 0 && lambda{set}(i) > 0)
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
lambdas = [d{index}(:)'; lambda{index} 0];
P1 = d{index}(1,:)*x1(:,:);



%19.07
function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

%Want to solve the equation [x2_hat*R*x1 x2_hat*T][lambda; gamma] = 0
% => Md = 0, by minimization of Md with respect to d and norm(d) = 1;

n = size(Korrespondenzen, 2);

%% Preprocess the KPs
%Put the pixel coordinates in homogenous form and use the calibration
%matrix K to find the homogenous coordinates.
x1 = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
x2 = [Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];
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
Pi_0 = [1 0 0 0;
        0 1 0 0;
        0 0 1 0];
for i = 1:n
    P1(:,i) = lambdas(i,1)*x1(:,i);
end
P1 = [P1(:,:); ones(1,n)];

end

function [sm] = skew(vec)
sm = [0 -vec(3) vec(2);
    vec(3) 0 -vec(1);
    -vec(2) vec(1) 0];
end


%% ------

function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

%Want to solve the equation [x2_hat*R*x1 x2_hat*T][lambda; gamma] = 0
% => Md = 0, by minimization of Md with respect to d and norm(d) = 1;

n = size(Korrespondenzen, 2);

%% Preprocess the KPs
%Put the pixel coordinates in homogenous form and use the calibration
%matrix K to find the homogenous coordinates.
x1 = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
x2 = [Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];
x1 = K\x1;
x2 = K\x2;

%% Constuct M for all four possible (R,T):

%We have four combinations of (R1,R2,T1,T2)
%(R1,T1) ==> E
%(R2,T2)
%(R1,T2) ==> -E
%(R2,T1)

% diagonal_elements = {zeros(3,n), zeros(3,n), zeros(3,n), zeros(3,n)};
% column_elements = {zeros(3,n), zeros(3,n), zeros(3,n), zeros(3,n)};
% R_cell = {R1, R1, R2, R2};
% T_cell = {T1, T2, T1, T2};
% M_diag_cell = {0, 0, 0, 0};
% M_cell = {0, 0, 0, 0};

%Calculate x2_hat*R*x1 and x2_hat*T for all KP
% for set = 1:4
%     for i = 1:n
%         x2_hat = skew(x2(:,i));
%         diagonal_elements{set}(:,i) = x2_hat*R_cell{set}*x1(:,i);
%         column_elements{set}(:,i) = x2_hat*T_cell{set};
%     end
% 
%     %Gather in matrix form
%     [r,c] = size(diagonal_elements{set});
%     i = 1:numel(diagonal_elements{set});
%     j = repmat(1:c,r,1);
%     M_diag_cell{set} = full(sparse(i',j(:),diagonal_elements{set}(:)));
%     M_cell{set} = [M_diag_cell{set} reshape(column_elements{set},3*n, 1)];
%     %M_cell{set} = full(M_cell{set});
%     
% end

R_cell = {R1, R2, R2, R1};
T_cell = {T1, T2, T1, T2};
M_diag = [zeros(3*n,n,4)];
M_column = [zeros(3*n, 1, 4)];
M = [zeros(3*n, n + 1, 4)];
for set = 1:4
    j = 1;
    for i = 1:n
        x2_hat = skew(x2(:,i));
        M_diag(j:j+2, i, set) = x2_hat*R_cell{set}*x1(:,i);
        M_column(j:j+2, 1, set) = x2_hat*T_cell{set};
        j = j + 3;
    end
    M(:,:,set) = [M_diag(:,:,set) M_column(:,1,set)];
end
        
    



%% Minimize M*d using svd
d = [zeros(n + 1, 4)];
for set = 1:4
    [~, ~, V] = svd(M(:,:,set));
    sz = size(V,2);
    d(:,set) = V(:,sz);
end


%% Find the correct (R,T) 
%by checking lambda2*x2 = lambda1*R*x1 + gamma*T for
%both lambda1, lambda2 > 0
counter = [0, 0, 0, 0];
for set = 1:4
    for i = 1:n
        if(d(i, set) > 0)
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


T = d(n + 1, index)*T_cell{index}; %T = gamma*Ti
R = R_cell{index};
lambdas = d(1:n, index);

P1 = zeros(3, n);
for i = 1:n
    P1(:,i) = lambdas(i,1)*x1(:,i);
end
P1 = [P1(:,:); ones(1,n)];

end

function [sm] = skew(vec)
sm = [0 -vec(3) vec(2);
    vec(3) 0 -vec(1);
    -vec(2) vec(1) 0];
 end
