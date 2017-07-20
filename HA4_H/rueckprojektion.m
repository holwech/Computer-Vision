function repro_error = rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler

n = size(Korrespondenzen,2);
%Put the pixel coordinates in homogenous form and use the calibration
%matrix K to find the homogenous coordinates.
x1 = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
x2 = [Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];
x1 = K\x1;
x2 = K\x2;

% Find projected points
P2 = [zeros(3,n); ones(1,n)];
x2_est = zeros(3,n);
M = [R T; zeros(1,3) 1];

for i = 1:n
    P2(:,i) = M*P1(:,i);
    %P2 = lambda2*x2 ==> find the lambda that sets the z-coord of x2_est to
    %one
    x2_est(:,i) = P2(1:3,i)/P2(3,i);
end

%Compute backprojection error

%Transform to pixel coordinates
%disp(x2_est)
%disp(x2)
x2 = K*x2;
x2_est = K*x2_est;

repro_error = 0;
for i = 1:n
    repro_error = repro_error + norm(x2_est(:,i) - x2(:,i));
end
repro_error = repro_error/n;

figure();
imshow(uint8(I2));
hold on;
%Plot estimates
plot(x2_est(1,:), x2_est(2,:), '*r');
%Plot KPs
plot(x2(1,:), x2(2,:), '*g');
title('Green points are the KPs, and red points are the projected points');


end