function repro_error = rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler
n=size(Korrespondenzen,2);

P2_est = R*P1 + repmat(T,1,n);
P2_est = P2_est * diag(1./P2_est(3,:));
x2_est = K*P2_est;
x2_est = x2_est(1:2,:);

x2=Korrespondenzen(3:4,:);
dx = x2-x2_est;
repro_error=0;
for i=1:n
    repro_error = repro_error + norm(dx(:,i));
end
repro_error = repro_error/n;

figure();
imshow(uint8(I2));
hold on;
%Plot estimates
plot(x2_est(1,:), x2_est(2,:), '*r','DisplayName', 'Projected Points');
%Plot KPs
plot(x2(1,:), x2(2,:), '*g', 'DisplayName', 'KP Points');
for i=1:n
    hold on
    a1 = [x2_est(1,i), x2(1,i)];
    a2 = [x2_est(2,i), x2(2,i)];
    line(a1,a2);
end
title('Green points are the KPs, and red points are the projected points');
hold off;
end