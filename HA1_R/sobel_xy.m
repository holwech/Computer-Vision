function [Fx,Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zurueckgibt.
n = 1;

% Berechne Normierungsfaktor C
sigma_quadr = 1/2*log(2);
C = 1;
for i = -n:n
    C = C + exp(-(i^2)/(2*sigma_quadr));
end
C = 1/C;
% Berechne Sobelfilter in x-Richtung
Sx = zeros(2*n+1,2*n+1);
for k = -n:n
    for l = -n:n
        Sx(k+n+1,l+n+1) = C*k/(sigma_quadr)*exp(-(k^2+l^2)/(2*sigma_quadr));
    end
end

% Berechne Sobelfilter in y-Richtung
Sy = zeros(2*n+1,2*n+1);
for k = -n:n
    for l = -n:n
        Sy(k+n+1,l+n+1) = C*l/(sigma_quadr)*exp(-(k^2+l^2)/(2*sigma_quadr));
    end
end

% Berechne Faltung mit Sobelfilter
Fx = conv2(Sx,Image);
Fy = conv2(Sy,Image);
    
end

