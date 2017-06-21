function [Fx, Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zur�ckgibt.

% %Sobel filter in x- an y-direction
% SFx = [
%        1  0 -1;
%        2  0 -2;
%        1  0 -1
%       ];
% SFy = [
%        1  2  1;
%        0  0  0;
%       -1 -2 -1
%       ];
% gaussianBlur = 1/16*[
%                      1  2  1;
%                      2  4  2;
%                      1  2  1
%                    ];
% 
% imgdim = size(Image);
% Fx = zeros(imgdim(1), imgdim(2));
% Fy = zeros(imgdim(1), imgdim(2));
% 
% % Add padding to image
% padding = floor(size(SFx, 1) / 2);
% paddedImage = double(padarray(Image, [padding padding], 'symmetric'));
% 
% 
% 
% %Flip kernels for convolution
% SFx = rot90(SFx,2);
% SFy = rot90(SFy,2);
% 
% for x = (padding + 1):(imgdim(2) + padding)
%     for y = (padding + 1):(imgdim(1) + padding)
%         accuX = 0.0;
%         accuY = 0.0;
%         for k = -padding:padding
%             for l = -padding:padding
%                 accuX = accuX + paddedImage(y - l, x - k)*SFx(2 - l, 2 - k) * gaussianBlur(2 - l, 2 - k);
%                 accuY = accuY + paddedImage(y - l, x - k)*SFy(2 - l, 2 - k) * gaussianBlur(2 - l, 2 - k);
%                 
%             end
%         end
%         Fx(y - padding,x - padding) = ceil(accuX);
%         Fy(y - padding,x - padding) = ceil(accuY);
%     end
% end

p = 1;

% Berechne Normierungsfaktor C
sigma_quadr = 1/2*log(2);
C = 1;
for i = -p:p
    C = C + exp(-(i^2)/(2*sigma_quadr));
end
C = 1/C;
% Berechne Sobelfilter in x-Richtung
Sx = zeros(2*p+1,2*p+1);
for k = -p:p
    for l = -p:p
        Sx(k+p+1,l+p+1) = C*k/(sigma_quadr)*exp(-(k^2+l^2)/(2*sigma_quadr));
    end
end

% Berechne Sobelfilter in y-Richtung
Sy = zeros(2*p+1,2*p+1);
for k = -p:p
    for l = -p:p
        Sy(k+p+1,l+p+1) = C*l/(sigma_quadr)*exp(-(k^2+l^2)/(2*sigma_quadr));
    end
end

% Berechne Faltung mit Sobelfilter
Fx = conv2(Sx,Image);
Fy = conv2(Sy,Image);
end