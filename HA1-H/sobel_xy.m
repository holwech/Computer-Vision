function [Fx, Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zur�ckgibt.

%Sobel filter in x- an y-direction
SFx = [1 0 -1; 2 0 -2; 1 0 -1];
SFy = [1 2 1; 0 0 0; -1 -2 -1];
gaussianBlur = 1/16*[1 2 1; 2 4 2; 1 2 1];

Fx = zeros(size(Image,1), size(Image,2));
Fy = zeros(size(Image,1), size(Image,2));

%Convert to double to do calculations
doubleImage = double(Image);

%Flip kernels for convolution
SFx = rot90(SFx,2);
SFy = rot90(SFy,2);

for x = 2:size(Image,2) - 1
    for y = 2:size(Image,1) - 1
        accuX = 0.0;
        accuY = 0.0;
        for k = -1:1
            for l = -1:1
                accuX = accuX + doubleImage(y - l, x - k)*SFx(2 - l, 2 - k)*gaussianBlur(2 - l, 2 - k);
                accuY = accuY + doubleImage(y - l, x - k)*SFy(2 - l, 2 - k)*gaussianBlur(2 - l, 2 - k);
                
            end
        end
        Fx(y,x) = ceil(accuX);
        Fy(y,x) = ceil(accuY);
    end
end

