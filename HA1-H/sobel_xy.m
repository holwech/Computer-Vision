function [Fx, Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zurückgibt.

%Sobel filter in x- an y-direction
SFx = [1 0 -1; 2 0 -2; 1 0 -1];
SFy = [1 2 1; 0 0 0; -1 -2 -1];

Fx = zeros(size(Image,1), size(Image,2));
Fy = zeros(size(Image,1), size(Image,2));

intImage = double(Image);

for x = 2:size(Image,2) - 1
    for y = 2:size(Image,1) - 1
        accuX = 0.0;
        accuY = 0.0;
        %disp(['x', num2str(x), 'y', num2str(y)]);
        for k = -1:1
            for l = -1:1
                accuX = accuX + intImage(y - l, x - k)*SFx(2 - l, 2 - k);
                accuY = accuY + intImage(y - l, x - k)*SFy(2 - l, 2 - k);
                
            end
        end
        %disp(['Store in:',' x ', num2str(x), ' y ', num2str(y)]);
        Fx(y,x) = ceil(accuX);
        Fy(y,x) = ceil(accuY);
        %Fx(y,x) = ceil(sqrt(accuX*accuX));
        %Fy(y,x) = ceil(sqrt(accuY*accuY));
    end
end

