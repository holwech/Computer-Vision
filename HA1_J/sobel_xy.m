function [Fx,Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zur�ckgibt.
sobel_x = [
           -1  0  1;
           -2  0  2;
           -1  0  1
          ];
sobel_y = [
           -1 -2 -1;
            0  0  0;
            1  2  1
          ];
imgdim = size(Image);
Image = double(Image);
Fx = zeros(imgdim(1), imgdim(2));
Fy = zeros(imgdim(1), imgdim(2));
for row = 2:(imgdim(1) - 1)
    for col = 2:(imgdim(2) - 1)
        for col2 = 1:3
            for row2 = 1:3
                Fx(row, col) = Fx(row, col) + sobel_x(row2, col2) * Image(row + row2 - 2, col + col2 - 2);
                Fy(row, col) = Fy(row, col) + sobel_y(row2, col2) * Image(row + row2 - 2, col + col2 - 2);
            end
        end
    end
end
end

