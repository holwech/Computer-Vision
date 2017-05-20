%function [Fx, Fy] = sobel_xy(Image)
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



%function  [Merkmale] = harris_detektor(Image, Fx, Fy, segment_length, k, tau, do_plot) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert

padding = floor(segment_length/2);
W = 1/16*[1 2 1; 2 4 2; 1 2 1]; %Gaussian blur
%W = [0 1 0; 1 2 1; 0 1 0];
Merkmale = zeros(size(Image,1), size(Image,2));
Hv = zeros(size(Image,1), size(Image,2));
%For each pixel
%Create window
G = zeros(2);
dI = zeros(2);
for y = 1 + padding:size(Image,1) - padding 
    for x = 1 + padding:size(Image,2) - padding
        G(:,:) = 0;
        for k = 0:segment_length - 1
            for l = 0:segment_length - 1
                
                Ixx = Fx(y - padding + l, x - padding + k)*Fx(y - padding + l, x - padding + k);
                Iyy = Fy(y - padding + l, x - padding + k)*Fy(y - padding + l, x - padding + k);
                Ixy = Fx(y - padding + l, x - padding + k)*Fy(y - padding + l, x - padding + k);
                dI(1,1) = Ixx;
                dI(1,2) = Ixy;
                dI(2,1) = Ixy;
                dI(2,2) = Iyy;
                G = G + W(l + 1,k + 1)*(dI);
                
            end
        end
        %G = single(G);
        H = det(G) - 0.05*(trace(G))^2;
        Hv(y,x) = H;
    end
end

%Find largest values in blocks
featureCandidates = zeros(size(Hv));
for y = 1 + padding:segment_length:size(Image,1) - padding 
    for x = 1 + padding:segment_length:size(Image,2) - padding
        max = -Inf;
        ymax = 0;
        xmax = 0;
        for k = 0:segment_length - 1
            for l = 0:segment_length - 1
                feature = Hv(y - padding + l, x - padding + k);
                if feature > max
                    max = feature;
                    ymax = y - padding + l;
                    xmax = x - padding + k;
                end 
            end
        end
        featureCandidates(ymax, xmax) = Hv(ymax, xmax);
    end
end

%Filter out N maximum H values
N = 1000;
[val , index] = sort(featureCandidates(:), 'descend');
for point = 1:N
    y = index(point) - size(Hv,1)*floor(index(point)/size(Hv,1));
    x = ceil(index(point)/size(Hv,1));
    if val(point) > tau 
        Merkmale(y,x) = 2;
    end
end
                

black = cat(3, zeros(size(Image)), zeros(size(Image)), zeros(size(Image)));
for y = 1:size(Merkmale,1)
    for x = 1:size(Merkmale,2)
        if Merkmale(y,x) == 2
            for k = 0:segment_length - 1
                for l = 0:segment_length - 1
                    Image(y - padding + l, x - padding + k,1) = 255;
                end
            end
        end
    end
end

if do_plot
    imshow(Image)
end
% I = imread(black);
% green = cat(3, zeros(size(Image)), ones(size(Image)), zeros(size(Image)));
% hold on
% h = imshow(green);
% hold off
% set(h, 'AlphaData', I);