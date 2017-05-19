padding = floor(segment_length/2);
%W = 1/16*[1 2 1; 2 4 2; 1 2 1]; %Gaussian blur
W = [0 1 0; 1 2 1; 0 1 0];
Merkmale = zeros(size(Image,1), size(Image,2));
Hv = zeros(size(Image,1), size(Image,2));
%For each pixel
for y = padding + 1:size(Image,1) - padding 
    for x = padding + 1:size(Image,2) - padding
        %Create window
        G = zeros(2);
        for k = 0:segment_length - 1
            for l = 0:segment_length - 1
                
                Ixx = Fx(y - padding + l, x - padding + k)*Fx(y - padding + l, x - padding + k);
                Iyy = Fy(y - padding + l, x - padding + k)*Fy(y - padding + l, x - padding + k);
                Ixy = Fx(y - padding + l, x - padding + k)*Fy(y - padding + l, x - padding + k);
                dI = [Ixx Ixy; Ixy Iyy];
                %disp(dI)
                G = G + W(l + 1,k + 1)*(dI);
                
            end
        end
        H = det(G) - k*trace(G)*trace(G);
        Hv(y,x) = H;
        if H > tau %|| H < -tau
            Merkmale(y,x) = 1;
        end
    end
end

black = cat(3, zeros(size(Image)), zeros(size(Image)), zeros(size(Image)));
for y = 1:size(Merkmale,1)
    for x = 1:size(Merkmale,2)
        if Merkmale(y,x) == 1
            black(y,x,1) = 255;
            black(y,x,2) = 0;
            black(y,x,3) = 0;
        end
    end
end
imshow(black);
% red = cat(3, ones(size(Image)), zeros(size(Image)), zeros(size(Image)));
% hold on
% h = imshow(red);
% hold off
% set(h, 'AlphaData', Image);