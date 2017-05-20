function  [Hvalues, Merkmale] = harris_detektor(Image, Fx, Fy, segment_length, k, tau, do_plot, min_distance, tile_size, N) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert

padding = floor(segment_length/2);
W = [0 1 0; 1 2 1; 0 1 0];
Merkmale = zeros(size(Image,1), size(Image,2));
Hvalues = zeros(size(Image,1), size(Image,2));
G = zeros(2);
dI = zeros(2);
traceWeight = k; %For some reason the expression for H does not work with k directly
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
        H = det(G) - traceWeight*power(trace(G),2);
        Hvalues(y,x) = H;
%         if H > 10000
%             Merkmale(y,x) = 255;
%         end
    end
end

%Find local maxima in blocks
tile_padding = floor(tile_size/2); %quadratic tile with with equal tile_size
max = -Inf;
ymax = 1;
xmax = 1;
draw_padding = 1; %quadratic square with width (2*draw_padding + 1)
tile_matrix = zeros(tile_padding*2 + 1);
placed_counter = 0;
length_counter = 1;
sz = size(Merkmale);
%IMPLEMENT: Improve the function so that it handles the case when the dimentions of the image
%does not fit the dimentions of the tiles
%Ex: picture width 2000, tile width 21, only room for 95 tiles => 5 pixel
%border not handeled
for y = 1 + tile_padding:(tile_padding*2 + 1):size(Hvalues,1) - tile_padding - 1
    for x = 1 + tile_padding:(tile_padding*2 + 1):size(Hvalues,2) - tile_padding - 1
        ys = y - tile_padding;
        ye = y + tile_padding;
        xs = x - tile_padding;
        xe = x + tile_padding;
        tile_matrix = Hvalues(ys:ye, xs:xe);
        [values, index] = sort(tile_matrix(:), 'descend');
        placed_counter = 1;
        length_counter = 1;
        while (placed_counter < N && length_counter < length(values))
            %Find 2D x, y coordinates relative the tile_matrix
            ySorted = index(length_counter) - size(tile_matrix,1)*(ceil(index(length_counter)/size(tile_matrix,1)) - 1);
            xSorted = ceil(index(length_counter)/size(tile_matrix,1));
            %Shift the 2D coordinates so they concure with the Hvalues matrix
            ySorted = ySorted + (ys - 1);
            xSorted = xSorted + (xs - 1);
            hit = checkVacinity(ySorted, xSorted, sz, Merkmale, min_distance);
            if (hit == 0 && values(length_counter) > 10000)
                %Draw a 3x3 square
                %Merkmale((ySorted-draw_padding:ySorted + draw_padding), (xSorted - draw_padding: xSorted + draw_padding)) = 255;
                Merkmale(ySorted, xSorted) = 255;
                placed_counter = placed_counter + 1;
            end
            length_counter = length_counter + 1;
        end
    end
end
imshow(Merkmale);
            
if do_plot
    imshow(Image)
    red = cat(3, ones(size(Image)), zeros(size(Image)), zeros(size(Image)));
    hold on
    h = imshow(red);
    hold off
    set(h, 'AlphaData', Merkmale);
end

%Function to check if there is already a discovered feature in min_distance
%from the newly discovered feature
%IMPLEMENT: Make the function work in areas close to the image boarder
function [hit] = checkVacinity(y, x, sz, Merkmale, min_distance)
%Case: The point (y,x) is min_distance pixels or more from the boarder
hit = 0;
if ((y > min_distance) && (y < sz(1) - min_distance) && (x > min_distance) && (x < sz(2) - min_distance))
    for k = 0:2*min_distance + 1
        for l = 0:2*min_distance + 1
            val = Merkmale(y - min_distance + l, x - min_distance + k);
            if (val == 255)
                hit = 1;
            end
        end
    end
end
            
            
                      
