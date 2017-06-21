function  [Merkmale] = harris_detektor(Image,varargin) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert

% initial values 
do_plot = false;
segment_length=3;
k=0.05;
tau=1;
min_distance = 1;
tile_size = 21;
N = 1;

% read optional arguments
n=1;
while n+1<=nargin
    switch varargin{n}
        case 'do_plot'
            do_plot = varargin{n+1};
        case 'segment_length'
            segment_length = varargin{n+1};
            if mod(segment_length,2) == 0
                disp 'Error: segment length is an even number! we should never got here';
                return;
            end
        case 'k'
            k = varargin{n+1};
        case 'tau'
            tau = varargin{n+1};
        case 'min_distance'
            min_distance = varargin{n+1};
        case 'tile_size'
            tile_size = varargin{n+1};
        case 'N'
            N = varargin{n+1};
        otherwise
            disp 'Error: unknown input in harris detector! we should never got here'
    end
    n= n+2;
end
fprintf('segment_length = %d, k = %d, tau = %d, do_plot = %d\n',segment_length,k,tau,do_plot);

% check whether the picture does only contain one colour dimension
if size(Image,3) ~= 1
    disp('Error: Picture is not gray!')
    return;
end

% preallocate memory
Merkmale = zeros(size(Image));
H = zeros(size(Image));
G = zeros(2);
[Fx,Fy]=sobel_xy(Image); %gradient
disp 'computing H for all pixels ...'
    p=(segment_length-1)/2;
    y=1+p;
    while(y<size(Image,2)-p)
        x=1+p;
        while(x<size(Image,1)-p)
            % Compute the approximated Harris-Matrix for the Pixel Image(x,y)
            G = zeros(2);
            for ix = -p:p
                for iy = -p:p
                    w = weight(ix, iy, segment_length);
                    G = G + w*[Fx(x+ix,y+iy),Fy(x+ix,y+iy)]'*[Fx(x+ix,y+iy),Fy(x+ix,y+iy)];
                end
            end
            H(x,y) = det(G)-k*(trace(G)^2);
            
            if H(x,y)>tau
                Merkmale(x,y)=1;% Ecke (black)
            %elseif H(x,y)<-tau
                %Merkmale(x,y)=2;% Kante
            %else
                %Merkmale(x,y)=3;% Flaeche
            end
            x=x+1;
        end
        y=y+1;
    end
% W = [0 1 0; 1 2 1; 0 1 0]; %weights

% p = (segment_length-1)/2; %segment bounds
% %fprintf('p=%d\n',p);
% for x = 1+p : size(Image,1)-p 
%     for y = 1+p : size(Image,2)-p
%         G(:,:) = 0;
%         for ix = -p : p
%             for iy = -p : p
%                 dI=[Fx(x+ix, y+iy),Fy(x+ix, y+iy)]'*[Fx(x+ix, y+iy),Fy(x+ix, y+iy)];
%                 G = G + W(ix+p+1, iy+p+1) * dI;
%             end
%         end
%         H(x,y) = det(G) - k*trace(G)^2;
%          if H(x,y) > tau
%              Merkmale(x,y) = 255;
%          end
%     end
% end
disp 'computation of H finished!';

if do_plot
    figure;
    subplot(1,3,1);
    imshow(Image);
    title('original image');
    subplot(1,3,2);
    imshow(Merkmale);
    title('all edges');
end
Merkmale(:,:)=0;

disp 'find local maxima in blocks...';
tile_p = floor(tile_size/2); %quadratic tile with equal tile_size

%IMPLEMENT: Improve the function so that it handles the case when the dimensions of the image
%do not fit the dimentions of the tiles
%Ex: picture width 2000, tile width 21, only room for 95 tiles => 5 pixel
%border not handeled
for y = 1+tile_p : tile_p*2+1 : size(H,1)-tile_p-1
    for x = 1+tile_p : tile_p*2+1 : size(H,2)-tile_p-1
        ys = y - tile_p;
        ye = y + tile_p;
        xs = x - tile_p;
        xe = x + tile_p;
        tile_matrix = H(ys:ye, xs:xe);
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
            hit = checkVacinity(ySorted, xSorted, size(Merkmale), Merkmale, min_distance);
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

if do_plot
    subplot(1,3,3);
    imshow(Merkmale);
    title('only best edges');
end
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
end


function w = weight(ix,iy,segment_length)
    % Computes the Weight, such that pixels that are closer to the the
    % central pixel (i.e. ix and iy are small) have a higher weight
    if (segment_length == 3)
    switch abs(ix*iy)
        case 0
            w = 4;
        case 1
            w = 2;
        case 2
            w = 1;
        otherwise
            disp('error');
    end
    end
end