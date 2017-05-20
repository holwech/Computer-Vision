%Isolate features
featureCandidates = zeros(size(Hv));
[val , index] = sort(Hv(:), 'descend');
area_size = 20;
for point = 1:length(val)
    %Find 2D-coordinates of possible feature point
    y = index(point) - size(Hv,1)*floor(index(point)/size(Hv,1));
    x = ceil(index(point)/size(Hv,1));
    %Check if a better feature candidate has already been chosen in the
    %area
    taken = 0;
    lLim = area_size;
    rLim = size(Hv,2) - area_size;
    uLim = area_size;
    dLim = size(Hv,1) - area_size;
    for k = 0:2*area_size
        for l = 0:2*area_size
            %make sure we don't get undefined indices
            if (y - area_size + l > uLim && y - area_size + l < dLim && x - area_size + k > lLim && x - area_size + k < rLim)
                if featureCandidates(y - area_size + l, x - area_size + k) == 1
                    taken = 1;
                end
            end
        end
    end
    if taken == 0 && val(point) > tau
        featureCandidates(y,x) = 1;
    end
end

%Draw all featureCandidates
rectangle_edge = 3;
for y = 1 + rectangle_edge:size(featureCandidates,1) - rectangle_edge
    for x = 1 + rectangle_edge:size(featureCandidates,2) - rectangle_edge
        if featureCandidates(y,x) == 1
            %draw a rectangle around feature
            for k = 0:2*rectangle_edge
                for l = 0:2*rectangle_edge
                    Merkmale(y - rectangle_edge + l, x - rectangle_edge + k) = 255;
                end
            end
        end
    end
end


%Filter out N maximum H values
% N = 2000;
% rectangle_edge = 10;
% [val , index] = sort(featureCandidates(:), 'descend');
% for point = 1:N
%     y = index(point) - size(Hv,1)*floor(index(point)/size(Hv,1));
%     x = ceil(index(point)/size(Hv,1));
%     if val(point) > tau
%         %draw a rectangle around feature
%         for k = 0:2*rectangle_edge
%             for l = 0:2*rectangle_edge
%                 Merkmale(y - rectangle_edge + l, x - rectangle_edge + k) = 255;
%             end
%         end
%     end
% end

%imshow(Merkmale);
                
% black = cat(3, zeros(size(Image)), zeros(size(Image)), zeros(size(Image)));
% for y = 1:size(Merkmale,1)
%     for x = 1:size(Merkmale,2)
%         if Merkmale(y,x) == 255
%             for k = 0:segment_length - 1
%                 for l = 0:segment_length - 1
%                     Image(y - padding + l, x - padding + k,1) = 255;
%                 end
%             end
%         end
%     end
% end
% 


%Find local maxima in blocks
tile_padding = floor(tile_size/2); %quadratic tile with with equal tile_size
max = -Inf;
ymax = 1;
xmax = 1;
draw_padding = 1; %quadratic square with width (2*draw_padding + 1)
%IMPLEMENT: Function that handles the case when the dimentions of the image
%does not fit the dimentions of the tiles (there will be an unhandeled
%area)
for y = 1 + tile_padding:(tile_padding*2 + 1):size(Merkmale,1) - tile_padding - 1
    for x = 1 + tile_padding:(tile_padding*2 + 1):size(Merkmale,2) - tile_padding - 1
        max = 0;
        for k = 0:2*tile_padding
            for l = 0:2*tile_padding
                val = Hvalues(y - tile_padding + l, x - tile_padding + k);
                if val > max
                    max = val;
                    ymax = y - tile_padding + l;
                    xmax = x - tile_padding + k;
                end
            end
        end
        if max > 10000 %tau
            Merkmale(ymax, xmax) = 255;
            for k = 0:2*draw_padding
                for l = 0:2*draw_padding
                    Merkmale(ymax - draw_padding + l, xmax - draw_padding + k) = 255;
                end
            end
        end
    end
end
      