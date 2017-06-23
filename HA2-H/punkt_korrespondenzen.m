function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
%% Set and read input variables
P = inputParser;

P.addOptional('window_length', 15, @isnumeric);
P.addOptional('do_plot', false, @islogical);
P.addOptional('min_corr', 0.9, @isnumeric);

P.parse(varargin{:});
window_length  = P.Results.window_length;
do_plot         = P.Results.do_plot;
min_corr        = P.Results.min_corr;

%We use windows with odd number width, so we have symmetry around the
%central pixel.
if mod(window_length, 2) == 0
    window_length = window_length + 1;
end

%% Allocate memory, calculate variables
dim_mpt1 = size(Mpt1);
dim_mpt2 = size(Mpt2);
I1_mean = mean2(I1);
I2_mean = mean2(I2);
V_mean = double(ones(window_length)*I1_mean);
W_mean = double(ones(window_length)*I2_mean);
matches = zeros(2, dim_mpt1(2));
max_ncc = NaN(1, dim_mpt1(2));

%% Image padding
% Padding is added to the images, in case of comparisons of points close to the
% the edge.
pdg = floor(window_length / 2);
twoPdg = 2 * pdg; % Speed optimization
paddedI1 = padarray(I1, [pdg pdg], 'symmetric');
paddedI2 = padarray(I2, [pdg pdg], 'symmetric');

%% Match search
% Search through each keypoint and find the point with max correlation
for p1 = 1:dim_mpt1(2)
    p1_seg = paddedI1(...
        Mpt1(2, p1):(Mpt1(2, p1) + twoPdg),  ...
        Mpt1(1, p1):(Mpt1(1, p1) + twoPdg) ...
    );
    % If there is no matching keypoint, the index will be stored
    for p2 = 1:dim_mpt2(2)
        p2_seg = paddedI2(...
            Mpt2(2, p2):(Mpt2(2, p2) + twoPdg),  ...            
            Mpt2(1, p2):(Mpt2(1, p2) + twoPdg) ...
        );
        % Calculate normalized cross correlation value
        ncc = NCC(p1_seg, p2_seg, V_mean, W_mean);
        % If NCC is higher than previous matches and over threshold, store
        % the values in matches and max_ncc
        if ((ncc > max_ncc(p1) || isnan(max_ncc(p1))) && (ncc > min_corr))
            max_ncc(p1) = ncc;
            matches(1:2, p1) = Mpt2(:, p2);
        end
    end
end

%% Filter non-matches
% Remove all points that do not have any matches
% Using the fact that max_ncc indices will be NaN if no match is found for
% that point. Take the inverse of isnan, and get the matrix of all matches.
Korrespondenzen = [Mpt1(:, ~isnan(max_ncc) == 1); matches(:, ~isnan(max_ncc) == 1)];
fprintf('Number of matches found: %i \n', size(Korrespondenzen,2));

%% Plotting
if(do_plot)
    %Plot feature points as blue dots and feature points that have a match
    %in red.
    imshow([I1 I2]);
    sx = size(I1,2);
    hold on;
    plot(Mpt1(1, :), Mpt1(2, :), 'b.');
    plot(Korrespondenzen(1, :), Korrespondenzen(2, :), 'r.');
    plot(Mpt2(1, :) + sx, Mpt2(2, :), 'b.');
    plot(Korrespondenzen(3, :) + sx, Korrespondenzen(4, :), 'r.');
    
    %Plot lines between the matching points
%     sx = size(I1,2);
%     imshow([I1 I2]);
%     x1 = Korrespondenzen(1,:)'; y1 = Korrespondenzen(2,:)';
%     x2 = Korrespondenzen(3,:)'; y2 = Korrespondenzen(4,:)';
%     hold on;
%     x2 = x2 + sx;
%     plot([x1 x2]', [y1 y2]', 'r', 'LineWidth', 1);
    
end
end

