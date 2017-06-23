function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.

%% Set and read input variables
P = inputParser;

P.addOptional('window_length', 40, @isnumeric);
P.addOptional('do_plot', false, @islogical);
P.addOptional('min_corr', 0.9, @isnumeric);

P.parse(varargin{:});
window_length  = P.Results.window_length;
do_plot         = P.Results.do_plot;
min_corr        = P.Results.min_corr;


%% Allocate memory calculate variables
dim_mpt1 = size(Mpt1);
dim_mpt2 = size(Mpt2);
I1_mean = mean2(I1);
I2_mean = mean2(I2);
matches = zeros(2, dim_mpt1(2));
max_ncc = NaN(1, dim_mpt1(2));
seg_area = (2 * floor(window_length / 2) + 1) ^ 2;

%% Image padding
% Added padding to the images, in case of comparisons of points close to the
% the edge.
pdg = floor(window_length / 2);
twoPdg = 2 * pdg; % Speed optimization
paddedI1 = padarray(I1, [pdg pdg], 'symmetric') - I1_mean;
paddedI2 = padarray(I2, [pdg pdg], 'symmetric') - I2_mean;

%% Match search
% Search through each keypoint and find the point with max correlation
for p1 = 1:dim_mpt1(2)
    p1
    p1_seg = double(reshape(paddedI1(...
        Mpt1(2, p1):(Mpt1(2, p1) + twoPdg),  ...
        Mpt1(1, p1):(Mpt1(1, p1) + twoPdg) ...
    ), seg_area, 1))';
    % If there is no matching keypoint, the index will be stored
    for p2 = 1:dim_mpt2(2)
        p2_seg = paddedI2(...
            Mpt2(2, p2):(Mpt2(2, p2) + twoPdg),  ...            
            Mpt2(1, p2):(Mpt2(1, p2) + twoPdg) ...
        );
        % Calculate normalized cross correlation value
        ncc = NCC(p1_seg, double(reshape(p2_seg, seg_area, 1)), seg_area);
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
filteredMatches = [max_ncc(~isnan(max_ncc) == 1); Mpt1(:, ~isnan(max_ncc) == 1); matches(:, ~isnan(max_ncc) == 1)];
% Remove dublicate matches
filteredMatches = sortrows(filteredMatches', 'descend');
[v, idx] = unique(filteredMatches(:, 4:5), 'rows');
Korrespondenzen = filteredMatches(idx, 2:5)';


%% Plotting
% TODO: Fix what they ask for at the end of the problem. Something with
% additional setting (?)
% TODO: Correspondence plotting has to be implemented manually (?)
showMatchedFeatures(I1, I2, Korrespondenzen(1:2, :)', Korrespondenzen(3:4, :)');
if(false) %do_plot
    subplot(1, 2, 1);
    imshow(I1);
    hold on;
    plot(Mpt1(1, :), Mpt1(2, :), 'b.');
    plot(Korrespondenzen(1, :), Korrespondenzen(2, :), 'r.');
    subplot(1, 2, 2);
    imshow(I2);
    hold on;
    plot(Mpt2(1, :), Mpt2(2, :), 'b.');
    plot(Korrespondenzen(3, :), Korrespondenzen(4, :), 'r.');
end

end

