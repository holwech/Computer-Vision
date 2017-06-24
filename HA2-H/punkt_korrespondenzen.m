function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
tic
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
I1_mean = mean2(I1); I2_mean = mean2(I2);
% whos I1_mean I2_mean
matches = zeros(2, size(Mpt1,2)); % for each interest point in Mpt1 store interest point of Mpt2 into matches
max_ncc = NaN(1, size(Mpt1,2)); % for each interest point in Mpt1 that has a match store the amount of correlation into max_ncc

%% Image padding
% Padding is added to the images, in case of comparisons of points close to the
% the edge.
pdg = floor(window_length / 2);
twoPdg = 2 * pdg; % Speed optimization
paddedI1 = padarray(I1, [pdg pdg], 'symmetric');
paddedI2 = padarray(I2, [pdg pdg], 'symmetric');

%% Match search
% Search through each keypoint and find the point with max correlation
for p1 = 1:size(Mpt1,2)
    p1_seg = paddedI1(...
        Mpt1(2, p1):(Mpt1(2, p1) + twoPdg),  ...
        Mpt1(1, p1):(Mpt1(1, p1) + twoPdg) ...
    );
    % If there is no matching keypoint, the index will be stored
    for p2 = 1:size(Mpt2,2)
        p2_seg = paddedI2(...
            Mpt2(2, p2):(Mpt2(2, p2) + twoPdg),  ...            
            Mpt2(1, p2):(Mpt2(1, p2) + twoPdg) ...
        );
        % Calculate normalized cross correlation value
        ncc = NCC(p1_seg, p2_seg, I1_mean, I2_mean);
        % If NCC is higher than previous matches and over threshold, store
        % the values in matches and max_ncc
        if (( ncc > max_ncc(p1) || isnan(max_ncc(p1)) ) && (ncc > min_corr))
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
toc

%% Plotting
if(do_plot)
    %Plot feature points as blue dots and feature points that have a match
    %in red.
    figure();
    sx = size(I1,2);
    x1 = Korrespondenzen(1,:); y1 = Korrespondenzen(2,:);
    x2 = Korrespondenzen(3,:); y2 = Korrespondenzen(4,:);
    x2 = x2 + sx;
    
    %find the right classified ones
    threshold = 300;
    bool = abs(y1-y2)<threshold & abs(x1-x2+sx)<threshold;
    x1_good = x1(bool); y1_good = y1(bool);
    x2_good = x2(bool); y2_good = y2(bool);
    x1_bad = x1(~bool); y1_bad = y1(~bool);
    x2_bad = x2(~bool); y2_bad = y2(~bool);
    
    %Plot labels to the dots
    imshow([I1 I2]);
    hold on;
    plot([Mpt1(1, :),Mpt2(1, :)+ sx] , [Mpt1(2, :),Mpt2(2, :)], 'b.', 'DisplayName','Alle Merkmalspunkte');
    plot([x1_bad,x2_bad], [y1_bad,y2_bad], 'r.', 'DisplayName', 'Merkmalspunkte mit Korrespondenzen');
    plot([x1_good,x2_good], [y1_good,y2_good], 'g.', 'DisplayName', 'Merkmalspunkte mit vermutlich richtigen Korrespondenzen');
    hold off;
    label_numbers = 1:size(Korrespondenzen,2);
    label_bad = cellstr(num2str(label_numbers(~bool)'));
    text(x1_bad, y1_bad, label_bad, 'Color','red','FontSize', 8);
	text(x2_bad, y2_bad, label_bad, 'Color','red','FontSize',8);
    label_good = cellstr(num2str(label_numbers(bool)'));
    text(x1_good, y1_good, label_good, 'Color','green','FontSize', 8);
	text(x2_good, y2_good, label_good, 'Color','green','FontSize',8);
    legend('Location', 'northoutside');
    
    %Plot lines between the matching points
    figure();
    imshow([I1 I2]);
    hold on;
    plot([x1_bad' x2_bad']', [y1_bad' y2_bad']', 'r', 'LineWidth', 1);
    plot([x1_good' x2_good']', [y1_good' y2_good']', 'g', 'LineWidth', 1);
    hold off;
    
    
    fprintf('Classification Quote: %.2f \n',sum(bool)/length(bool));
end
end

