function plot_korrespondenzen( I1,I2,Korrespondenzen )
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

