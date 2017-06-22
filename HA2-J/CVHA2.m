%  Gruppennummer:
%  Gruppenmitglieder:

%% Hausaufgabe 2
%  Bestimmung von Punktkorrespondenzen zwischen Merkmalspunkten einer Stereo
%  Aufnahme.

%  F�r die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter �ber den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden k�nnen.


%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);

%% Korrespondenzsch�tzung
Korrespondenzen = punkt_korrespondenzen(IGray1, IGray2, Merkmale1, Merkmale2, 'do_plot', true);

%% Comparison using inbuilt matlab functions
figure;
points1 = detectHarrisFeatures(IGray1);
points2 = detectHarrisFeatures(IGray2);
[f1, vpts1] = extractFeatures(IGray1, points1);
[f2, vpts2] = extractFeatures(IGray2, points2);
indexPairs = matchFeatures(f1, f2) ;
matchedPoints1 = vpts1(indexPairs(1:20, 1));
matchedPoints2 = vpts2(indexPairs(1:20, 2));
figure; ax = axes;
showMatchedFeatures(IGray1,IGray2,matchedPoints1,matchedPoints2,'Parent',ax);
title(ax, 'Candidate point matches with inbuilt matlab functions');
legend(ax, 'Matched points 1','Matched points 2');