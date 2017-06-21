%  Gruppennummer:
%  Gruppenmitglieder:

%% Hausaufgabe 1
%  Einlesen und Konvertieren von Bildern sowie Bestimmung von 
%  Merkmalen mittels Harris-Detektor. 

%  F�r die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter �ber den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden k�nnen.


%% Bild laden
Image = imread('szene.jpg');
%Image = imread('checkerboard.png'); %Used for testing

IGray = rgb_to_gray(Image);
[Fx, Fy] = sobel_xy(IGray);

% Plotten
figure;
subplot(1,3,1);
imshow(Image);
title('Original Image')
subplot(1,3,2);
imshow(Fx);
title('Sobel Fx')
subplot(1,3,3);
imshow(Fy)
title('Sobel Fy');

%% Harris-Merkmale berechnen
tic;
%min_distance   : Minimum distance between two features
%tile_size      : Divides the image in quadratic tiles, this variable describes the width/height
%N              : Maximal number of features within each tile
[Merkmale] = harris_detektor(IGray, 'segment_length', 3, 'k', 0.05, 'tau', 1, 'do_plot', true, 'min_distance', 10, 'tile_size', 11, 'N',  3);
toc;
