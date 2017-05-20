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
checkerImage = imread('checkerboard.png'); %Used for testing
IGray = rgb_to_gray(checkerImage);

[Fx, Fy] = sobel_xy(IGray);

%% Harris-Merkmale berechnen
%  tic;
k = 0.05;
tau = 10000;
segmenth_length = 3;
min_distance = 10; %Minimum distance between two features
tile_size = 11; %Divides the image in quadratic tiles, this variable describes the width/height
N = 3; %Maximal number of features within each tile
[Hvalues, Merkmale] = harris_detektor(IGray, Fx, Fy, segmenth_length, k, tau, 'do_plot', min_distance, tile_size, N);
%  toc;
