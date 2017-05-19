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
IGray = rgb_to_gray(Image);
[sobel_x, sobel_y] = sobel_xy(IGray);
%sobel_mag = gradientMagnitude(sobel_x, sobel_y);
%sobel_mag = normalize(sobel_mag);
%imshow(uint8(sobel_mag))

%% Harris-Merkmale berechnen
%  tic;
Merkmale = harris_detektor(IGray, Fx, Fy);
imshow(Merkmale);
%  toc;
