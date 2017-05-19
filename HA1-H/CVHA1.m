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
%tImage = imread('engine.png');
tImage = imread('test.jpg');
checkerImage = imread('szene.jpg');
IGray = rgb_to_gray(checkerImage);
%RIGraySegment = IGray(800:1300,800:1300);

[Fx, Fy] = sobel_xy(IGray);
mag = gradientMagnitude(Fx,Fy);
 nMag = normalize(mag);
 imshow(uint8(nMag));

%% Harris-Merkmale berechnen
%  tic;
%  Merkmale = harris_detektor(IGray,'do_plot',true);
%SFx = [1 0 -1; 2 0 -2; 1 0 -1];
%SFy = [1 2 1; 0 0 0; -1 -2 -1];
%Fx = conv2(IGray, SFx, 'same');
%Fy = conv2(IGray, SFy, 'same');
%[ev, H, Merkmale] = harris_detektor(IGray, Fx, Fy, 3, 0.05, 100000000000);
%  toc;
