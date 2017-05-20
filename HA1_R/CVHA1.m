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
figure();
subplot(2,2,1);
imshow(Image);
IGray = rgb_to_gray(Image);
subplot(2,2,2);
imshow(IGray);
[Fx,Fy]=sobel_xy(IGray);
subplot(2,2,3);
imshow(Fx);
subplot(2,2,4);
imshow(Fy);

%% Harris-Merkmale berechnen
tic;
Merkmale = harris_detektor(IGray,'do_plot',true);
toc;
