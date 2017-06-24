%  Gruppennummer: M13
%  Gruppenmitglieder: Joachim Holwech, Regine Hartwig, Bjørn Hoffmann

%% Hausaufgabe 2
%  Bestimmung von Punktkorrespondenzen zwischen Merkmalspunkten einer Stereo
%  Aufnahme.

%  Fuer die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter ueber den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden koennen.

clc
%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);

%% Korrespondenzschaetzung
Korrespondenzen = punkt_korrespondenzen(IGray1, IGray2, Merkmale1, Merkmale2, 'do_plot', true, 'window_length', 15, 'min_corr', 0.995);