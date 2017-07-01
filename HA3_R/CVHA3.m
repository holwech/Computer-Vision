%  Gruppennummer:
%  Gruppenmitglieder:

%% Hausaufgabe 3
%  Bestimmung von robusten Korrespondenzpunktpaaren mittels RANSAC und
%  Berechnung der Essentiellen Matrix.

%  Für die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter über den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden können.


%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);



%% Korrespondenzschätzung
tic;
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,'min_corr',0.92,'do_plot',true);
zeit_korrespondenzen = toc;
disp(['Es wurden ' num2str(size(Korrespondenzen,2)) ' Korrespondenzpunktpaare in ' num2str(zeit_korrespondenzen) 's gefunden.'])

%% Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus 
Korrespondenzen_robust = F_ransac(Korrespondenzen)
if ~isempty(Korrespondenzen_robust)
    plot_korrespondenzen(IGray1,IGray2,Korrespondenzen_robust);
end
%% Berechnung der Essentiellen Matrix
% load('K.mat');
% E = achtpunktalgorithmus(Korrespondenzen_robust,K);

