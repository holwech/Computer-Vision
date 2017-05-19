function [Gray_image] = rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es direkt zurückgegeben werden.

sz = size(Image);
if size(sz,2) == 3
    val = 0.299 * Image(:,:,1) + 0.587 * Image(:,:,2) + 0.114 * Image(:,:,3);
    Gray_image = ceil(val);
else
    Gray_image = Image;

end
