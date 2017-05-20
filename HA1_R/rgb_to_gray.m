function [Gray_image] = rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es direkt zurückgegeben werden.
s = size(Image);
    if s(1,3) == 1
        Gray_image = Image;
        return;
    end
Gray_image = 0.299 * Image(:,:,1) + 0.587 * Image(:,:,2) + 0.114 * Image(:,:,3);
end
