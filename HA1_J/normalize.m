function [normA] = normalize(A)

minimum = min(A(:));
maximum = max(A(:));

range = maximum - minimum;
add = abs(minimum);

normA = A;

for y = 1:size(A,1)
    for x = 1:size(A,2)
        normA(y,x) = ceil((A(y,x) + add)/range * 255);
    end
end