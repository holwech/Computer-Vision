function [res] = gradientMagnitude(Fx, Fy)

sz = size(Fx);
res = zeros(sz);
for y = 1:size(Fx,1)
    for x = 1:size(Fx,2)
        res(y,x) = sqrt(Fx(y,x)*Fx(y,x) + Fy(y,x)*Fy(y,x));
    end
end