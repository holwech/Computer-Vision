function [ norm_corr ] = NCC( p1_seg, p2_seg, I1_mean, I2_mean )
% Calculates the correlation between two image segments
% TODO: Check if normalization is correctly implemented
% TODO: Optimize for speed? (Matrix multiplication)
% see: http://paulbourke.net/miscellaneous/correlate/
    corr = sum(sum((int32(p1_seg) - I1_mean) .* (int32(p2_seg) - I2_mean)));
    norm1 = sum(sum((int32(p1_seg) - I1_mean) .^ 2));
    norm2 = sum(sum((int32(p2_seg) - I2_mean) .^ 2));
    norm_corr = corr / (sqrt(double(norm1)) * sqrt(double(norm2)));
end

