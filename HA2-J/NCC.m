function [ norm_corr ] = NCC( V1_s, V2_s, N)
% Calculates the correlation between two image segments
% TODO: Check if normalization is correctly implemented
% TODO: Optimize for speed? (Matrix multiplication)
% see: http://paulbourke.net/miscellaneous/correlate/
    %corr = V1_s' * V2_s;
    %norm1 = V1_s' * V1_s;
    %norm2 = V2_s' * V2_s;
    %norm_corr = corr / (sqrt(double(norm1) * double(norm2)));
    corr = V1_s / (sqrt((1 / (N - 1)) * (norm(V1_s, 'fro')^2))) * ...
           V2_s / (sqrt((1/(N - 1)) * (norm(V2_s, 'fro')^2)));

    norm_corr = corr/(N-1);
end

