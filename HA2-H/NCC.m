function [ norm_corr ] = NCC( p1_seg, p2_seg, V_mean, W_mean)

%Normalize image segments
N = size(p1_seg,1)^2;
V = double(p1_seg);
W = double(p2_seg);
Vn = (V - V_mean)/( sqrt((1/(N-1))*(norm(V - V_mean, 'fro')^2)) );
Wn = (W - W_mean)/( sqrt((1/(N-1))*(norm(W - W_mean, 'fro')^2)) );

%Vectorize the normalized image segments
Vn_vec = double(reshape(Vn, N, 1));
Wn_vec = double(reshape(Wn, N, 1));

%Calculate the NCC
corr = Wn_vec'*Vn_vec;
norm_corr = corr/(N-1);

end

