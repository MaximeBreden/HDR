function uv = convo_cos(u,v)

% Computing the discrete convolution corresponding to cosine series (could
% be made faster by using the FFT instead of building the multiplication
% matrix).

uv = convo_cos_mat(u) * v;