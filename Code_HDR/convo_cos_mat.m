function Mu = convo_cos_mat(u)

% Building the multiplication matrix for cosine series.

Hu = hankel(u);
Tu = toeplitz(u,u);
Tu(:,1) = 0;
Mu = Hu + Tu;
