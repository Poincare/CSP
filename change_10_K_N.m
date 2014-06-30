function [vect]=change_10_K_N(dec,K,N)
% dec<=N^K
% map [1: N^K] to [0:N-1]^K *******start from 0*********
dec=dec-1;
vect = zeros(1,K);
for k=1:K
    vect(K-k+1)=mod(dec,N);
    dec=floor(dec./(N));
    if dec==0
        break
    end
end
