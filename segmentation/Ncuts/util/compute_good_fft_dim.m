function p_best = compute_good_fft_dim(p);
% p_best=2^ceil(log2(p));
% return;
% Timothee Cour, 29-Aug-2006 07:49:15

pmax = 2^ceil(log2(p));
p_best = p;
for k=p:pmax
    z = factor(k);
    if (all(z<=5) && sum(z~=2)<=2)
%     if (all(z<=3) && sum(z~=2)<=2)
        p_best = k;
        return;
    end
end

