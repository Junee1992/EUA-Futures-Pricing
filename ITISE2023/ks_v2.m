function wt_1 = ks_v2(n_t, k, C, G, att, att_1, Ptt, Ptt_1)

as(:,:,n_t) = att(:,end);

for p = 1:k
    J = Ptt(:,:,n_t-p) * G' * inv(Ptt_1(:,:,n_t-p+1));
    as(:,:,n_t-p) = att(:,n_t-p) + J * (as(:,:,n_t-p+1) - att_1(:,n_t-p+1));
    wt_1(:,p) = as(:,:,n_t-p+1) - C - G * as(:,:,n_t-p);
 end
