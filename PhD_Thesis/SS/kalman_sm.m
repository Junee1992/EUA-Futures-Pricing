function [save_xt_1n, save_Pt_1n] = kalman_sm(save_Ptt, save_Ptt_1, save_att, save_att_1, G, a0, P0)

nobsn = size(save_att_1,1);
% Kalman_smoothing
for i = nobsn: -1 : 2
    if i == nobsn
        J = save_Ptt(:,:,i-1) * G' * inv(save_Ptt_1(:, :, i));
        xt_1n = save_att(i-1, :)' + J * (save_att(i, :)' - save_att_1(i, :)');
        Pt_1n = save_Ptt(:, :, i-1) + J * (save_Ptt(:, :, i) - save_Ptt_1(:, :, i)) * J';
    elseif i > 1
        J = save_Ptt(:, :, i-1) * G' * inv(save_Ptt_1(:, :, i));
        xt_1n = save_att(i-1, :)' + J * (save_xt_1n(i,:)' - save_att_1(i, :)');
        Pt_1n = save_Ptt(:, :, i-1) + J * (save_Pt_1n(:, :, i) - save_Ptt_1(:, :, i)) * J';
    end
    save_xt_1n(i-1, :) = xt_1n';
    save_Pt_1n(:, :, i-1) = Pt_1n;
    save_xt_1n(end,:) = save_att(end,:);
    save_Pt_1n(:,:,end) = save_Ptt(:,:,end);
end

