function [gamma_best, p_best, MSE_best] = GridMLSSVR(trnX, trnY, fold)
%
% [gamma_best, lambda_best, p_best, MSE_best] = GridMLSSVR(trnX, trnY, fold);
%
% author: XU, Shuo (pzczxs@gmail.com)
% date: 2011-12-25
%
gamma = 2.^(-5: 0.5: 10);
% lambda = 2.^(-10: 2: 10);
p = 2.^(-5: 0.5: 10);

m = size(trnY, 2);

% random permutation
% [trnX, trnY] = random_perm(trnX, trnY);

MSE_best = inf;
MSE = zeros(fold, m);
curR2 = zeros(1, m);
R2 = zeros(1, m);
for i = 1: length(gamma)
    % for j = 1: length(lambda)
    for k = 1: length(p)
        

        for v = 1: fold
            predictY = [];
            [train_inst, train_lbl, test_inst, test_lbl] = folding(trnX, trnY, fold, v);

            % [alpha, b] = MLSSVRTrain(train_inst, train_lbl, gamma(i), lambda(j), p(k));
            [alpha, b] = MLSSVRTrain(train_inst, train_lbl, gamma(i), p(k));

            % [tmpY, MSE(v, :)] = MLSSVRPredict(test_inst, test_lbl, train_inst, alpha, b, lambda(j), p(k));
            [tmpY, MSE(v, :)] = MLSSVRPredict(test_inst, test_lbl, train_inst, alpha, b, p(k));

            predictY = [predictY; tmpY];
        end

        curMSE = sum(sum(MSE)) / numel(trnY);

        if MSE_best > curMSE
            gamma_best = gamma(i);
            % lambda_best = lambda(j);
            p_best = p(k);
            MSE_best = curMSE;
        end

        % fprintf('gamma = %g, lambda = %g, p = %g, mean_MSE = %g (%g, %g, %g, %g)\n', ...
            % log2(gamma(i)), log2(lambda(j)), log2(p(k)), sum(sum(MSE))/numel(trnY), ...
            % log2(gamma_best), log2(lambda_best), log2(p_best), MSE_best);
        fprintf('gamma = %g, p = %g, mean_MSE = %g (%g, %g, %g)\n', ...
            log2(gamma(i)), log2(p(k)), sum(sum(MSE))/numel(trnY), ...
            log2(gamma_best), log2(p_best), MSE_best);
    end
% end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random permutation by swapping i and j instance for each class
function [svm_inst, svm_lbl] = random_perm(svm_inst, svm_lbl)
n = size(svm_inst, 1);
rand('state', 0);
for i = 1: n
    k = round(i + (n - i)*rand());   % [i, n]
    svm_inst([k, i], :) = svm_inst([i, k], :);
    svm_lbl([k, i], :) = svm_lbl([i, k], :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [train_inst, train_lbl, test_inst, test_lbl] = folding(svm_inst, svm_lbl, fold, k)
    % n is the total number of instances
    n = size(svm_inst, 1);
    
    % Compute the size of each fold
    fold_size = floor(n / (fold + 1));
    
    % Define the indices for the training and testing sets
    start_train = 1;
    end_train = start_train + k * fold_size - 1;
    start_test = end_train + 1;
    end_test = min(start_test + round(0.2*fold_size), n); % Ensure we do not exceed the dataset
    
    % Extract training instances and labels
    train_inst = svm_inst(start_train:end_train, :);
    train_lbl = svm_lbl(start_train:end_train, :);
    
    % Extract testing instances and labels
    test_inst = svm_inst(start_test:end_test, :);
    test_lbl = svm_lbl(start_test:end_test, :);


    % function [train_inst, train_lbl, test_inst, test_lbl] = folding(svm_inst, svm_lbl, fold, k)
    %     n = size(svm_inst, 1);
    % 
    %     % folding instances
    %     start_index = round((k - 1)*n/fold) + 1;
    %     end_index = round(k*n/fold);
    %     test_index = start_index: end_index;
    % 
    %     % extract test instances and corresponding labels
    %     test_inst = svm_inst(test_index, :);
    %     test_lbl = svm_lbl(test_index, :);
    % 
    %     % extract train instances and corresponding labels
    %     train_inst = svm_inst;
    %     train_inst(test_index, :) = [];
    %     train_lbl = svm_lbl;
    %     train_lbl(test_index, :) = [];