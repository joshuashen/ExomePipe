%**************************************************************************
% This function returns the one-sided p value and its corresponding test
% score by performing permutations.
% Note that if specified, freq should be the frequencies of alleles
% (additive model) or genotypes (01 for dominant model; 11 for recessive).
% All elements in freq must be greater than zero.
%**************************************************************************

function [p, score, w, p_twoSided] = getPermPval(ctrls, cases, num_perm, ...
        method, mode, twoSided, freq)
    num_ctrls = size(ctrls, 1);
    num_cases = size(cases, 1);
    permCtrls = ctrls;
    permCases = cases;
    combined = [ctrls; cases];
    
    if ~exist('freq', 'var')
        freq = [];
    end
    
    for i = 1: (num_perm + 1)
        switch lower(method)
            case 'sigmamidp'
                [s(i), w1] = getScoreSigmaP(permCtrls, permCases, ...
                        mode, true, twoSided, freq);
            case 'sigmap'
                [s(i), w1] = getScoreSigmaP(permCtrls, permCases, ...
                        mode, false, twoSided, freq);
            case 'calpha'
                [s(i), w1] = getScoreCAlpha(permCtrls, permCases, mode);
            case 'vt'
                [s(i), w1] = getScoreVT(permCtrls, permCases, ...
                        mode, twoSided);
            case 'poisson'
                [s(i), w1] = getScorePoisson(permCtrls, permCases, ...
                        mode, twoSided);
            case 'weightedsum'
                [s(i), w1] = getScoreWeightedSum(permCtrls, permCases, ...
                        mode, false, false, freq);
            case 'weightedsumscore'
                [s(i), w1] = getScoreWeightedSum(permCtrls, permCases, ...
                        mode, true, false, freq);
            case 'modifiedweightedsumscore'
                [s(i), w1] = getScoreWeightedSum(permCtrls, permCases, ...
                        mode, true, true, freq);
        end
        
        if (i == 1)
            w = w1(:)';
        end
        
        % Permute the controls and cases.
        idx = randperm(size(combined, 1));
        permCtrls = combined(idx(1: num_ctrls), :);
        permCases = combined(idx(num_ctrls + 1: size(combined, 1)), :);
    end
    
    score = s(1);
    s = s(2: num_perm + 1);
    p = sum(s >= score) / num_perm;

    if (strcmpi(method, 'weightedsum') || strcmpi(method, 'weightedsumscore'))
        p_twoSided = sum(s >= score) / num_perm;
        if (p_twoSided > 0.5)
         p_twoSided = 1 - p_twoSided;
        end
        p_twoSided = p_twoSided * 2;
    elseif strcmpi(method, 'calpha')
        % The method may have negative scores, but risk and protectiveness
        % are not distinguished and both identified by a statistic larger
        % than expected using a one-tailed test.
        p_twoSided = p;
    else
        p_twoSided = sum(abs(s) >= abs(score)) / num_perm;
    end
    
end

