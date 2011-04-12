%**************************************************************************
% This function computes the C-alpha statistic in Neale et al.,
% "Testing for an Unusual Distribution of Rare Variants". Since permutation
% is applied to assess the p value, normalization of the statistic is not
% required.
%**************************************************************************

function [score, w] = getScoreCAlpha(ctrls, cases, mode)
    % Check whether there are missing values in the genotype data.
    casesNonMissPtr = (cases >= 0);
    ctrlsNonMissPtr = (ctrls >= 0);
    if (any(ctrls(:) < 0) || any (cases(:) < 0))
        hasMiss = true;
    end
    
    % Modify genotype values based on the mode of inheritance.
    if exist('mode', 'var')
        switch lower(mode)
            case 'dominant'
                ctrls = (ctrls >= 1);
                cases = (cases >= 1);
            case 'recessive'
                ctrls = (ctrls == 2);
                cases = (cases == 2);
        end
    end
    
    % For each variant, count the number of minor alleles in controls (ku)
    % and the number of minor alleles in cases (ka) respectively.
    ku = sum(ctrls .* ctrlsNonMissPtr, 1);
    ka = sum(cases .* casesNonMissPtr, 1);
    n = (ku + ka);
    
    % Compute p_0 (the probability of observing a variant in the cases) for
    % each SNP.
    n_cases = sum(casesNonMissPtr, 1);
    n_ctrls = sum(ctrlsNonMissPtr, 1);
    p_0 = n_cases ./ (n_cases + n_ctrls);
    
    % Exclude SNPs without minor alleles.
    ptr = (n > 0);
    ku = ku(ptr); ka = ka(ptr); n = n(ptr);
    n_cases = n_cases(ptr); n_ctrls = n_ctrls(ptr); p_0 = p_0(ptr);
    
    % Group singletons.
    ptr = (n == 1);
    if (sum(ptr) > 0)
        ka1 = sum(ka(ptr)); n1 = sum(n(ptr));
        p_01 = sum(n_cases(ptr)) / (sum(n_cases(ptr)) + sum(n_ctrls(ptr)));
        w1 = (ka1 - n1 * p_01) .^ 2 - n1 * p_01 * (1 - p_01);
    else
        w1 = 0;
    end
    
    % Non-singletons.
    ptr = (n > 1);
    ka2 = ka(ptr); n2 = n(ptr); p_02 = p_0(ptr);
    w2 = (ka2 - n2 .* p_02) .^ 2 - n2 .* p_02 .* (1 - p_02);
    
    score = w1 + sum(w2);
    w = [w1, w2];
end

