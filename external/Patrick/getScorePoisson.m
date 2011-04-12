%**************************************************************************
% This function computes the weighted sum statistic in Ionita-Laza et al.,
% "A new testing strategy to identify rare variants with either risk or
% protective effect on disease".
%**************************************************************************

function [score, w] = getScorePoisson(ctrls, cases, mode, twoSided)
    % Check whether there are missing values in the genotype data.
    casesNonMissPtr = (cases >= 0);
    ctrlsNonMissPtr = (ctrls >= 0);
    if (any(ctrls(:) < 0) || any (cases(:) < 0))
        hasMiss = true;
    else
        hasMiss = false;
    end
        
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
    
    n_s = 1;
    if (exist('twoSided', 'var') && twoSided)
        n_s = 2;
    end
    
    % For each variant, count the number of minor alleles in controls (ku)
    % and the number of minor alleles in cases (ka) respectively.
    ku = sum(ctrls .* ctrlsNonMissPtr, 1)';
    ka = sum(cases .* casesNonMissPtr, 1)';
    
    score = []; w = {};
    for j = 1: n_s;
        if (j == 1)
            % Consider only risk variants.
            ptr = (ka > ku);
            k1 = ku(ptr); k2 = ka(ptr);
        else
            % Consider only protective variants.
            ptr = (ku > ka);
            k1 = ka(ptr); k2 = ku(ptr);
        end
        
        if all(ptr == 0)
            score(j) = 0; w{j} = [];
            continue;
        end
        
        % Compute the weights.
        if (hasMiss)
            % If there are missing values, the expected number of variants
            % for each SNP is different.
            n_ctrls = sum(ctrlsNonMissPtr(:, ptr), 1)';
            n_cases = sum(casesNonMissPtr(:, ptr), 1)';
            q = n_ctrls ./ (n_ctrls + n_cases);
            if (j == 2)
                q = 1 - q;
            end
            f = (k1 + k2) .* q;
            
            % Identify all unique groups (ku, ka) and their corresponding sizes.
            % Then compute the weight of each group.
            [gp, ptr1, ptr2] = unique([k1, k2, f], 'rows');
            f1 = gp(:, 3);
            f2 = gp(:, 1) + gp(:, 2) - f1;            
        else
            [gp, ptr1, ptr2] = unique([k1, k2], 'rows');
            f1 = sum(gp, 2) / 2;
            f2 = f1;
        end
        
        n = [];
        for i = 1: size(gp, 1)
            n(i, 1) = sum(ptr2 == i);
        end
        
        if ~isempty(gp)
            w{j} = poisscdf(gp(:, 1), f1) .* (1 - poisscdf(gp(:, 2) - 1, f2));
            w{j} = -log(w{j});
            score(j) = sum(w{j} .* n);
        else
            w{j} = 0;
            score(j) = 0;
        end
    end
    
    [score, idx] = max(score);
    w = w{idx};
end

