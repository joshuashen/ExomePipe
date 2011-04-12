%**************************************************************************
% This function implements the sum of ranks statistics in Madsen et al.,
% "A Groupwise Association Test for Rare Mutations Using a Weighted Sum
% Statistic".
%**************************************************************************

function [score, w] = getScoreWeightedSum(ctrls, cases, mode, ...
        sumScore, modified, freq)
    
    % Check whether there are missing values in the genotype data.
    casesNonMissPtr = (cases >= 0);
    ctrlsNonMissPtr = (ctrls >= 0);
    if (any(ctrls(:) < 0) || any (cases(:) < 0))
        hasMiss = true;
    else
        hasMiss = false;
    end    
    
    % Assume additive mode of inheritance if not specified.
    if ~exist('mode', 'var') || isempty(mode)
        mode = 'additive';
    end
    
    % Modified weighted-sum method by Price et al.
    if ~exist('modified', 'var')
        modified = false;
    end
    
    % For each variant, count the number mutant alleles in controls (mu).
    % The weight w is the estimated standard deviation of the total number
    % of mutations in the sample under the null hypothesis of no frequency
    % differences between affected and unaffected individuals.
    switch lower(mode)
        case ('recessive')
            % 0 for genotypes 00 or 01, 1 for genotype 11
            Iu = (ctrls == 2); Ia = (cases == 2);
            if (hasMiss)
                nu = sum(ctrlsNonMissPtr, 1);
            else
                nu = size(ctrls, 1);
            end
        case ('dominant')
            % 0 for genotype 00, 1 for genotypes 01 or 11
            Iu = (ctrls >= 1); Ia = (cases >= 1);
            if (hasMiss)
                nu = sum(ctrlsNonMissPtr, 1);
            else
                nu = size(ctrls, 1);
            end
        case {'additive', 'casemaf'}
            % 0, 1, 2 for genotypes 00, 01 and 11 respectively.
            Iu = ctrls .* ctrlsNonMissPtr;
            Ia = cases .* casesNonMissPtr;
            if (hasMiss)
                nu = sum(ctrlsNonMissPtr, 1) * 2;
            else
                nu = size(ctrls, 1) * 2;
            end
    end
    
    % Use maf in all samples to compute weights.
    % mu = sum([Iu; Ia], 1);
    % q = (mu + 1) / (2*size([Iu; Ia], 1) + 2);
    % n = size(ctrls, 1) + size(cases, 1);
    % w = sqrt(q .* (1 - q));
    
    % Use maf in controls only to compute weights.
    mu = sum(Iu, 1);
    if (exist('freq', 'var') && ~isempty(freq))
        q = freq;
    else
        q = (mu + 1) ./ (nu + 2);
    end
    n = sum(ctrlsNonMissPtr, 1) + sum(casesNonMissPtr, 1);
    w = sqrt(n .* q .* (1 - q));
    w = 1 ./ w;
    
    % The score of each individual
    vu = sum(Iu .* w(ones(1, size(Iu, 1)), :), 2);
    va = sum(Ia .* w(ones(1, size(Ia, 1)), :), 2);
    
    if (exist('sumScore', 'var') && sumScore)
        if (modified)
            score = sum(va) - sum(vu);
        else
            score = sum(va);
        end
    else
        % Rank of all individuals.
        [v, sIdx] = sort([va; vu]);
        [v, idx, rIdx] = unique(v);
        r(sIdx) = rIdx;
        
        % If n individuals at rank k have the same score, then the next
        % individual at a higher score is at rank (k + n).
        r_o = r;
        for i = 1: max(rIdx)
            count = sum(r_o == i);
            if (count > 1)
                r(r_o > i) = r(r_o > i) + count - 1;
            end
        end
        
        ra = r(1: size(cases, 1));
        score = sum(ra);
    end
end

