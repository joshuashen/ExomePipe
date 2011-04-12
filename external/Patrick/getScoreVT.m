%**************************************************************************
% This function implements the variable threshold (VT) statistics in Price
% et al., "Pooled Association Tests for Rare Variants in Exon-Resequencing
% Studies".
%**************************************************************************

function [score, w] = getScoreVT(ctrls, cases, mode, twoSided)
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
    
    % For each variant, count the number mutant alleles in controls (mu).
    % The weight w is the estimated standard deviation of the total number
    % of mutations in the sample under the null hypothesis of no frequency
    % differences between affected and unaffected individuals.
    switch lower(mode)
        case ('recessive')
            % 0 for genotypes 00 or 01, 1 for genotype 11
            Iu = (ctrls == 2); Ia = (cases == 2);
        case ('dominant')
            % 0 for genotype 00, 1 for genotypes 01 or 11
            Iu = (ctrls >= 1); Ia = (cases >= 1);
        case {'additive', 'casemaf'}
            % 0, 1, 2 for genotypes 00, 01 and 11 respectively.
            Iu = ctrls .* ctrlsNonMissPtr;
            Ia = cases .* casesNonMissPtr;
    end
    
    Ca = sum(Ia, 1); Cu = sum(Iu, 1); Csq = sum([Iu; Ia] .^ 2, 1);
    meanPhen = (size(cases, 1) - size(ctrls, 1)) / (size(cases, 1) + size(ctrls, 1));
    
    alleleCount = sum([ctrls; cases] .* [ctrlsNonMissPtr; casesNonMissPtr], 1);
    if (hasMiss)
        alleleFreq = alleleCount ./ ...
            (sum([ctrlsNonMissPtr; casesNonMissPtr], 1) .* 2);
        s_af = sort(unique(alleleFreq));
        for i = 1: length(s_af)
            ptr = (alleleFreq == s_af(i));
            phenCount(i) = sum(Ca(ptr) - Cu(ptr));
            countMeanPhen(i) = sum(Ca(ptr) + Cu(ptr)) * meanPhen;
            countSq(i) = sum(Csq(ptr));
        end
    else
        s_ac = sort(unique(alleleCount));
        for i = 1: length(s_ac)
            ptr = (alleleCount == s_ac(i));
            phenCount(i) = sum(Ca(ptr) - Cu(ptr));
            countMeanPhen(i) = sum(Ca(ptr) + Cu(ptr)) * meanPhen;
            countSq(i) = sum(Csq(ptr));
        end
    end
    
    csPhenCount = cumsum(phenCount);
    csCountMeanPhen = cumsum(countMeanPhen);
    sqrtCsCountSq = sqrt(cumsum(countSq));
    
    if (sqrtCsCountSq(1) == 0)
        csPhenCount(1) = [];
        csCountMeanPhen(1) = [];
        sqrtCsCountSq(1) = [];
        if (hasMiss)
            s_af(1) = [];
        else
            s_ac(1) = [];
        end
    end
    
    z = (csPhenCount - csCountMeanPhen) ./ sqrtCsCountSq;
    
    if (twoSided)
        [score, idx] = max(abs(z));
        score = z(idx(1));
    else
        [score, idx] = max(z);
    end
    
    if (hasMiss)
        w = sum(alleleFreq <= s_af(idx(1)));
    else
        w = sum(alleleCount <= s_ac(idx(1)));
    end
end

