function subset_indices = GetSubsetIndices(subsetNum)
% Usage: subset_indices = GetSubsetIndices(subsetNum)
% 
% This function returns the indices of a specified subset
% NOTE: This assumes that there are 21 subsets of 252 projection angles

nAnglesPerSubset = 12;
nSubsets = 21;

if subsetNum > 0
    % This is to ensure that each run of this function returns the same random
    % permutation of index numbers:
    s = RandStream('mt19937ar','Seed',0);

    % - generate random permutation of index numbers
    rand_idx_list = randperm(s,252,252);

    % - determine which indices to return for this subset
    idx_start = 1 + nAnglesPerSubset*(subsetNum-1);
    idx_end = nAnglesPerSubset*subsetNum;

    subset_indices = rand_idx_list(idx_start:idx_end);
else
    subset_indices = 1:(nAnglesPerSubset*nSubsets);
end

    

end