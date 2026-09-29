function features = extract_dFC_features(wpli_mats, nClusters)
% EXTRACT_DFC_FEATURES extracts four dynamic functional connectivity features from the dynamic wPLI matrix.
%
% Input:
%   wpli_mats - A channels×channels×nWindows wPLI matrix, with the third dimension representing the time window.
%   nClusters - Optional. 
%
% Output：
%   features  - 1×4 Vector：
%       [GlobalVariability, Metastability, ComplexityMean, ComplexityStd]
%

    if nargin < 2
        nClusters = []; 
    end

    nChannels = size(wpli_mats, 1);
    nWindows  = size(wpli_mats, 3);

    global_strength = zeros(nWindows, 1);
    upperVecs = zeros(nWindows, nChannels*(nChannels-1)/2);

    for w = 1:nWindows
        mat = wpli_mats(:, :, w);
        upperVecs(w, :) = mat(triu(true(nChannels), 1));
        global_strength(w) = sum(upperVecs(w, :));
    end

    % Global variability: The average Euclidean distance of the triangular connection vectors over each time window
    distMat = squareform(pdist(upperVecs, 'euclidean'));
    upperTriDist = distMat(triu(true(nWindows), 1));
    if ~isempty(upperTriDist)
        GlobalVariability = mean(upperTriDist);
    else
        GlobalVariability = NaN;
    end

    % Metastable state: Standard deviation of the global connection strength
    Metastability = std(global_strength);

    % Von Neumann entropy
    entropy_vn = zeros(nWindows, 1);
    for w = 1:nWindows
        lambda = eig(wpli_mats(:, :, w));
        lambda = lambda(lambda > 1e-10);
        lambda = lambda / sum(lambda);
        if ~isempty(lambda)
            entropy_vn(w) = -sum(lambda .* log(lambda));
        end
    end
    ComplexityMean = mean(entropy_vn);
    ComplexityStd  = std(entropy_vn);

    features = [GlobalVariability, Metastability, ComplexityMean, ComplexityStd];
end
