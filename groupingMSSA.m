function [ restored ] = groupingMSSA(ppgSample1, ppgSample2) 
%workDraftGrouping Decompose a signal, find harmonics pairs, then restore.

% General SSA parameters.
% 20 to allow for more factor vects needed on noisier signals.

pcNumMax = 20; % robust for 12-20, 14/16 is good middle ground.
M = 1000;
L = 334;% for VMSSA instead of M/2 -> from VMSSA theory.
K = M-L+1;
N=4096;

% SSA Decomposition into U,S,V eigentriples.

% % Replacing statement above for trajectory matrix.
mi = repmat(1:K,[L,1]) + repmat((1) * (0:(L-1))',[1,K]);
trajMat1 = ppgSample1(mi);
trajMat2 = ppgSample2(mi);
trajMat = [trajMat1;trajMat2];

% Normalizing data, necessary for MSSA.
trajMat = zscore(trajMat);

[L1, K1] = size(trajMat);

% Matrix used in eigen-calculation to determine U,S, and V.
% Use covariance instead, more accurate.
[U, eigV] = eigs(cov(trajMat'), pcNumMax);
S = sqrt(eigV);

% Finding how many pc's to grab to get ~90% of signal energy. 
% Server to prevent emphasis on noise.
pcNumIdeal = 0;
eigSum = diag(eigV);
eigSum = log10(eigSum);
eigSum = sum(eigSum);
eigNum = 0;
currentEigV = 1;

while (eigNum/eigSum <= 0.75)
    pcNumIdeal = pcNumIdeal + 1;
    eigNum = eigNum + log10(eigV(currentEigV, currentEigV));
end
pairNumIdeal = pcNumIdeal-1;

% Building the factor vectors from left singular vectors and trajMat. 
V = zeros(K1,pcNumIdeal);
for i = 1:pcNumIdeal
    V(:,i) = trajMat'*( U(:,i) );
end

% fftPC will hold the N-point fft of the principal components.
% Use freq spectrum used to find pairs by max bin value in freq dom.

fftPC = zeros(N,pcNumIdeal); %!!!!!! testing PCNumIdeal

for i = 1:pcNumIdeal
    fftPC(:,i) = fft(V(:,i),N); % !!! testingVPC fft(allPCs(:,1,i), N);
    fftPC(:,i) = abs(fftPC(:, i));
    fftPC(:,i) = fftPC(:,i)/max(fftPC(:,i));
end

% Find frequency pairs no more than two samples away.
% Finding which pairs to pass through.
harmonicsInd = zeros(1, pairNumIdeal);
pairsToSum = zeros(1,pairNumIdeal);

% Pairs should only be in possible heartbeat range.
thr = 2; 
for i = 1:pairNumIdeal
    % Initializing max indeces, and PC's to compare.
    maxIndex1 = 0;
    maxIndex2 = 0;
    tol = 0.1; % Tolerance for finding maximum in frequency domain.
    a=fftPC(:,i);
    b=fftPC(:,i+1);
    for j = 33:236
        if abs( a(i) - 1 ) < tol
            maxIndex1 = j;
        end
        if abs( b(i) - 1) < tol
            maxIndex2 = j;
        end
    end
    if abs(maxIndex1-maxIndex2) <= thr
        % Store the index of frequency for oscillatory pairs.
        harmonicsInd(i) = round((maxIndex1+maxIndex2)/2);
        pairsToSum(i) = 1;
    end
end

% Will hold sum of the chosen components.
groupedPCs = zeros(L1, K1);

for i = 1:pairNumIdeal
    if (pairsToSum(i)) == 1
        groupedPCs = groupedPCs + U(:,i)*S(i,i)*V(:,i)' + U(:,i+1)*S(i+1,i+1)*V(:,i+1)'; % !!! testingVPC allPCs(:,:,i) + allPCs(:,:,i+1);
        i=i+1; % to prevent double counting a factVect.
    end
end

% Diagonal averaging.
groupedPCs = fliplr(groupedPCs);
restored = arrayfun(@(i) mean(diag(groupedPCs,i)), size(groupedPCs,2)-1:-1:-size(groupedPCs,1)+1).';
restored = restored(1:M);

% % Temporal Difference
temp = zeros(1,M);
restored = diff(restored,2);
for i = 1:length(restored)
    temp(i) = restored(i);
end

restored = temp;
restored = restored - mean(restored);
restored = restored/std(restored);

end