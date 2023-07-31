function mask = rateProfile(N, K, dsnrdB, snrType, construction)
    if strcmpi(construction, 'dega')
        mask = DEGABuild(N, K, dsnrdB, snrType);
    elseif strcmpi(construction, 'rmpolar')
        mask = RMPolarBuild(N, K, dsnrdB, snrType);
    elseif strcmpi(construction, 'rm')
        mask = RM(N, K);
    elseif strcmpi(construction, 'test')
        mask = test(N, K);
    else
        disp('None such type')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%Density Evolution With Gaussian Approximation Construction%%%%%%%%%%%%%%%%%%%%%%%%%

function mask = DEGABuild(N, K, dsnrdB, snrType)
    mu = GA(N, K, dsnrdB, snrType);
    [~, I] = maxk(mu, K);
    mask = zeros(1, N);
    mask(I) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%Reed-Muller Polar Construction%%%%%%%%%%%%%%%%%%%%%%%%%

function mask = RMPolarBuild(N, K, dsnrdB, snrType)
    mask = ones(1, N);
    weights = GRowsWeight(N);
    mu = GA(N, K, dsnrdB, snrType);
    weightCount = zeros(1, log2(N) + 1);
    for weight = weights
        weightCount(weight + 1) = weightCount(weight + 1) + 1;
    end

    bitCount = 0;
    k = 1;
    while (bitCount + weightCount(k)) <= (N - K)
        for i = 1:N
            if weights(i) == (k - 1)
                mask(i) = 0;
                bitCount = bitCount + 1;
            end
        end
        k = k + 1;
    end

    mask2 = [];
    for i = 1:N
        if weights(i) == (k - 1)
            mask2(end + 1) = i;
        end
    end
    [~, I] = sort(mu(mask2));
    mask2 = mask2(I);
    remainder = (N - K) - bitCount;
    for i = 1:remainder
        mask(mask2(i)) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%Reed-Muller Construction%%%%%%%%%%%%%%%%%%%%%%%%%

function mask = RM(N, K)
    weights = sum(dec2bin(0:N-1)-'0',2);
    [~, I] = sort(weights,'descend');
    mask = zeros(1, N);
    mask(I(1:K)) = 1;
end

function mask = test(N, K)
    w = sum(dec2bin(0:N-1)-'0',2);
    [~, index] = sort(w,'descend');
    kindex = index(1:K);
    mask = zeros(1,N);
    mask(kindex) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%Auxiliary Functions%%%%%%%%%%%%%%%%%%%%%%%%%

function weights = GRowsWeight(N)
    weights = zeros(1, N);
    n = log2(N);
    for i = 1:N
        ones = 0;
        binExp = dec2bin(i - 1, n);
        for num = binExp
            if num == '1'
                ones = ones + 1;
            end
        end
        weights(i) = ones;
    end
end

function mu = GA(N, K, dsnrdB, snrType)
    if strcmpi(snrType, 'snr')
        sigmaSquare = 1 / (2 * power(10, dsnrdB / 10));
    elseif strcmpi(snrType, 'snrb')
        R = K / N;
        sigmaSquare = 1 / (2 * R * power(10, dsnrdB / 10));
    end

    mu = 2 / sigmaSquare;
    while length(mu) ~= N
        tmp = [];
        for i = 1:length(mu)
            tmp(end + 1) = phiInv(mu(i));
            tmp(end + 1) = 2 * mu(i);
        end
        mu = tmp;
    end
end

function y = phiInv(x)
    if x > 12
        y = 0.9861 * x - 2.3152;
    elseif x <= 12 && x > 3.5
        y = x * (9.005e-3 * x + 0.7694) - 0.9507;
    elseif x <= 3.5 && x > 1
        y = x * (6.2883e-2 * x + 0.3678) - 0.1627;
    else
        y = x * (2.202e-1 * x + 0.06448);
    end
end