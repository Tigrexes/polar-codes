function [decoded, iterations] = decoder(y, llr, mask, crcPolynomial, convG, parameter, dsnrdB, isCRC, snrType, decoderType, msg)
    if strcmpi(decoderType, 'sphere')
        [decoded, iterations] = polarSphere(y, llr, mask, inf);
    elseif strcmpi(decoderType, 'casd')
        [decoded, iterations] = CA_SD(y, llr, mask, crcPolynomial, inf, zeros(1, length(y)));
    elseif strcmpi(decoderType, 'sc')
        decoded = SC(llr, mask);
    elseif strcmpi(decoderType, 'genie')
        decoded = oracleAssisted(llr, mask, crcPolynomial, parameter, isCRC, msg);
    elseif strcmpi(decoderType, 'scl') || strcmpi(decoderType, 'cascl')
        [decoded, iterations] = SCL(llr, mask, crcPolynomial, parameter, isCRC);
    elseif strcmpi(decoderType, 'scs') || strcmpi(decoderType, 'cascs')
        decoded = SCS(llr, mask, parameter, isCRC);
    elseif strcmpi(decoderType, 'polarfano')
        [decoded, iterations] = polarFano(llr, mask, parameter, dsnrdB, snrType);
    elseif strcmpi(decoderType, 'pbf')
        decoded = progressiveBitFlipping(llr, mask, crcPolynomial, parameter, dsnrdB, snrType);
    elseif strcmpi(decoderType, 'cahd')
        [decoded, iterations] = CA_HD(y, llr, mask, crcPolynomial, parameter, msg);
    elseif strcmpi(decoderType, 'pacfano')
        decoded = PACFano(llr, mask, convG, parameter, dsnrdB, snrType);
    elseif strcmpi(decoderType, 'pscl')
        [decoded, iterations] = PSCL(llr, mask, crcPolynomial, parameter, isCRC);
    %elseif strcmpi(decoderType, 'piao')
     %   [decoded, iterations] = PiaoRealization(y, llr, mask, crcPolynomial, inf);
    else
        disp('None such type')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%CA/Polar-Sphere Decoder%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function decoded = myCA_SD(y, mask, g, squareRadius)

    N = length(mask);
    K = length(find(mask)); % K = KI + crcWidth
    crcWidth = length(g) - 1;
    KI = K - crcWidth;
    n = log2(N);
    i = N;
    A = find(mask); 
    counter = ones(1, N);
    decoded = zeros(1, N);
    pathMetric = zeros(1, N);
    omiga = [0 1];
    uOptimal = zeros(1, N);
    %squareRadius = inf;
    GN = getFn([1 0; 1 1], n);
    fixedLowerBound = getFixedLowerBound(y);
    [I, D] = getIandD(GN, mask);
    Q = TransformPCRs(PCRs(KI, g, A), N);
    P = zeros(1, crcWidth);
    for l = 1:length(Q)
        P(l) = min(Q{l});
    end
    
    while true
       
        if mask(i) == 1 && ~ismember(i, P)
            decoded(i) = omiga(counter(i));
            if i == N
                pathMetric(i) = 0 + computeMetric(decoded, i, y, GN);
            else            
                pathMetric(i) = pathMetric(i + 1) + computeMetric(decoded, i, y, GN);
            end

        elseif mask(i) == 0

            decoded(i) = 0;
            if i == N
                pathMetric(i) = 0 + computeMetric(decoded, i, y, GN);
            else
                pathMetric(i) = pathMetric(i + 1) + computeMetric(decoded, i, y, GN);
            end

        elseif ismember(i, P)

            for l = 1:length(P)
                if P(l) == i
                    break
                end
            end
            tmp = 0;
            for t = Q{l}
                if t == P(l)
                    tmp = tmp + 0;
                else
                    tmp = tmp + decoded(t);
                end
            end
            decoded(i) = mod(tmp, 2);
            if i == N
                pathMetric(i) = 0 + computeMetric(decoded, i, y, GN);
            else
                pathMetric(i) = pathMetric(i + 1) + computeMetric(decoded, i, y, GN);
            end

        else
            disp('wtf')
        end

        if i > 1
            %metric = pathMetric(i) + fixedLowerBound(i - 1);
            if D(i) > 0
                metric = pathMetric(i) + dynamicLowerBound(i, y, I, GN, decoded);
            else
                metric = pathMetric(i) + fixedLowerBound(i - 1);
            end
        else
            metric = pathMetric(i);
        end

        if metric > squareRadius
            if mask(i) == 1 && ~ismember(i, P)
                counter(i) = counter(i) + 1;
                if counter(i) == length(omiga) + 1
                    i = i + 1;
                    if i == N + 1
                        decoded = uOptimal(mask == 1);
                        return
                    end
                    while mask(i) == 0 || (counter(i) + 1 == (length(omiga) + 1)) || ismember(i, P)
                        i = i + 1;
                        if i == N + 1
                            decoded = uOptimal(mask == 1);
                            return
                        end
                    end
                    counter(1:i - 1) = 1;
                    counter(i) = counter(i) + 1;
                end
           else
                i = i + 1;
                if i == N + 1
                    disp(111)
                    decoded = uOptimal(mask == 1);            
                    return
                end
                while mask(i) == 0 || (counter(i) + 1 == (length(omiga) + 1)) || ismember(i, P)
                      i = i + 1;
                      if i == N + 1
                         decoded = uOptimal(mask == 1);
                         return
                      end
                end
                counter(1:i - 1) = 1;
                counter(i) = counter(i) + 1;
            end
        else
            if i == 1
                uOptimal = decoded;
                squareRadius = pathMetric(1);
                %disp(squareRadius)
                i = A(1);
                counter(i) = counter(i) + 1;
                if counter(i) == length(omiga) + 1
                    i = i + 1;
                    if i == N + 1                    
                        disp(222)
                        decoded = uOptimal(mask == 1);
                        return
                    end
                    while mask(i) == 0 || (counter(i) + 1 == (length(omiga) + 1)) || ismember(i, P)
                        i = i + 1;
                        if i == N + 1
                            decoded = uOptimal(mask == 1);
                            return
                        end
                    end
                    counter(1:i - 1) = 1;
                    counter(i) = counter(i) + 1;
                end
            else
                i = i - 1;
            end
        end

    end
end
%}

%{
function decoded = myPolarSphere(y, mask, squareRadius)

    N = length(mask);
    n = log2(N);
    i = N;
    A = find(mask); 
    counter = ones(1, N);
    decoded = zeros(1, N);
    pathMetric = zeros(1, N);
    omiga = [0 1];
    uOptimal = zeros(1, N);
    %squareRadius = inf;
    GN = getFn([1 0; 1 1], n);
    fixedLowerBound = getFixedLowerBound(y);
    [I, D] = getIandD(GN, mask);
    
    while true
       
        if mask(i) == 1
            decoded(i) = omiga(counter(i));
            if i == N
                pathMetric(i) = 0 + computeMetric(decoded, i, y, GN);
            else            
                pathMetric(i) = pathMetric(i + 1) + computeMetric(decoded, i, y, GN);
            end
        elseif mask(i) == 0
            decoded(i) = 0;
            if i == N
                pathMetric(i) = 0 + computeMetric(decoded, i, y, GN);
            else
                pathMetric(i) = pathMetric(i + 1) + computeMetric(decoded, i, y, GN);
            end
        else
            disp('wtf')
        end

        if i > 1
            %metric = pathMetric(i) + fixedLowerBound(i - 1);
            if D(i) > 0
                metric = pathMetric(i) + dynamicLowerBound(i, y, I, GN, decoded);
            else
                metric = pathMetric(i) + fixedLowerBound(i - 1);
            end
        else
            metric = pathMetric(i);
        end

        if metric > squareRadius
            if mask(i) == 1
                counter(i) = counter(i) + 1;
                if counter(i) == length(omiga) + 1
                    i = i + 1;
                    if i == N + 1
                        decoded = uOptimal(mask == 1);
                        return
                    end
                    while mask(i) == 0 || (counter(i) + 1 == (length(omiga) + 1))
                        i = i + 1;
                        if i == N + 1
                            decoded = uOptimal(mask == 1);
                            return
                        end
                    end
                    counter(1:i - 1) = 1;
                    counter(i) = counter(i) + 1;
                end
           else
                i = i + 1;
                if i == N + 1
                    disp(111)
                    decoded = uOptimal(mask == 1);            
                    return
                end
                while mask(i) == 0 || (counter(i) + 1 == (length(omiga) + 1))
                      i = i + 1;
                      if i == N + 1
                         decoded = uOptimal(mask == 1);
                         return
                      end
                end
                counter(1:i - 1) = 1;
                counter(i) = counter(i) + 1;
            end
        else
            if i == 1
                uOptimal = decoded;
                squareRadius = pathMetric(1);
                %disp(squareRadius)
                i = A(1);
                counter(i) = counter(i) + 1;
                if counter(i) == length(omiga) + 1
                    i = i + 1;
                    if i == N + 1
                        %decoded = uOptimal;
                        %return
                        disp(222)
                        decoded = uOptimal(mask == 1);
                        return
                    end
                    while mask(i) == 0 || (counter(i) + 1 == (length(omiga) + 1))
                        i = i + 1;
                        if i == N + 1
                            decoded = uOptimal(mask == 1);
                            return
                        end
                    end
                    counter(1:i - 1) = 1;
                    counter(i) = counter(i) + 1;
                end
            else
                i = i - 1;
            end
        end

    end
end

function branchMetric = computeMetric(decoded, i, y, GN)
   
    N = length(decoded);
    tmp = 0;
    for j = i:N
        tmp = tmp + GN(j, i) * decoded(j);
    end
    tmp = mod(tmp, 2);
    branchMetric = (((1 - y(i)) / 2) - tmp) ^ 2;
end
%}

function [fixedLowerBound, LB] = getFixedLowerBound(y)
    yBar = (1 - y) / 2;
    fixedLowerBound = zeros(1, length(y));
    for i = 1:length(y)
        if (yBar(i) - 0) ^ 2 < (yBar(i) - 1) ^ 2
            fixedLowerBound(i) = (yBar(i) - 0) ^ 2;
        else
            fixedLowerBound(i) = (yBar(i) - 1) ^ 2;
        end
    end
    LB = fixedLowerBound;
    for i = 1:length(y)
        if i == 1
            continue
        else
            fixedLowerBound(i) = sum(fixedLowerBound(i - 1:i));
        end
    end
end

function [I, D] = getIandD(GN, mask)
    GHat = zeros(length(mask));
    GHat(mask == 1, :) = GN(mask == 1, :);   
    I = transform(GHat, mask);
    D = zeros(1, length(mask));
    for i = 1:length(I)
        if iscell(I(i))
            D(i) = length(I{i});
        end
    end
end

function I = transform(GHat, mask)
    I = {};
    for i = 1:length(mask)
        if mask(i) == 1
            
            G = GHat(1:i - 1, 1:i - 1);
            if sum(sum(G)) == 0
                I{end + 1} = [];
                continue
            else
                visits = ones(1, length(G));
                visits(sum(G) == 0) = 0;
                I1 = {};
            end

            while sum(visits) ~= 0
                a = find(visits);
                column = G(:, a(1));
                visits(a(1)) = 0;
                tmp = [];
                for j = (a(1) + 1):length(G)
                    if all(column == G(:, j))
                        if ~ismember(a(1), tmp)
                            tmp(end + 1) = a(1);
                        end 
                        tmp(end + 1) = j;
                        visits(j) = 0;
                    end
                end
                if ~isempty(tmp)
                    I1{end + 1} = sort(tmp, 'descend');
                    %I1{end + 1} = tmp;             
                end
            end
            I{end + 1} = fliplr(I1);
            %I{end + 1} = I1;
        else
            I{end + 1} = [];
        end
    end
end

function tighterBound = dynamicLowerBound(l, y, I, GN, decoded)
    N = length(y);
    IL = I{l};
    m = 0;
    for k = 1:length(IL)
        m0 = 0;
        m1 = 0;
        for i = IL{k}
            t = 0;
            for j = l:N
                t = t + GN(j, i) * decoded(j);
            end
            m0 = m0 + ((1 - y(i)) / 2 - mod(0 + t, 2)) ^ 2;
            m1 = m1 + ((1 - y(i)) / 2 - mod(1 + t, 2)) ^ 2;
        end
        m = m + min(m0, m1);
    end
    tighterBound = m;
end

function Ru = PCRs(KI, g, A)
    crcWidth = length(g) - 1;
    GC = getGC(KI, g);
    Rs = cell(1, crcWidth);
    Ru = cell(1, crcWidth);
    for l = 1:crcWidth
        for i = 1:KI
            if GC(i, l + KI) == 1
                Rs{l}(end + 1) = i;
            end
        end
        Rs{l}(end + 1) = (KI + l);
    end
    for l = 1:crcWidth
        for i = Rs{l}
            Ru{l}(end + 1) = A(i);
        end
    end
end

function Q = TransformPCRs(Ru, N)
    crcWidth = length(Ru);
    D = zeros(crcWidth, N);
    Q = cell(1, crcWidth);
    for l = 1:crcWidth
        for j = 1:N
            if ismember(j, Ru{l})
                D(l, j) = 1;
            end
        end
    end
    D = rowEchelonForm(D); % mod2 form
    for l = 1:crcWidth
        for j = 1:N
            if D(l, j) == 1
                Q{l}(end + 1) = j;
            end
        end
    end
end

%{
                    check = mod(decoded * getGC(K, g), 2);
                    if sum(check(K + 1:end)) ~= 0
                        disp('err occ')
                        pause()
                    end
%}

function [decoded, iterations] = CA_SD(y, softMess, mask, g, squareRadius, uOptimal)
    N = length(mask);
    K = length(find(mask)); % K = KI + crcWidth
    crcWidth = length(g) - 1;
    KI = K - crcWidth;
    n = log2(N);
    i = N;
    A = find(mask);
    decoded = zeros(1, N);
    %uOptimal = zeros(1, N);
    pathMetric = zeros(1, N);
    GN = getFn([1 0; 1 1], n);
    [fixedLowerBound, LB] = getFixedLowerBound(y);
    [I, D] = getIandD(GN, mask);
    Q = TransformPCRs(PCRs(KI, g, A), N);
    P = zeros(1, crcWidth);
    T = [];
    for l = 1:length(Q)
        P(l) = min(Q{l});
    end
    iterations = 0;
    initRadius = false;
    %firstReach = true;

    if initRadius
        decoded1 = SC(softMess, mask);
        %{
        check = mod(decoded1 * getGC(K, g), 2);
        if sum(check(K + 1:end)) == 0
            uOptimal(mask == 1) = decoded1;
            %decoded = decoded1;
            %return
        end
        %}
        b = decoded1(1:KI);
        uOptimal(mask == 1) = mod(b * getGC(KI, g), 2);
        xTilde = CRCPolarEncode(b, mask, g);
        eucDistance = 0;
        for ii = 1:N
            eucDistance = eucDistance + (((1 - y(ii)) / 2) - xTilde(ii)) ^ 2;
        end
        squareRadius = eucDistance;
        %fprintf('initial radius = %.2f\n', squareRadius)
    end
    %fprintf('initial radius = %.2f\n', squareRadius)

    while i <= N
        iterations = iterations + 1;
        if (i + 1) == (N + 1)
            partialPathMetric = 0;
        else
            partialPathMetric = pathMetric(i + 1);
        end

        %if i == 0 && d(decoded, i + 1, y, GN, pathMetric) <= squareRadius
        if i == 0 && partialPathMetric <= squareRadius
            %{
            if firstReach
                firstReach = false;
            end
            %}
            %squareRadius = d(decoded, i + 1, y, GN, pathMetric);
            squareRadius = partialPathMetric;
            uOptimal = decoded;
            i = i + 1;
            while mask(i) == 0 || ismember(i, P) || ismember(i, T)
                i = i + 1;
                if i == (N + 1)
                    decoded = uOptimal(mask == 1);
                    %decoded = decoded(mask == 1);
                    return
                end
            end
            decoded(i) = mod(decoded(i) + 1, 2);
            T(end + 1) = i;
            T(T < i) = [];
            pathMetric(i) = d(decoded, i, y, GN, pathMetric);
            i = i - 1;

        %elseif (d(decoded, i + 1, y, GN, pathMetric) + alpha(i, y, I, D, GN, decoded, fixedLowerBound)) <= squareRadius
        elseif (partialPathMetric + alpha(i, y, I, D, GN, decoded, fixedLowerBound, LB)) <= squareRadius

            %iterations = iterations + 1; %forward move?
            if ismember(i, P)

                for l = 1:length(P)
                    if P(l) == i
                        break
                    end
                end
                tmp = 0;
                for t = Q{l}
                    if t == P(l)
                        tmp = tmp + 0;
                    else
                        tmp = tmp + decoded(t);
                    end
                end
                decoded(i) = mod(tmp, 2);

            elseif mask(i) == 1

                tmp0 = 0;
                tmp1 = 0;
                for j = i:N
                    if j == i
                        tmp0 = tmp0 + GN(j, i) * 0;
                        tmp1 = tmp1 + GN(j, i) * 1;
                    else
                        tmp0 = tmp0 + GN(j, i) * decoded(j);
                        tmp1 = tmp1 + GN(j, i) * decoded(j);
                    end
                end
                tmp0 = mod(tmp0, 2);
                b0 = (((1 - y(i)) / 2) - tmp0) ^ 2;
                tmp1 = mod(tmp1, 2);
                b1 = (((1 - y(i)) / 2) - tmp1) ^ 2;
                if b0 < b1
                    decoded(i) = 0;
                elseif b1 < b0
                    decoded(i) = 1;
                else
                    disp('b0 = b1')
                    pause()
                end
                %{
                decoded(i) = 0;
                partialPathMetric0 = d(decoded, i, y, GN, pathMetric);
                decoded(i) = 1;
                partialPathMetric1 = d(decoded, i, y, GN, pathMetric);
                if partialPathMetric0 < partialPathMetric1
                    decoded(i) = 0;
                elseif partialPathMetric1 < partialPathMetric0
                    decoded(i) = 1;
                else
                    disp('p0 = p1')
                    pause()
                end
                %}
            elseif mask(i) == 0
                decoded(i) = 0;
            end
            pathMetric(i) = d(decoded, i, y, GN, pathMetric);
            i = i - 1;

        %elseif (d(decoded, i + 1, y, GN, pathMetric) + alpha(i, y, I, D, GN, decoded, fixedLowerBound)) > squareRadius
        elseif (partialPathMetric + alpha(i, y, I, D, GN, decoded, fixedLowerBound, LB)) > squareRadius

            i = i + 1;
            while mask(i) == 0 || ismember(i, P) || ismember(i, T)
                i = i + 1;
                if i == (N + 1)
                    decoded = uOptimal(mask == 1);
                    %decoded = decoded(mask == 1);
                    %firstReach
                    return
                end
            end
            decoded(i) = mod(decoded(i) + 1, 2);
            T(end + 1) = i;
            T(T < i) = [];
            pathMetric(i) = d(decoded, i, y, GN, pathMetric);
            i = i - 1;

        end
    end
    decoded = uOptimal(mask == 1);
    %decoded = decoded(mask == 1);
end

function partialPathMetric = d(decoded, i, y, GN, pathMetric)
    N = length(decoded);
    if i == (N + 1)
        partialPathMetric = 0;
    else
        tmp = 0;
        for j = i:N
            tmp = tmp + GN(j, i) * decoded(j);
        end
        tmp = mod(tmp, 2);
        if i == N
            partialPathMetric = 0 + (((1 - y(i)) / 2) - tmp) ^ 2;
        else
            partialPathMetric = pathMetric(i + 1) + (((1 - y(i)) / 2) - tmp) ^ 2;
        end
        %{
        for k = i:N
            tmp = 0;
            for j = k:N
                tmp = tmp + GN(j, k) * decoded(j);
            end
            tmp = mod(tmp, 2);
            partialPathMetric = partialPathMetric + (((1 - y(k)) / 2) - tmp) ^ 2;
        end
        %}
    end
end

function DLB = alpha(i, y, I, D, GN, decoded, fixedLowerBound, LB) % DLB = dynamic lower bound
    N = length(y);
    i = i + 1;
    if i == N + 1
        DLB = fixedLowerBound(i - 1);
    else
        if D(i) > 0
            DLB = dynamicLowerBound(i, y, I, GN, decoded);
            tmp = [];
            for ii = 1:length(I{i})
                tmp = [tmp I{i}{ii}];
            end
            for ii = i - 1:-1:1
                if ~ismember(ii, tmp)
                    DLB = DLB + LB(ii);
                end
            end
            %{
            if abs(DLB - fixedLowerBound(i - 1)) > 1e-10 && (DLB < fixedLowerBound(i - 1))
                %DLB - fixedLowerBound(i - 1)
                disp('err')
                pause()
            end
            %}
        else
            if i == 1
                DLB = 0;
            else
                DLB = fixedLowerBound(i - 1);
            end
        end
    end
end

function [decoded, iterations] = polarSphere(y, softMess, mask, squareRadius)
    N = length(mask);
    K = length(find(mask)); % K = KI + crcWidth
    n = log2(N);
    i = N;
    decoded = zeros(1, N);
    uOptimal = zeros(1, N);
    pathMetric = zeros(1, N);
    %squareRadius = inf;
    GN = getFn([1 0; 1 1], n);
    [fixedLowerBound, LB] = getFixedLowerBound(y);
    [I, D] = getIandD(GN, mask);
    T = [];
    iterations = 0;
    isInit = true;
    if isInit
        msg = SC(softMess, mask);
        uOptimal(mask == 1) = msg;
        xTilde = polarEncode(msg, mask);
        eucDistance = 0;
        for ii = 1:N
            eucDistance = eucDistance + (((1 - y(ii)) / 2) - xTilde(ii)) ^ 2;
        end
        squareRadius = eucDistance;
    end

    while i <= N
        iterations = iterations + 1;
        if (i + 1) == (N + 1)
            partialPathMetric = 0;
        else
            partialPathMetric = pathMetric(i + 1);
        end

        if i == 0 && partialPathMetric <= squareRadius
            
            %squareRadius = d(decoded, i + 1, y, GN, pathMetric);
            squareRadius = partialPathMetric;
            uOptimal = decoded;
            i = i + 1;
            while mask(i) == 0 || ismember(i, T)
                i = i + 1;
                if i == (N + 1)
                    decoded = uOptimal;
                    decoded = decoded(mask == 1);
                    %squareRadius
                    return
                end
            end
            decoded(i) = mod(decoded(i) + 1, 2);
            T(end + 1) = i;
            T(T < i) = [];
            pathMetric(i) = d(decoded, i, y, GN, pathMetric);
            i = i - 1;

        elseif (partialPathMetric + alpha(i, y, I, D, GN, decoded, fixedLowerBound, LB)) <= squareRadius
        %elseif partialPathMetric <= squareRadius

            if mask(i) == 1

                tmp0 = 0;
                tmp1 = 0;
                for j = i:N
                    if j == i
                        tmp0 = tmp0 + GN(j, i) * 0;
                        tmp1 = tmp1 + GN(j, i) * 1;
                    else
                        tmp0 = tmp0 + GN(j, i) * decoded(j);
                        tmp1 = tmp1 + GN(j, i) * decoded(j);
                    end
                end
                tmp0 = mod(tmp0, 2);
                b0 = (((1 - y(i)) / 2) - tmp0) ^ 2;
                tmp1 = mod(tmp1, 2);
                b1 = (((1 - y(i)) / 2) - tmp1) ^ 2;
                if b0 < b1
                    decoded(i) = 0;
                elseif b1 < b0
                    decoded(i) = 1;
                else
                    disp('b0 = b1')
                    pause()
                end
            elseif mask(i) == 0
                decoded(i) = 0;
            end
            pathMetric(i) = d(decoded, i, y, GN, pathMetric);
            i = i - 1;

        elseif (partialPathMetric + alpha(i, y, I, D, GN, decoded, fixedLowerBound, LB)) > squareRadius
        %elseif partialPathMetric > squareRadius

            i = i + 1;
            while mask(i) == 0 || ismember(i, T)
                i = i + 1;
                if i == (N + 1)
                    decoded = uOptimal;
                    decoded = decoded(mask == 1);
                    %squareRadius
                    %iterations
                    return
                end
            end
            decoded(i) = mod(decoded(i) + 1, 2);
            T(end + 1) = i;
            T(T < i) = [];
            pathMetric(i) = d(decoded, i, y, GN, pathMetric);
            i = i - 1;

        end
    end
    decoded = uOptimal;
    decoded = decoded(mask == 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%Polar-SC/SCL/SCS/Genie Decoder%%%%%%%%%%%%%%%%%%%%%%%%%

function decoded = oracleAssisted(softMess, mask, g, parameter, isCRC, msg)
    N = length(mask);
    LLRs = zeros(1, 2 * N - 1);
    LLRs(N:end) = softMess;
    Bits = zeros(2, N - 1);
    decoded = zeros(1, N);
    correct = zeros(1, N);
    TMax = parameter;
    T = 0;
    if isCRC
        correct(mask == 1) = mod(msg * getGC(length(msg), g), 2);
    else
        correct(mask == 1) = msg;
    end

    for i = 1:N
        LLRs = updateLLRs(i, N, LLRs, Bits);
        if mask(i) == 1
            
            if LLRs(1) > 0
                decoded(i) = 0;
            elseif LLRs(1) < 0
                decoded(i) = 1;
            else
                disp('warning')
            end

        else
            decoded(i) = 0;
        end
        if decoded(i) ~= correct(i) && T < TMax
            decoded(i) = correct(i);
            T = T + 1;
        end
        Bits = updateBits(i, N, Bits, decoded(i));
    end
    decoded = decoded(mask == 1);
end

function decoded = SC(softMess, mask)
    N = length(mask);
    %n = log2(N);
    LLRs = zeros(1, 2 * N - 1);
    LLRs(N:end) = softMess;
    Bits = zeros(2, N - 1);
    decoded = zeros(1, N);
    for i = 1:N
        LLRs = updateLLRs(i, N, LLRs, Bits);
        if mask(i) == 1
            
            if LLRs(1) > 0
                decoded(i) = 0;
            elseif LLRs(1) < 0
                decoded(i) = 1;
            else
                disp('warning')
            end

        else
            decoded(i) = 0;
        end
        Bits = updateBits(i, N, Bits, decoded(i));
    end
    decoded = decoded(mask == 1);
end

function [decoded, iterations] = SCL(softMess, mask, g, L, isCRC)
    N = length(mask);
    LLRs = zeros(L, 2 * N - 1);
    LLRs(1, N:end) = softMess;
    Bits = zeros(2, N - 1, L);
    decodeds = zeros(L, N);
    pathMetrics = zeros(1, L);
    pathNum = 1;
    iterations = 0;
    for i = 1:N
        for j = 1:pathNum
            LLRs(j, :) = updateLLRs(i, N, LLRs(j, :), Bits(:, :, j));
            iterations = iterations + 1;
        end

        if mask(i) == 1

            if 2 * pathNum <= L
    
                for i1 = 1:pathNum
                    %{
                    penalty = abs(LLRs(i1, 1));
                    if penalty == 0
                        disp('llr warning')
                        pause()
                    end
                    branchMetric0 = (0 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                    branchMetric1 = (1 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                    %}
                    Li = LLRs(i1, 1) / log(2);
                    branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                    branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                    pathMetrics(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                    pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                    %pathMetrics(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                    %decodeds(i1, i) = 0;
                    decodeds(i1 + pathNum, :) = decodeds(i1, :);
                    decodeds(i1 + pathNum, i) = 1;
                    LLRs(i1 + pathNum, :) = LLRs(i1, :);
                    Bits(:, :, i1 + pathNum) = Bits(:, :, i1);
                end
                pathNum = 2 * pathNum;
            else
         
                PM = zeros(1, 2 * pathNum);
                discardTag = ones(1, 2 * pathNum);
                pathStates = zeros(1, pathNum); % 1 for both of paths that should be discarded; 2 for one of paths that should be discarded; 3 for both of paths that should be retained
                for i1 = 1:pathNum
                    %{
                    penalty = abs(LLRs(i1, 1));
                    if penalty == 0
                        disp('llr warning')
                        pause()
                    end
                    branchMetric0 = (0 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                    branchMetric1 = (1 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                    %}
                    Li = LLRs(i1, 1) / log(2);
                    branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                    branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                    PM(i1) = pathMetrics(i1) + branchMetric0;
                    PM(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                    %pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                end

                %iterations = iterations + 1;
                %[~, indices] = sort(PM);
                [~, indices] = sort(PM, 'descend');
                discardTag(indices(1:L)) = 0;
                disIdxRecord = [];
                for i1 = 1:pathNum
                    if discardTag(i1) == 1 && discardTag(i1 + pathNum) == 1
                        pathStates(i1) = 1;
                        disIdxRecord(end + 1) = i1;
                    elseif discardTag(i1) == 1 || discardTag(i1 + pathNum) == 1
                        pathStates(i1) = 2;
                    elseif discardTag(i1) == 0 && discardTag(i1 + pathNum) == 0
                        pathStates(i1) = 3;
                    else
                        disp('state warning')
                    end
                end
                
                for i1 = 1:length(pathStates)
                    if pathStates(i1) == 2
                        if discardTag(i1) == 1
                            pathMetrics(i1) = PM(i1 + pathNum);
                            decodeds(i1, i) = 1;
                        elseif discardTag(i1 + pathNum) == 1
                            pathMetrics(i1) = PM(i1);
                        end
                    elseif pathStates(i1) == 3
                        pathMetrics(i1) = PM(i1);
                        pathMetrics(disIdxRecord(1)) = PM(i1 + pathNum);
                        decodeds(disIdxRecord(1), :) = decodeds(i1, :);
                        decodeds(disIdxRecord(1), i) = 1;
                        LLRs(disIdxRecord(1), :) = LLRs(i1, :);
                        Bits(:, :, disIdxRecord(1)) = Bits(:, :, i1);
                        disIdxRecord(1) = [];
                    end
                end
            end

        else

            for j = 1:pathNum
                decodeds(j, i) = 0;
                %{
                penalty = abs(LLRs(j, 1));
                if penalty == 0
                    disp('llr warning1')
                    pause()
                end
                branchMetric0 = (decodeds(j, i) ~= 1 / 2 * (1 - sign(LLRs(j, 1)))) * penalty;
                %}
                Li = LLRs(j, 1) / log(2);
                branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                pathMetrics(j) = pathMetrics(j) + branchMetric0;
            end

        end
        for j = 1:pathNum
            Bits(:, :, j) = updateBits(i, N, Bits(:, :, j), decodeds(j, i));
        end
    end

    %[~, I] = sort(pathMetrics);
    [~, I] = sort(pathMetrics, 'descend');

    if isCRC
        %disp('hi')
        decodeds = decodeds(I, :);
        crcWidth = length(g) - 1;
        K = length(find(mask)); % K = KI + crcWidth
        KI = K - crcWidth;
        for l = 1:size(decodeds, 1)
            best = decodeds(l, :);
            best = best(mask == 1);
            check = mod(best * getGC(K, g), 2);
            if sum(check(K + 1:end)) == 0
                %decoded = best(1:KI);
                decoded = best;
                return
            end
        end
        decoded = decodeds(1, :);
        decoded = decoded(mask == 1);
        %decoded = decoded(1:KI);
    else
        decoded = decodeds(I(1), :);
        decoded = decoded(mask == 1);
    end
end

function decoded = SCS(softMess, mask, D, isCRC)
    N = length(mask);
    LLRs = zeros(D, 2 * N - 1);
    LLRs(1, N:end) = softMess;
    Bits = zeros(2, N - 1, D);
    decodeds = zeros(D, N);
    pathMetrics = zeros(1, D);
    i = zeros(1, D);
    T = 1;
    flag = true;
    while flag
        
        if mask(i(T) + 1) == 1
           
           i(T) = i(T) + 1;
           LLRs(T, :) = updateLLRs(i(T), N, LLRs(T, :), Bits(:, :, T));
           penalty = abs(LLRs(T, 1));
           if penalty == 0
               disp('llr warning')
               pause()
           end
           branchMetric0 = (0 ~= 1 / 2 * (1 - sign(LLRs(T, 1)))) * penalty;
           branchMetric1 = (1 ~= 1 / 2 * (1 - sign(LLRs(T, 1)))) * penalty;

           if T <= (D - 1)
               
               i(T + 1) = i(T);
               pathMetrics(T + 1) = pathMetrics(T) + branchMetric1;
               pathMetrics(T) = pathMetrics(T) + branchMetric0;
               decodeds(T + 1, :) = decodeds(T, :);
               decodeds(T + 1, i(T)) = 1;
               LLRs(T + 1, :) = LLRs(T, :);
               Bits(:, :, T) = updateBits(i(T), N, Bits(:, :, T), decodeds(T, i(T)));
               Bits(:, :, T + 1) = updateBits(i(T), N, Bits(:, :, T), decodeds(T + 1, i(T)));
               T = T + 1;

               [pathMetrics(1:T), I] = sort(pathMetrics(1:T), 'descend');
               i(1:T) = i(I);
               decodeds(1:T, :) = decodeds(I, :);
               LLRs(1:T, :) = LLRs(I, :);
               Bits(:, :, 1:T) = Bits(:, :, I);
           else
               
               i(T + 1) = i(T);
               pathMetrics(T + 1) = pathMetrics(T) + branchMetric1;
               pathMetrics(T) = pathMetrics(T) + branchMetric0;
               decodeds(T + 1, :) = decodeds(T, :);
               decodeds(T + 1, i(T)) = 1;
               LLRs(T + 1, :) = LLRs(T, :);
               Bits(:, :, T) = updateBits(i(T), N, Bits(:, :, T), decodeds(T, i(T)));
               Bits(:, :, T + 1) = updateBits(i(T), N, Bits(:, :, T), decodeds(T + 1, i(T)));
               T = T + 1;
               
               [pathMetrics(1:T), I] = sort(pathMetrics(1:T), 'descend');
               
               pathMetrics(1) = [];
               i(1:T) = i(I);
               i(1) = [];
               decodeds(1:T, :) = decodeds(I, :);
               decodeds(1, :) = [];
               LLRs(1:T, :) = LLRs(I, :);
               LLRs(1, :) = [];
               Bits(:, :, 1:T) = Bits(:, :, I);
               Bits(:, :, 1) = [];
               T = T - 1;
           end

        else
            i(T) = i(T) + 1;
            LLRs(T, :) = updateLLRs(i(T), N, LLRs(T, :), Bits(:, :, T));
            decodeds(T, i(T)) = 0;
            penalty = abs(LLRs(T, 1));
            if penalty == 0
                disp('llr warning1')
                pause()
            end
            branchMetric0 = (decodeds(T, i(T)) ~= 1 / 2 * (1 - sign(LLRs(T, 1)))) * penalty;
            pathMetrics(T) = pathMetrics(T) + branchMetric0;
            Bits(:, :, T) = updateBits(i(T), N, Bits(:, :, T), decodeds(T, i(T)));
            [pathMetrics(1:T), I] = sort(pathMetrics(1:T), 'descend');
            
            i(1:T) = i(I);
            decodeds(1:T, :) = decodeds(I, :);
            LLRs(1:T, :) = LLRs(I, :);
            Bits(:, :, 1:T) = Bits(:, :, I);
        end
        
        if i(T) == N
            if isCRC
                disp('hi')
            else
                decoded = decodeds(T, :);
                decoded = decoded(mask == 1);
                %flag = false
                return
            end
        end
        
    end
end
%{
function [decoded, iterations] = PSCL(softMess, mask, g, L, isCRC)
    N = length(mask);
    LLRs = zeros(L, 2 * N - 1);
    LLRs(1, N:end) = softMess;
    Bits = zeros(2, N - 1, L);
    decodeds = zeros(L, N);
    pathMetrics = zeros(1, L);
    pathNum = 1;
    mT = -inf;
    iterations = 0;

    for i = 1:N
        for j = 1:pathNum
            LLRs(j, :) = updateLLRs(i, N, LLRs(j, :), Bits(:, :, j));
        end

        if mask(i) == 1
            duplicateIdx = [];
            swappingIdx = [];
            retainIdx = [];
            discardIdx = [];
            for i1 = 1:pathNum
                Li = LLRs(i1, 1) / log(2);
                branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                if branchMetric0 >= mT && branchMetric1 >= mT
                    duplicateIdx(end + 1) = i1;
                elseif branchMetric0 < mT && branchMetric1 >= mT
                    swappingIdx(end + 1) = i1;
                elseif branchMetric0 >= mT && branchMetric1 < mT
                    retainIdx(end + 1) = i1;
                elseif branchMetric0 < mT && branchMetric1 < mT
                    discardIdx(end + 1) = i1;
                end
            end
            
            if (2 * length(duplicateIdx) + length(swappingIdx) + length(retainIdx) - length(discardIdx)) <= L
                ii = 1;
                for i1 = 1:pathNum
                    Li = LLRs(i1, 1) / log(2);
                    branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                    branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                    if ismember(i1, duplicateIdx) && ~isempty(discardIdx)
                        pathMetrics(discardIdx(1)) = pathMetrics(i1);
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        pathMetrics(discardIdx(1)) = pathMetrics(discardIdx(1)) + branchMetric1;
                        decodeds(discardIdx(1), :) = decodeds(i1, :);
                        decodeds(discardIdx(1), i) = 1;
                        LLRs(discardIdx(1), :) = LLRs(i1, :);
                        Bits(:, :, discardIdx(1)) = Bits(:, :, i1);
                        discardIdx(1) = [];
                    elseif ismember(i1, duplicateIdx)
                        %{
                        pathMetrics(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        %pathMetrics(i1 + pathNum) = pathMetrics(i1 + pathNum) + branchMetric1;
                        decodeds(i1 + pathNum, :) = decodeds(i1, :);
                        decodeds(i1 + pathNum, i) = 1;
                        LLRs(i1 + pathNum, :) = LLRs(i1, :);
                        Bits(:, :, i1 + pathNum) = Bits(:, :, i1);
                        %}
                        
                        pathMetrics(ii + pathNum) = pathMetrics(i1) + branchMetric1;
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        decodeds(ii + pathNum, :) = decodeds(i1, :);
                        decodeds(ii + pathNum, i) = 1;
                        LLRs(ii + pathNum, :) = LLRs(i1, :);
                        Bits(:, :, ii + pathNum) = Bits(:, :, i1);
                        ii = ii + 1;
                        
                    elseif ismember(i1, swappingIdx)
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric1;
                        decodeds(i1, i) = 1;
                    elseif ismember(i1, retainIdx)
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                    end
                end
                pathNum = 2 * length(duplicateIdx) + length(swappingIdx) + length(retainIdx) - length(discardIdx);
            else
            
                tmpMetrics = zeros(1, 2 * length(duplicateIdx) + length(swappingIdx) + length(retainIdx) - length(discardIdx));
                for i1 = 1:pathNum
                    Li = LLRs(i1, 1) / log(2);
                    branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                    branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                    if ismember(i1, duplicateIdx)
                        pathMetric0 = pathMetrics(i1) + branchMetric0;
                        pathMetric1 = pathMetrics(i1) + branchMetric1;
                        tmpMetrics(i1) = pathMetric0;
                        tmpMetrics(i1 + pathNum) = pathMetric1;
                    elseif ismember(i1, swappingIdx)
                        pathMetric1 = pathMetrics(i1) + branchMetric1;
                        tmpMetrics(i1) = pathMetric1;
                    elseif ismember(i1, retainIdx)
                        pathMetric0 = pathMetrics(i1) + branchMetric0;
                        tmpMetrics(i1) = pathMetric0;
                    end
                end

                iterations = iterations + 1;
                [~, I] = sort(tmpMetrics, 'descend');
                survive = I(1:L);
                kill = I(L + 1:end);
                %kill(kill > L) = [];
                kill(kill > pathNum) = [];
                %{
                for i1 = survive
                    if i1 > L && ismember(i1 - pathNum, duplicateIdx)
                        pathMetrics(kill(1)) = tmpMetrics(i1);
                        decodeds(kill(1), :) = decodeds(i1 - pathNum, :);
                        decodeds(kill(1), i) = 1;
                        LLRs(kill(1), :) = LLRs(i1 - pathNum, :);
                        Bits(:, :, kill(1)) = Bits(:, :, i1 - pathNum);
                        kill(1) = [];
                    else
                        if ismember(i1, swappingIdx)
                            pathMetrics(i1) = tmpMetrics(i1);
                            decodeds(i1, i) = 1;
                        else
                            pathMetrics(i1) = tmpMetrics(i1);
                        end
                    end
                end
                %}
                ii = 1;
                for i1 = survive
                    if i1 > pathNum
                        if ismember(i1 - pathNum, survive)
                            if isempty(kill)
                                pathMetrics(pathNum + ii) = tmpMetrics(i1);
                                decodeds(pathNum + ii, :) = decodeds(i1 - pathNum, :);
                                decodeds(pathNum + ii, i) = 1;
                                LLRs(pathNum + ii, :) = LLRs(i1 - pathNum, :);
                                Bits(:, :, pathNum + ii) = Bits(:, :, i1 - pathNum);
                                ii = ii + 1;
                            else
                                pathMetrics(kill(1)) = tmpMetrics(i1);
                                decodeds(kill(1), :) = decodeds(i1 - pathNum, :);
                                decodeds(kill(1), i) = 1;
                                LLRs(kill(1), :) = LLRs(i1 - pathNum, :);
                                Bits(:, :, kill(1)) = Bits(:, :, i1 - pathNum);
                                kill(1) = [];
                            end
                        else
                            pathMetrics(i1 - pathNum) = tmpMetrics(i1);
                            decodeds(i1 - pathNum, i) = 1;
                        end
                    else
                        if ismember(i1, swappingIdx)
                            pathMetrics(i1) = tmpMetrics(i1);
                            decodeds(i1, i) = 1;
                        elseif ismember(i1, retainIdx)
                            pathMetrics(i1) = tmpMetrics(i1);
                        elseif ismember(i1, duplicateIdx)
                            pathMetrics(i1) = tmpMetrics(i1);
                        end
                    end
                end
                pathNum = L;
            end
        else

            for j = 1:pathNum
                decodeds(j, i) = 0;
                %penalty = abs(LLRs(j, 1));
                %branchMetric0 = (decodeds(j, i) ~= 1 / 2 * (1 - sign(LLRs(j, 1)))) * penalty;
                Li = LLRs(j, 1) / log(2);
                branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                pathMetrics(j) = pathMetrics(j) + branchMetric0;
            end

        end
        for j = 1:pathNum
            Bits(:, :, j) = updateBits(i, N, Bits(:, :, j), decodeds(j, i));
        end
    end

    [~, I] = sort(pathMetrics, 'descend');
    if isCRC
        decodeds = decodeds(I, :);
        crcWidth = length(g) - 1;
        K = length(find(mask)); % K = KI + crcWidth
        KI = K - crcWidth;
        for l = 1:size(decodeds, 1)
            best = decodeds(l, :);
            best = best(mask == 1);
            check = mod(best * getGC(K, g), 2);
            if sum(check(K + 1:end)) == 0
                decoded = best;
                return
            end
        end
        decoded = decodeds(1, :);
        decoded = decoded(mask == 1);
    else
        decoded = decodeds(I(1), :);
        decoded = decoded(mask == 1);
    end
end
%}

function [decoded, iterations] = PSCL(softMess, mask, g, L, isCRC)
    N = length(mask);
    LLRs = zeros(L, 2 * N - 1);
    LLRs(1, N:end) = softMess;
    Bits = zeros(2, N - 1, L);
    decodeds = zeros(L, N);
    pathMetrics = zeros(1, L);
    pathNum = 1;
    mT = -inf;
    iterations = 0;

    for i = 1:N
        for j = 1:pathNum
            LLRs(j, :) = updateLLRs(i, N, LLRs(j, :), Bits(:, :, j));
        end

        if mask(i) == 1
            duplicateIdx = [];
            swappingIdx = [];
            retainIdx = [];
            discardIdx = [];
            for i1 = 1:pathNum
                Li = LLRs(i1, 1) / log(2);
                branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                if branchMetric0 >= mT && branchMetric1 >= mT
                    duplicateIdx(end + 1) = i1;
                elseif branchMetric0 < mT && branchMetric1 >= mT
                    swappingIdx(end + 1) = i1;
                elseif branchMetric0 >= mT && branchMetric1 < mT
                    retainIdx(end + 1) = i1;
                elseif branchMetric0 < mT && branchMetric1 < mT
                    discardIdx(end + 1) = i1;
                end
            end
            
            if (2 * length(duplicateIdx) + length(swappingIdx) + length(retainIdx) - length(discardIdx)) <= L
                ii = 1;
                for i1 = 1:pathNum
                    Li = LLRs(i1, 1) / log(2);
                    branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                    branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                    if ismember(i1, duplicateIdx) && ~isempty(discardIdx)
                        pathMetrics(discardIdx(1)) = pathMetrics(i1);
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        pathMetrics(discardIdx(1)) = pathMetrics(discardIdx(1)) + branchMetric1;
                        decodeds(discardIdx(1), :) = decodeds(i1, :);
                        decodeds(discardIdx(1), i) = 1;
                        LLRs(discardIdx(1), :) = LLRs(i1, :);
                        Bits(:, :, discardIdx(1)) = Bits(:, :, i1);
                        discardIdx(1) = [];
                    elseif ismember(i1, duplicateIdx)
                        %{
                        pathMetrics(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        %pathMetrics(i1 + pathNum) = pathMetrics(i1 + pathNum) + branchMetric1;
                        decodeds(i1 + pathNum, :) = decodeds(i1, :);
                        decodeds(i1 + pathNum, i) = 1;
                        LLRs(i1 + pathNum, :) = LLRs(i1, :);
                        Bits(:, :, i1 + pathNum) = Bits(:, :, i1);
                        %}
                        
                        pathMetrics(ii + pathNum) = pathMetrics(i1) + branchMetric1;
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        decodeds(ii + pathNum, :) = decodeds(i1, :);
                        decodeds(ii + pathNum, i) = 1;
                        LLRs(ii + pathNum, :) = LLRs(i1, :);
                        Bits(:, :, ii + pathNum) = Bits(:, :, i1);
                        ii = ii + 1;
                        
                    elseif ismember(i1, swappingIdx)
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric1;
                        decodeds(i1, i) = 1;
                    elseif ismember(i1, retainIdx)
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                    end
                end
                pathNum = 2 * length(duplicateIdx) + length(swappingIdx) + length(retainIdx) - length(discardIdx);
            else
                totalNum = 2 * length(duplicateIdx) + length(swappingIdx) + length(retainIdx) - length(discardIdx);
                tmpMetrics = zeros(1, totalNum);
                tmpDecodeds = zeros(totalNum, N);
                tmpLLRs = zeros(totalNum, 2 * N - 1);
                tmpBits = zeros(2, N - 1, totalNum);
                for i1 = 1:pathNum
                    Li = LLRs(i1, 1) / log(2);
                    branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                    branchMetric1 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 1));
                    if ismember(i1, duplicateIdx)

                        tmpMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        tmpMetrics(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                        tmpDecodeds(i1, :) = decodeds(i1, :);
                        tmpDecodeds(i1 + pathNum, :) = decodeds(i1, :);
                        tmpDecodeds(i1 + pathNum, i) = 1;
                        tmpLLRs(i1, :) = LLRs(i1, :);
                        tmpLLRs(i1 + pathNum, :) = LLRs(i1, :);
                        tmpBits(:, :, i1) = Bits(:, :, i1);
                        tmpBits(:, :, i1 + pathNum) = Bits(:, :, i1);
                    elseif ismember(i1, swappingIdx)
 
                        tmpMetrics(i1) = pathMetrics(i1) + branchMetric1;
                        tmpDecodeds(i1, :) = decodeds(i1, :);
                        tmpDecodeds(i1, i) = 1;
                        tmpLLRs(i1, :) = LLRs(i1, :);
                        tmpBits(:, :, i1) = Bits(:, :, i1);
                    elseif ismember(i1, retainIdx)
 
                        tmpMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        tmpDecodeds(i1, :) = decodeds(i1, :);
                        tmpLLRs(i1, :) = LLRs(i1, :);
                        tmpBits(:, :, i1) = Bits(:, :, i1);
                    end
                end

                iterations = iterations + 1;
                [~, I] = sort(tmpMetrics, 'descend');
                pathMetrics(1:L) = tmpMetrics(I(1:L));
                decodeds(1:L, :) = tmpDecodeds(I(1:L), :);
                LLRs(1:L, :) = tmpLLRs(I(1:L), :);
                Bits(:, :, 1:L) = tmpBits(:, :, I(1:L));
                pathNum = L;
            end
        else

            for j = 1:pathNum
                decodeds(j, i) = 0;
                Li = LLRs(j, 1) / log(2);
                branchMetric0 = 1 - log2(1 + 2 ^ (-Li * (-1) ^ 0));
                pathMetrics(j) = pathMetrics(j) + branchMetric0;
            end

        end
        for j = 1:pathNum
            Bits(:, :, j) = updateBits(i, N, Bits(:, :, j), decodeds(j, i));
        end
    end

    [~, I] = sort(pathMetrics, 'descend');
    if isCRC
        decodeds = decodeds(I, :);
        crcWidth = length(g) - 1;
        K = length(find(mask)); % K = KI + crcWidth
        KI = K - crcWidth;
        for l = 1:size(decodeds, 1)
            best = decodeds(l, :);
            best = best(mask == 1);
            check = mod(best * getGC(K, g), 2);
            if sum(check(K + 1:end)) == 0
                decoded = best;
                return
            end
        end
        decoded = decodeds(1, :);
        decoded = decoded(mask == 1);
    else
        decoded = decodeds(I(1), :);
        decoded = decoded(mask == 1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%Polar-Fano Decoder%%%%%%%%%%%%%%%%%%%%%%%%%

function [decoded, iterations] = polarFano(softMess, mask, parameter, dsnrdB, snrType) %caution!
    N = length(mask);
    K = length(find(mask)); % K = KI + crcWidth
    i = 0;
    LLRs = zeros(1, 2 * N - 1);
    LLRs(N:end) = softMess;
    Bits = zeros(2, N - 1);
    decoded = zeros(1, N);
    pathMetrics = zeros(1, N);
    followOtherBranch = zeros(1, N);
    delta = parameter;
    %dsnrdB = 2.0;
    %pe = PE(N, K, dsnrdB, snrType);
    I = symmetricCapacity(N, K, dsnrdB, snrType);
    T = 0;
    iterations = 0;
    flag = true;
    while flag
        %iterations = iterations + 1;
        pointer = i + 1;
        if mask(pointer) == 1
            LLRs = updateLLRs(pointer, N, LLRs, Bits);
            u0 = 0;
            u1 = 1;
            %branchMetric0 = log((exp(LLRs(1)) / (exp(LLRs(1)) + 1)) / (1 - pe(pointer)));
            %branchMetric1 = log((1 / (exp(LLRs(1)) + 1)) / (1 - pe(pointer)));
            branchMetric0 = 1 - log2(1 + exp(LLRs(1)) ^ (-(1 - 2 * u0))) - I(pointer);
            branchMetric1 = 1 - log2(1 + exp(LLRs(1)) ^ (-(1 - 2 * u1))) - I(pointer);
            if pointer == 1
                pathMetric0 = 0 + branchMetric0;
                pathMetric1 = 0 + branchMetric1;
            else
                pathMetric0 = pathMetrics(pointer - 1) + branchMetric0;
                pathMetric1 = pathMetrics(pointer - 1) + branchMetric1;
            end
            if followOtherBranch(pointer) == 0
                pathMetrics(pointer) = max([pathMetric0, pathMetric1]);
            else
                pathMetrics(pointer) = min([pathMetric0, pathMetric1]);
            end
        else
            LLRs = updateLLRs(pointer, N, LLRs, Bits);
            u0 = 0;
            %branchMetric0 = log((exp(LLRs(1)) / (exp(LLRs(1)) + 1)) / (1 - pe(pointer)));
            branchMetric0 = 1 - log2(1 + exp(LLRs(1)) ^ (-(1 - 2 * u0))) - I(pointer);
            if pointer == 1
                pathMetric0 = 0 + branchMetric0;
            else
                pathMetric0 = pathMetrics(pointer - 1) + branchMetric0;
            end
            pathMetrics(pointer) = pathMetric0;
        end

        if pathMetrics(pointer) >= T
            iterations = iterations + 1;
            i = i + 1;
            if i == 0
                muPre = -inf;
            elseif i == 1
                muPre = 0;
            else
                muPre = pathMetrics(i - 1);
            end
            if muPre < (T + delta)
                while (T + delta) <= pathMetrics(i)
                    T = T + delta;
                end
            end
            if pathMetrics(i) == pathMetric0
                decoded(i) = 0;
            else
                decoded(i) = 1;
            end
            Bits = updateBits(i, N, Bits, decoded(i));
            if i == N
                decoded = decoded(mask == 1);
                %flag = false
                return
            end
            followOtherBranch(i + 1) = 0;
        else
            j = i;
            [i, T, followOtherBranch] = moveBack(i, T, delta, mask, pathMetrics, followOtherBranch);
            if j ~= i
                for ii = 1:i
                    LLRs = updateLLRs(ii, N, LLRs, Bits);
                    Bits = updateBits(ii, N, Bits, decoded(ii));
                end
            end
        end

    end
end

function [i, T, followOtherBranch] = moveBack(i, T, delta, mask, pathMetrics, followOtherBranch)
    while true
        if i == 0
            muPre = -inf;
        elseif i == 1
            muPre = 0;
        else
            muPre = pathMetrics(i - 1);
        end
        if muPre < T
            T = T - delta;
            followOtherBranch(i + 1) = 0;
            return
        else
            if (followOtherBranch(i) + 1) == 2 || mask(i) == 0
                i = i - 1;
            else
                followOtherBranch(i) = followOtherBranch(i) + 1;
                i = i - 1;
                return
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%Progressive Bit-Flipping Decoder%%%%%%%%%%%%%%%%%%%%%%%%%

function [decoded, llrs] = scFlip(softMess, mask, flippingPosition)
    N = length(mask);
    %n = log2(N);
    LLRs = zeros(1, 2 * N - 1);
    LLRs(N:end) = softMess;
    llrs = zeros(1, N);
    Bits = zeros(2, N - 1);
    decoded = zeros(1, N);
    for i = 1:N
        LLRs = updateLLRs(i, N, LLRs, Bits);
        llrs(i) = LLRs(1);
        if mask(i) == 1
            
            if LLRs(1) > 0
                decoded(i) = 0;
            elseif LLRs(1) < 0
                decoded(i) = 1;
            else
                disp('warning')
            end

        else
            decoded(i) = 0;
        end
        if ismember(i, flippingPosition)
            decoded(i) = mod(decoded(i) + 1, 2);
        end
        Bits = updateBits(i, N, Bits, decoded(i));
    end
    decoded = decoded(mask == 1);
end

function decoded = progressiveBitFlipping(softMess, mask, g, parameter, dsnrdB, snrType)
    N = length(mask);
    K = length(find(mask));
    crcWidth = length(g) - 1;
    KI = K - crcWidth;
    maxLevel = parameter;
    l = 0;
    w = [0.35, 0.5, 0.25];
    gamaLeft = 3.6;
    gamaRight = 2;
    mu = GA(N, K, dsnrdB, snrType);
    S = cell(1, maxLevel);
    isPruningTechnique = false;
    while (l <= maxLevel)
        if l == 0
            [decoded, llrs] = scFlip(softMess, mask, 0);
            check = mod(decoded * getGC(K, g), 2);
            if sum(check(K + 1:end)) == 0
                return
            elseif sum(check(K + 1:end)) ~=0 && (l + 1) <= maxLevel             
                CS = modifyCriticalSet(mask, 0);
                M = abs(llrs(CS) ./ sqrt(mu(CS)));
                [~, I] = sort(M);
                CS = CS(I);
                for u = CS
                    S{l + 1}(end + 1) = {u};
                end           
            end        
        else
            currNode = 1;
            while currNode <= length(S{l})
                [decoded, llrs] = scFlip(softMess, mask, S{l}{currNode});
                check = mod(decoded * getGC(K, g), 2);
                if sum(check(K + 1:end)) == 0
                    return
                elseif sum(check(K + 1:end)) ~=0 && (l + 1) <= maxLevel

                    if isPruningTechnique
                        iMax = max(S{l}{currNode});
                        tmpMask = mask;
                        tmpMask(criticalSet(mask)) = 0;
                        metricLeft = mu - gamaLeft * sqrt(2 * mu);
                        N1 = length(find(tmpMask(iMax: end)));
                        %llrs(tmpMask(iMax: end) == 1) < metric(tmpMask(iMax: end) == 1);
                        N2 = length(find(llrs(tmpMask(iMax: end) == 1) < metricLeft(tmpMask(iMax: end) == 1)));
                        noChild = ((N2 / N1) >= w(l));
                        if ~noChild
                            CS = modifyCriticalSet(mask, max(S{l}{currNode}));
                            metricRight = mu + gamaRight * sqrt(2 * mu);
                            recordIdx = [];
                            for i = 1:length(CS)
                                notSelect = (llrs(CS(i)) > metricRight(CS(i)));
                                if notSelect
                                    recordIdx(end + 1) = i;
                                end
                            end
                            CS(recordIdx) = [];
                            M = abs(llrs(CS) ./ sqrt(mu(CS)));
                            [~, I] = sort(M);
                            CS = CS(I);
                            for u = CS
                                S{l + 1}(end + 1) = {S{l}{currNode}};
                                S{l + 1}{end}(end + 1) = u;
                            end
                        end
                    else
                        %iMax = max(S{l}{currNode});
                        CS = modifyCriticalSet(mask, max(S{l}{currNode}));
                        M = abs(llrs(CS) ./ sqrt(mu(CS)));
                        [~, I] = sort(M);
                        CS = CS(I);
                        for u = CS
                            S{l + 1}(end + 1) = {S{l}{currNode}};
                            S{l + 1}{end}(end + 1) = u;
                        end
                    end
                end
                currNode = currNode + 1;
            end
        end
        l = l + 1;
        %{
        if l == 2
            S
            pause()
        end
        %}
    end
end

function CS = modifyCriticalSet(mask, iMax)
    if iMax == 0
        CS = criticalSet(mask);
    else
        mask(1:iMax) = 0;
        CS = criticalSet(mask);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%CA-Hybrid Decoder%%%%%%%%%%%%%%%%%%%%%%%%%

function [decoded, iterations] = CA_HD(y, softMess, mask, g, parameter, msg)
    L = 1;%1;
    LMax = parameter;
    N = length(mask);
    crcWidth = length(g) - 1;
    K = length(find(mask)); % K = KI + crcWidth
    KI = K - crcWidth;
    squareRadius = inf;
    uOptimal = zeros(1, N);
    while L <= LMax
        %fprintf('L = %d, LMax = %d, AD-SCL\n', L, LMax)
        LLRs = zeros(L, 2 * N - 1);
        LLRs(1, N:end) = softMess;
        Bits = zeros(2, N - 1, L);
        decodeds = zeros(L, N);
        pathMetrics = zeros(1, L);
        pathNum = 1;
        iterations = 0;
        for i = 1:N
            for j = 1:pathNum
                LLRs(j, :) = updateLLRs(i, N, LLRs(j, :), Bits(:, :, j));
                iterations = iterations + 1;
            end

            if mask(i) == 1

                if 2 * pathNum <= L

                    for i1 = 1:pathNum
                        penalty = abs(LLRs(i1, 1));
                        if penalty == 0
                            disp('llr warning')
                            pause()
                        end
                        branchMetric0 = (0 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                        branchMetric1 = (1 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                        pathMetrics(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                        pathMetrics(i1) = pathMetrics(i1) + branchMetric0;
                        decodeds(i1 + pathNum, :) = decodeds(i1, :);
                        decodeds(i1 + pathNum, i) = 1;
                        LLRs(i1 + pathNum, :) = LLRs(i1, :);
                        Bits(:, :, i1 + pathNum) = Bits(:, :, i1);
                    end
                    pathNum = 2 * pathNum;

                else

                    PM = zeros(1, 2 * pathNum);
                    discardTag = ones(1, 2 * pathNum);
                    pathStates = zeros(1, pathNum); % 1 for both of paths that should be discarded; 2 for one of paths that should be discarded; 3 for both of paths that should be retained
                    for i1 = 1:pathNum
                        penalty = abs(LLRs(i1, 1));
                        if penalty == 0
                            disp('llr warning')
                            pause()
                        end
                        branchMetric0 = (0 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                        branchMetric1 = (1 ~= 1 / 2 * (1 - sign(LLRs(i1, 1)))) * penalty;
                        PM(i1) = pathMetrics(i1) + branchMetric0;
                        PM(i1 + pathNum) = pathMetrics(i1) + branchMetric1;
                    end

                    [~, indices] = sort(PM);
                    discardTag(indices(1:L)) = 0;
                    disIdxRecord = [];
                    for i1 = 1:pathNum
                        if discardTag(i1) == 1 && discardTag(i1 + pathNum) == 1
                            pathStates(i1) = 1;
                            disIdxRecord(end + 1) = i1;
                        elseif discardTag(i1) == 1 || discardTag(i1 + pathNum) == 1
                            pathStates(i1) = 2;
                        elseif discardTag(i1) == 0 && discardTag(i1 + pathNum) == 0
                            pathStates(i1) = 3;
                        else
                            disp('state warning')
                        end
                    end

                    for i1 = 1:length(pathStates)
                        if pathStates(i1) == 2
                            if discardTag(i1) == 1
                                pathMetrics(i1) = PM(i1 + pathNum);
                                decodeds(i1, i) = 1;
                            elseif discardTag(i1 + pathNum) == 1
                                pathMetrics(i1) = PM(i1);
                            end
                        elseif pathStates(i1) == 3
                            pathMetrics(i1) = PM(i1);
                            pathMetrics(disIdxRecord(1)) = PM(i1 + pathNum);
                            decodeds(disIdxRecord(1), :) = decodeds(i1, :);
                            decodeds(disIdxRecord(1), i) = 1;
                            LLRs(disIdxRecord(1), :) = LLRs(i1, :);
                            Bits(:, :, disIdxRecord(1)) = Bits(:, :, i1);
                            disIdxRecord(1) = [];
                        end
                    end
                end

            else

                for j = 1:pathNum
                    decodeds(j, i) = 0;
                    penalty = abs(LLRs(j, 1));
                    if penalty == 0
                        disp('llr warning1')
                        pause()
                    end
                    branchMetric0 = (decodeds(j, i) ~= 1 / 2 * (1 - sign(LLRs(j, 1)))) * penalty;
                    pathMetrics(j) = pathMetrics(j) + branchMetric0;
                end

            end
            for j = 1:pathNum
                Bits(:, :, j) = updateBits(i, N, Bits(:, :, j), decodeds(j, i));
            end
        end

        [~, I] = sort(pathMetrics);
        decodeds = decodeds(I, :);
        
        for l = 1:size(decodeds, 1)
            best = decodeds(l, :);
            best = best(mask == 1);
            check = mod(best * getGC(K, g), 2);
            if sum(check(K + 1:end)) == 0
                %decoded = best(1:KI);
                fprintf('AD-SCL succeeded in L = %d\n', L)
                decoded = best;
                return
            end
        end
        
        L = L * 2;
    end
    
    for l = 1:LMax
        u = decodeds(l, :);
        u = u(mask == 1);
        b = u(1:KI);
        %s = mod(b * getGC(KI, g), 2);
        %u(KI + 1, :) = s(KI + 1, :);
        xTilde = CRCPolarEncode(b, mask, g);
        eucDistance = 0;
        for i = 1:N
            eucDistance = eucDistance + (((1 - y(i)) / 2) - xTilde(i)) ^ 2;
        end
        if eucDistance < squareRadius
            squareRadius = eucDistance;
            uOptimal(mask == 1) = mod(b * getGC(KI, g), 2);
        end
    end

    fprintf('For CA-SCL, square radius = %f\n', squareRadius)
    correctRadius = 0;
    xTilde1 = CRCPolarEncode(msg, mask, g);
    for i = 1:N
        correctRadius = correctRadius + (((1 - y(i)) / 2) - xTilde1(i)) ^ 2;
    end
    if squareRadius > correctRadius
        squareRadius = correctRadius;
        uOptimal(mask == 1) = mod(msg * getGC(KI, g), 2);
    end
    %fprintf('ready to start CA-SD, square radius = %f\n', squareRadius)
    fprintf('ready to start CA-SD, correct square radius = %f\n', squareRadius)
    [decoded, iterations] = CA_SD(y, softMess, mask, g, squareRadius, uOptimal);
end

%%%%%%%%%%%%%%%%%%%%%%%%%PAC-Fano Decoder%%%%%%%%%%%%%%%%%%%%%%%%%

function decoded = PACFano(softMess, mask, g, parameter, dsnrdB, snrType)
    N = length(mask);
    K = length(find(mask)); % K = KI + crcWidth
    LLRs = zeros(1, 2 * N - 1);
    LLRs(N:end) = softMess;
    Bits = zeros(2, N - 1);
    decoded = zeros(1, N);
    polarDecoded = zeros(1, N);
    pathMetrics = zeros(1, N);
    followOtherBranch = zeros(1, N);
    %pe = PE(N, K, dsnrdB, snrType);
    I = symmetricCapacity(N, K, dsnrdB, snrType);
    i = 0;
    T = 0;
    m = length(g) - 1;
    delta = parameter;
    currState = zeros(1, m);
    state = cell(1, N);
    flag = true;
    while flag
        pointer = i + 1;
        if mask(pointer) == 1
            LLRs = updateLLRs(pointer, N, LLRs, Bits);
            u0 = conv1Bit(0, currState, g);
            nextState0 = getNextState(0, currState, m);
            u1 = conv1Bit(1, currState, g);
            nextState1 = getNextState(1, currState, m);
            %branchMetric0 = log((exp(LLRs(1)) / (exp(LLRs(1)) + 1)) / (1 - pe(pointer)));
            %branchMetric1 = log((1 / (exp(LLRs(1)) + 1)) / (1 - pe(pointer)));
            branchMetric0 = 1 - log2(1 + exp(LLRs(1)) ^ (-(1 - 2 * u0))) - I(pointer);
            branchMetric1 = 1 - log2(1 + exp(LLRs(1)) ^ (-(1 - 2 * u1))) - I(pointer);
            if pointer == 1
                pathMetric0 = 0 + branchMetric0;
                pathMetric1 = 0 + branchMetric1;
            else
                pathMetric0 = pathMetrics(pointer - 1) + branchMetric0;
                pathMetric1 = pathMetrics(pointer - 1) + branchMetric1;
            end
            if followOtherBranch(pointer) == 0
                pathMetrics(pointer) = max([pathMetric0, pathMetric1]);
            else
                pathMetrics(pointer) = min([pathMetric0, pathMetric1]);
            end
        else
            LLRs = updateLLRs(pointer, N, LLRs, Bits);
            u0 = conv1Bit(0, currState, g);
            nextState0 = getNextState(0, currState, m);
            %branchMetric0 = log((exp(LLRs(1)) / (exp(LLRs(1)) + 1)) / (1 - pe(pointer)));
            branchMetric0 = 1 - log2(1 + exp(LLRs(1)) ^ (-(1 - 2 * u0))) - I(pointer);
            if pointer == 1
                pathMetric0 = 0 + branchMetric0;
            else
                pathMetric0 = pathMetrics(pointer - 1) + branchMetric0;
            end
            pathMetrics(pointer) = pathMetric0;
        end

        if pathMetrics(pointer) >= T
            i = i + 1;
            if i == 0
                muPre = -inf;
            elseif i == 1
                muPre = 0;
            else
                muPre = pathMetrics(i - 1);
            end
            if muPre < (T + delta)
                while (T + delta) <= pathMetrics(i)
                    T = T + delta;
                end
            end
            if pathMetrics(i) == pathMetric0
                decoded(i) = 0;
                polarDecoded(i) = u0;
                currState = nextState0;
                state{i} = nextState0;
            else
                decoded(i) = 1;
                polarDecoded(i) = u1;
                currState = nextState1;
                state{i} = nextState1;
            end
            Bits = updateBits(i, N, Bits, polarDecoded(i));
            if i == N
                decoded = decoded(mask == 1);
                %flag = false
                return
            end
            followOtherBranch(i + 1) = 0;
        else
            j = i;
            [i, T, followOtherBranch] = moveBackPAC(i, T, delta, mask, pathMetrics, followOtherBranch);
            if i == 0
                currState = zeros(1, m);
            elseif j ~= i
                currState = state{i};
            end
            if j ~= i
                for ii = 1:i
                    LLRs = updateLLRs(ii, N, LLRs, Bits);
                    Bits = updateBits(ii, N, Bits, polarDecoded(ii));
                end
            end
        end

    end
end

function [i, T, followOtherBranch] = moveBackPAC(i, T, delta, mask, pathMetrics, followOtherBranch)
    while true
        if i == 0
            muPre = -inf;
        elseif i == 1
            muPre = 0;
        else
            muPre = pathMetrics(i - 1);
        end
        if muPre < T
            T = T - delta;
            followOtherBranch(i + 1) = 0;
            return
        else
            if (followOtherBranch(i) + 1) == 2 || mask(i) == 0
                i = i - 1;
            else
                followOtherBranch(i) = followOtherBranch(i) + 1;
                i = i - 1;
                return
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%Auxiliary Functions%%%%%%%%%%%%%%%%%%%%%%%%%

function Fn = getFn(F, n)
    if n == 1
        Fn = F;
    else
        Fn = kron(F, getFn(F, n - 1));
    end
end

function GC = getGC(KI, g)
    KC = length(g) - 1;
    K = KC + KI;
    identityMatrix = eye(KI);
    matrixRight = zeros(KI, KC);
    for i = 1:KI
        x = zeros(K - i + 1);
        x(1) = 1;
        matrixRight(i, :) = polynomialMod(x, g);
    end
    GC = [identityMatrix, matrixRight];
end

function remainder = polynomialMod(a, b)
    lenA = length(a);
    lenB = length(b);
    lenC = lenA - lenB + 1;
    if lenC < 0
        disp('lenC negative warning')
    else
        c = zeros(lenC);
        tmp = a(1:lenB);
        for i = 1:lenC
            if tmp(1) == 0
                c(i) = 0;
            else
                c(i) = 1;
            end
            tmp = mod(tmp + b * c(i), 2);
            if i ~= lenC
                tmp(end + 1) = a(lenB + i);
                tmp = tmp(2:end);
            end
        end
    end
    remainder = tmp(2:end);
end

function A = rowEchelonForm(A) % design for modulo 2
    [rowNum, columnNum] = size(A);
    indices = [];
    for column = 1:columnNum
        lock = false;
        for row = 1:rowNum
            if ~sum(A(row, 1:(column - 1))) && A(row, column) ~= 0 && lock == false
                tmp = A(row, :);
                indices(end + 1) = row;
                lock = true;
            elseif ~sum(A(row, 1:(column - 1))) && A(row, column) ~= 0 && lock == true               
                common = lcm(A(row, column), tmp(column));
                A(row, :) = A(row, :) * common / A(row, column);
                A(row, :) = A(row, :) + tmp * (-1) * common / tmp(column);
                A(row, :) = mod(A(row, :), 2);
            else
                continue
            end
        end
    end
    if length(indices) ~= rowNum
        for row = 1:rowNum
            if ~ismember(row, indices)
                indices(end + 1) = row;
            end
        end
    end
    %{
    E = zeros(size(A));
    for row = 1:rowNum
        E(row, :) = A(indices(row), :);
    end
    %}
    A = A(indices, :);
end

function LLRs = updateLLRs(i, N, LLRs, Bits)
    
    f = @(upperDecision, upperLLR, lowerLLR) (1 - 2 * upperDecision) * upperLLR + lowerLLR;
    g = @(upperLLR, lowerLLR) sign(upperLLR) * sign(lowerLLR) * min(abs(upperLLR), abs(lowerLLR));
    
    n = log2(N);
    position = bin2dec(fliplr(dec2bin(i - 1, n)));
    if position == 0
        nextLevel = n;
    else
        indices1 = dec2bin(position, n);
        for i = 1:length(indices1)
            if indices1(i) == '1'
                lastLevel = i;
                break
            end
        end
        %indices1 = find(dec2bin(position, n));
        %lastLevel = indices1(1);
        start = 2 ^ (lastLevel - 1);
        final = 2 ^ (lastLevel) - 1;
        for i = start:final
            exp1 = final + (i - start);
            LLRs(i) = f(Bits(1, i), LLRs(exp1 + 1), LLRs(exp1 + 1 + 2 ^ (lastLevel - 1)));
        end
        nextLevel = lastLevel - 1;
    end

    for lev = nextLevel:-1:1
        start = 2 ^ (lev - 1);
        final = 2 ^ (lev) - 1;
        for i = start:final
            exp1 = final + (i - start);
            LLRs(i) = g(LLRs(exp1 + 1), LLRs(exp1 + 1 + 2 ^ (lev - 1)));
        end
    end
end

function Bits = updateBits(i, N, Bits, latestBit)
    n = log2(N);
    position = bin2dec(fliplr(dec2bin(i - 1, n)));
    if position == N - 1
        return
    elseif position < N / 2
        Bits(1, 1) = latestBit;
    else
        indices0 = dec2bin(position, n);
        for i = 1:length(indices0)
            if indices0(i) == '0'
                lastLevel = i;
                break
            end
        end
        Bits(2, 1) = latestBit;
        %indices0 = find(dec2bin(position, n));
        %lastLevel = indices0(1);
        for lev = 1:(lastLevel - 2)
            start = 2 ^ (lev - 1);
            final = 2 ^ (lev) - 1;
            for i = start:final
                exp1 = final + (i - start);
                Bits(2, exp1 + 1) = mod((Bits(1, i) + Bits(2, i)), 2);
                Bits(2, exp1 + 1 + 2 ^ (lev - 1)) = Bits(2, i);
            end
        end
        lev = lastLevel - 1;
        start = 2 ^ (lev - 1);
        final = 2 ^ (lev) - 1;
        for i = start:final
            exp1 = final + (i - start);
            Bits(1, exp1 + 1) = mod((Bits(1, i) + Bits(2, i)), 2);
            Bits(1, exp1 + 1 + 2 ^ (lev - 1)) = Bits(2, i);
        end
    end
end

function I = symmetricCapacity(N, K, dsnrdB, snrType)
    if strcmpi(snrType, 'snr')
        sigmaSquare = 1 / (2 * power(10, dsnrdB / 10));
    elseif strcmpi(snrType, 'snrb')
        R = K / N;
        sigmaSquare = 1 / (2 * R * power(10, dsnrdB / 10));
    end

    I = JFunction(2 / sigmaSquare);
    while length(I) ~= N
        tmp = [];
        for i = 1:length(I)
            tmp(end + 1) = 1 - JFunction(sqrt(2) * JFunctionInv(1 - I(i)));
            tmp(end + 1) = JFunction(sqrt(2) * JFunctionInv(I(i)));
        end
        I = tmp;
    end
end

function pe = PE(N, K, dsnrdB, snrType)
    mllr = GA(N, K, dsnrdB, snrType);
    %pe = zeros(1, N);
    pe = 1 / 2 * erfc(1 / sqrt(2) * sqrt(mllr(1:N) / 2));
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

function y = JFunction(t)
    y = (1 - 2 ^ (-0.3073 * (t ^ (2 * 0.8935)))) ^ 1.1064;
end

function y = JFunctionInv(t)
    y = (-1 / 0.3073 * log2(1 - (t ^ (1 / 1.1064)))) ^ (1 / (2 * 0.8935));
end

function CS = criticalSet(mask)
    N = length(mask);
    n = log2(N);
    tree = zeros(n + 1, N);
    tree(n + 1, mask == 1) = 1;
    CS = [];
    for i = (n + 1):-1:1
        if i == (n + 1)
            continue
        else
            for j = 1:2 ^ (i - 1)
                if tree(i + 1, 2 * j - 1) == 1 && tree(i + 1, 2 * j) == 1
                    tree(i, j) = 1;
                %else
                 %   continue
                end
            end
        end
    end
   
    for i = 1:(n + 1)
        if tree(1, 1) == 1
            CS(end + 1) = 1;
            return
        else
            for j = 1:2 ^ (i - 1)
                if tree(i, j) == 1 && tree(i - 1, ceil(j / 2)) ~= 1
                    tmp = j;
                    for k = i:n
                        tmp = tmp * 2 - 1;
                    end
                    CS(end + 1) = tmp;
                end
            end
        end
    end
    CS = sort(CS);
end

%%%%%%%%%%%%%%%%%%%%%%%%%Polar Codes Encoder%%%%%%%%%%%%%%%%%%%%%%%%%

function x = polarEncode(msg, mask)
    N = length(mask);
    n = log2(N);
    u = zeros(1, N);
    u(mask == 1) = msg;
    x = u;
    for i = 1:n
        if i == 1
           x(1:2:N) = mod(x(1:2:N) + x(2:2:N), 2);
        elseif i == n
            x(1:N / 2) = mod(x(1:(N / 2)) + x((N / 2 + 1):N), 2);
        else
            encStep = power(2, i - 1);
        for j = 1:encStep
            x(j:(2 * encStep):N) = mod(x(j:(2 * encStep):N) + x((j + power(2, i - 1)):(2 * encStep):N), 2);
        end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%CRC-Polar Codes Encoder%%%%%%%%%%%%%%%%%%%%%%%%%
function x = CRCPolarEncode(msg, mask, g)
    s = mod(msg * getGC(length(msg), g), 2);
    x = polarEncode(s, mask);
end

%%%%%%%%%%%%%%%%%%%%%%%%%Convolutional Codes Encoder%%%%%%%%%%%%%%%%%%%%%%%%%

function gBit = conv1Bit(inBit, currState, g)
    gLen = length(g);
    gBit = inBit * g(1);
    for i = 2:gLen
        if g(i) == 1
            gBit = bitxor(gBit, currState(i - 1));
        end
    end
end

function nextState = getNextState(inBit, currState, m)
    if inBit == 0
        nextState = 0;
        nextState(2:m) = currState(1:m - 1);
    else
        nextState = 1;
        nextState(2:m) = currState(1:m - 1);
    end
end