function x = encoder(msg, mask, crcPolynomial, convG, encoderType)
    if strcmpi(encoderType, 'polar')
        x = polarEncode(msg, mask);
    elseif strcmpi(encoderType, 'crcpolar')
        x = CRCPolarEncode(msg, mask, crcPolynomial);
    elseif strcmpi(encoderType, 'pac')
        x = PACEncode(msg, mask, convG);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%Polar Encoder%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%CRC-Polar Encoder%%%%%%%%%%%%%%%%%%%%%%%%%
function x = CRCPolarEncode(msg, mask, g)
    s = mod(msg * getGC(length(msg), g), 2);
    x = polarEncode(s, mask);
end

%%%%%%%%%%%%%%%%%%%%%%%%%PAC Encoder%%%%%%%%%%%%%%%%%%%%%%%%%

function x = PACEncode(msg, mask, g)
    N = length(mask);
    m = length(g) - 1;
    v = zeros(1, N);
    v(mask == 1) = msg;
    u = convEncode(v, g, m);
    x = mulMatrix(u);
end

%%%%%%%%%%%%%%%%%%%%%%%%%Auxiliary Functions%%%%%%%%%%%%%%%%%%%%%%%%%

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

function x = mulMatrix(u)
    N = length(u);
    n = log2(N);
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

function convCode = convEncode(inCode, g, m)
    currState = zeros(1, m);
    inCodeLen = length(inCode);
    convCode = zeros(1, inCodeLen);
    for i = 1:inCodeLen
        inBit = inCode(i);
        output = conv1Bit(inBit, currState, g);
        currState = getNextState(inBit, currState, m);
        convCode(i) = output;
    end
end
