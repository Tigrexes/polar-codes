clear;clc; close all

isParallel = false;
isCRC = true;
if isParallel
    CoreNum = 6;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(CoreNum);
    else
        disp('matlab pool already started');
    end
end

N = 2 ^ 7;
R = 1 / 2;
KI = round(N * R);
parameter = 1024; % if list size do not equal to the power of 2, it may go wrong
parameter2 = 1;
parameter3 = 16;
convG = [1 0 1 1 0 1 1];
if isCRC
    %crcPolynomial = [1 1 1 0 0 0 1 0 0 0 0 1]; % for CA-SD, N = 64, R = 1 / 4; 1 / 2
    %crcPolynomial = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1]; % for progressive bit-flipping
    %crcPolynomial = [1 0 0 0 0 0 1 1 1]; % for short CA-SCL polar codes, CRC width = 8
    %crcPolynomial = [1 1 0 0 0 0 1]; % for short CA-SCL polar codes, CRC width = 6
    %crcPolynomial = [1 0 0 0 1 1 1 1 1 0 0 0 1 0 1 0 1]; % for CA-HD, N = 128, R = 1 / 3 0001 0001 1111 0001 0101
    crcPolynomial = [1 0 1 0 0 0 1 0 0 0 1 0 1 1 1 0 0 1 0 1]; % for CA-HD, N = 128, R = 1 / 2
    %crcPolynomial = [1 0 0 0 1 0 1 1 1 1 0 1 1 0 1 1 1]; % for CA-HD, N = 128, R = 2 / 3 0001 0001 0111 1011 0111
    %crcPolynomial = [1 1 1 1 0 1 0 1 0 1 0 1 0 1]; % for CA-HD, N = 64, R = 1 / 2 0011 1101 0101 0101
else
    crcPolynomial = NaN;
end
K = KI + length(crcPolynomial) - 1;
K2 = KI;
%K3 = KI;
%dsnrdB = 1.0;
SNRdB = -2.0:0.5:1.5;
snrType = 'SNR'; % SNR for EsN0; SNRb for EbN0
encoderType = 'crcpolar'; % CRC-Polar/polar for CRC-Polar/polar codes encoder; pac for polarization-adjusted convolutional (PAC) codes encoder;
encoderType2 = 'polar';
%encoderType3 = 'polar';
decoderType = 'cahd';
decoderType2 = 'scl';
%decoderType3 = 'scl';
construction = 'dega';
maxIterations = 1e9;
errCnt = 100;
curvesNum = 1;
FER = zeros(curvesNum, length(SNRdB));
Iterations = zeros(curvesNum, length(SNRdB));
clustNum = ones(1, length(SNRdB)) * 1000;
SimulatorParameter = struct('encoderType', {'crcpolar', 'crcpolar', 'polar'}, ...
                            'decoderType', {'casd', 'scl', 'scl'}, ...
                            'construction', {'dega', 'dega', 'dega'}, ...
                            'parameter', {NaN, 32, 32}, ...
                            'isCRC', {true, true, false}, ...
                            'crcPolynomial', {[1 1 1 0 0 0 1 0 0 0 0 1], [1 1 0 0 0 0 1], NaN});

for i = 1:length(SNRdB)
    snrdB = SNRdB(i);
    if strcmpi(snrType, 'snr')
        sigmaSquare = 1 / (2 * power(10, snrdB / 10));
    elseif strcmpi(snrType, 'snrb')
        R = K / N;
        sigmaSquare = 1 / (2 * R * power(10, snrdB / 10));
    end
    fprintf('Starting! SNR = %.2fdB, snr type = %s\n', snrdB, snrType);

    dsnrdB = snrdB;
    mask = rateProfile(N, K, dsnrdB, snrType, construction);
    mask2 = rateProfile(N, K2, dsnrdB, snrType, construction);
    %mask3 = rateProfile(N, K3, dsnrdB, snrType, construction);
    
    frameError = 0;
    AVN = 0;
    frameError2 = 0;
    AVN2 = 0;
    frameError3 = 0;
    AVN3 = 0;
    if isParallel

        for j = 1:maxIterations / clustNum(i)

            parfor k = 1:clustNum(i)
                msg = randi([0 1], 1, KI);
                x = encoder(msg, mask, crcPolynomial, convG, encoderType);
                %x = encoder(msg, mask, [1 1 1 0 0 0 1 0 0 0 0 1], convG, encoderType);
                %x2 = encoder(msg, mask2, crcPolynomial, convG, encoderType2);
                %x3 = encoder(msg, mask3, [1 1 0 0 0 0 1], convG, encoderType3);
                modulated = 1 - 2 * x;
                %modulated2 = 1 - 2 * x2;
                %modulated3 = 1 - 2 * x3;
                noise = sqrt(sigmaSquare) * randn(1, N);
                %y = modulated + sqrt(sigmaSquare) * randn(1, N);
                y = modulated + noise;
                %y2 = modulated2 + noise;
                %y3 = modulated3 + noise;
                llr = 2 / sigmaSquare * y;
                %llr2 = 2 / sigmaSquare * y2;
                %llr3 = 2 / sigmaSquare * y3;
                [decoded, iterations] = decoder(y, llr, mask, crcPolynomial, convG, parameter, dsnrdB, isCRC, snrType, decoderType, msg);
                [decoded2, iterations2] = decoder(y, llr, mask, crcPolynomial, convG, parameter2, dsnrdB, isCRC, snrType, 'polarfano', msg);
                [decoded3, iterations3] = decoder(y, llr, mask, crcPolynomial, convG, parameter3, dsnrdB, isCRC, snrType, 'scl', msg);
                %if isCRC
                 %   decoded = decoded(1:KI);
                %end
                decoded = decoded(1:KI);
                decoded2 = decoded2(1:KI);
                decoded3 = decoded3(1:KI);

                frameError = frameError + any(decoded~=msg);
                AVN = AVN + iterations;
                frameError2 = frameError2 + any(decoded2~=msg);
                AVN2 = AVN2 + iterations2;
                frameError3 = frameError3 + any(decoded3~=msg);
                AVN3 = AVN3 + iterations3;
            end
            if min([frameError, frameError2]) >= errCnt && (j * clustNum(i)) >= 2 * 1e4
                fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, decoder type = Sphere\n', snrdB, snrType, frameError / (j * clustNum(i)))
                fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, decoder type = Fano\n', snrdB, snrType, frameError2 / (j * clustNum(i)))
                fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, decoder type = SCL\n\n', snrdB, snrType, frameError3 / (j * clustNum(i)))
                break;
            end
%             if frameError >= errCnt && (j * clustNum(i)) >= 2 * 1e4
%                 fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, decoder type = CA-SD\n', snrdB, snrType, frameError / (j * clustNum(i)))
%                 break;
%             end
            fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, decoder type = Sphere\n', snrdB, snrType, j * clustNum(i), frameError, frameError / (j * clustNum(i)))
            fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, decoder type = Fano\n', snrdB, snrType, j * clustNum(i), frameError2, frameError2 / (j * clustNum(i)))
            fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, decoder type = SCL\n', snrdB, snrType, j * clustNum(i), frameError3, frameError3 / (j * clustNum(i)))
        end
        temp = j * clustNum(i);
        FER(1, i) = frameError / temp;
        Iterations(1, i) = AVN / temp;
        FER(2, i) = frameError2 / temp;
        Iterations(2, i) = AVN2 / temp;
        FER(3, i) = frameError3 / temp;
        Iterations(3, i) = AVN3 / temp;
    else

        for j = 1:maxIterations

            
            msg = randi([0 1], 1, KI);
            %noise = sqrt(sigmaSquare) * randn(1, N);
            x = encoder(msg, mask, crcPolynomial, convG, encoderType);
            %x2 = encoder(msg, mask2, [1 1 0 0 0 0 1], convG, encoderType);
            modulated = 1 - 2 * x;
            %modulated2 = 1 - 2 * x2;
            y = modulated + sqrt(sigmaSquare) * randn(1, N);
            %y = modulated + noise;
            %y2 = modulated2 + noise;
            llr = 2 / sigmaSquare * y;
            %llr2 = 2 / sigmaSquare * y2;
            
            eucDistance = 0;
            for h = 1:N
                eucDistance = eucDistance + (((1 - y(h)) / 2) - x(h)) ^ 2;
            end
            fprintf('Correct Radius = %.2f\n', eucDistance)
            
            [decoded, iterations] = decoder(y, llr, mask, crcPolynomial, convG, parameter, dsnrdB, isCRC, snrType, decoderType, msg);
            %[decoded2, iterations2] = decoder(y2, llr2, mask2, [1 1 0 0 0 0 1], convG, parameter2, dsnrdB, isCRC, snrType, decoderType2, msg);
            %{
            if isCRC
                decoded = decoded(1:KI);
                decoded2 = decoded2(1:KI);
            end
            %}
            decoded = decoded(1:KI);
            %decoded2 = decoded2(1:KI);
            
            eucDistance = 0;
            eucDistance2 = 0;
            error = encoder(decoded, mask, crcPolynomial, convG, encoderType);
            for h = 1:N
                eucDistance = eucDistance + (((1 - y(h)) / 2) - x(h)) ^ 2;
                eucDistance2 = eucDistance2 + (((1 - y(h)) / 2) - error(h)) ^ 2;
            end
            %disp([eucDistance, eucDistance2])
            if eucDistance < eucDistance2
                disp('error warning12')
                disp(eucDistance)
                disp(eucDistance2)
                pause()
            end
            

            frameError = frameError + any(decoded~=msg);
            AVN = AVN + iterations;
            %frameError2 = frameError2 + any(decoded2~=msg);
            %AVN2 = AVN2 + iterations2;
%             if min(frameError, frameError2) >= errCnt && j >= 2 * 1e4
%                 fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, \n', snrdB, snrType, frameError / j)
%                 fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, \n\n', snrdB, snrType, frameError2 / j)
%                 break;
%             end
            if frameError >= errCnt && j >= 2 * 1e4
                fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, \n\n', snrdB, snrType, frameError / j)
                break;
            end
            if any(decoded~=msg)
               fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, AVN = %.2e, Sphere\n', snrdB, snrType, j, frameError, frameError / j, AVN / j)
            end
%             if any(decoded2~=msg)
%                 fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, AVN = %.2e, CA-SCL\n', snrdB, snrType, j, frameError2, frameError2 / j, AVN2 / j)
%             end
        end
        FER(1, i) = frameError / j;
        Iterations(1, i) = AVN / j;
        FER(2, i) = frameError2 / j;
        Iterations(2, i) = AVN2 / j;
    end
end

%{
N = 2 ^ 6;
R = 1 / 4;
KI = round(N * R);
parameter = NaN; % if list size do not equal to the power of 2, it may go wrong
parameter2 = 32;
parameter3 = 32;
convG = [1 0 1 1 0 1 1];
if isCRC
    crcPolynomial = [1 1 1 0 0 0 1 0 0 0 0 1]; % for CA-SD, N = 64, R = 1 / 4; 1 / 2
    %crcPolynomial = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1]; % for progressive bit-flipping
    %crcPolynomial = [1 0 0 0 0 0 1 1 1]; % for short CA-SCL polar codes, CRC width = 8
    %crcPolynomial = [1 1 0 0 0 0 1]; % for short CA-SCL polar codes, CRC width = 6
    %crcPolynomial = [1 0 0 0 1 1 1 1 1 0 0 0 1 0 1 0 1]; % for CA-HD, N = 128, R = 1 / 3 0001 0001 1111 0001 0101
    %crcPolynomial = [1 0 1 0 0 0 1 0 0 0 1 0 1 1 1 0 0 1 0 1]; % for CA-HD, N = 128, R = 1 / 2
    %crcPolynomial = [1 0 0 0 1 0 1 1 1 1 0 1 1 0 1 1 1]; % for CA-HD, N = 128, R = 2 / 3 0001 0001 0111 1011 0111
    %crcPolynomial = [1 1 1 1 0 1 0 1 0 1 0 1 0 1]; % for CA-HD, N = 64, R = 1 / 2 0011 1101 0101 0101
else
    crcPolynomial = NaN;
end
K = KI + length([1 1 1 0 0 0 1 0 0 0 0 1]) - 1;
K2 = KI + length([1 1 0 0 0 0 1]) - 1;
K3 = KI;
%dsnrdB = 1.0;
SNRdB = -5.0:0.5:-1.5;
snrType = 'SNR'; % SNR for EsN0; SNRb for EbN0
encoderType = 'crcpolar'; % CRC-Polar/polar for CRC-Polar/polar codes encoder; pac for polarization-adjusted convolutional (PAC) codes encoder;
encoderType2 = 'crcpolar';
encoderType3 = 'polar';
decoderType = 'casd';
decoderType2 = 'scl';
decoderType3 = 'scl';
construction = 'dega';
maxIterations = 1e9;
errCnt = 200;
curvesNum = 3;
FER = zeros(curvesNum, length(SNRdB));
Iterations = zeros(curvesNum, length(SNRdB));
clustNum = ones(1, length(SNRdB)) * 1000;
SimulatorParameter = struct('encoderType', {'crcpolar', 'crcpolar', 'polar'}, ...
                            'decoderType', {'casd', 'scl', 'scl'}, ...
                            'construction', {'dega', 'dega', 'dega'}, ...
                            'parameter', {NaN, 32, 32}, ...
                            'isCRC', {true, true, false}, ...
                            'crcPolynomial', {[1 1 1 0 0 0 1 0 0 0 0 1], [1 1 0 0 0 0 1], NaN});

%mask = rateProfile(N, K, dsnrdB, snrType, construction);

for i = 1:length(SNRdB)
    snrdB = SNRdB(i);
    if strcmpi(snrType, 'snr')
        sigmaSquare = 1 / (2 * power(10, snrdB / 10));
    elseif strcmpi(snrType, 'snrb')
        R = K / N;
        sigmaSquare = 1 / (2 * R * power(10, snrdB / 10));
    end
    fprintf('Starting! SNR = %.2fdB, snr type = %s\n', snrdB, snrType);

    dsnrdB = snrdB;
    mask = rateProfile(N, K, dsnrdB, snrType, construction);
    mask2 = rateProfile(N, K2, dsnrdB, snrType, construction);
    mask3 = rateProfile(N, K3, dsnrdB, snrType, construction);
    
    %for ii = 1:curvesNum
     %   K = KI + length(SimulatorParameter(ii).crcPolynomial) - 1;
      %  SimulatorParameter(ii).mask = rateProfile(N, K, dsnrdB, snrType, SimulatorParameter(ii).construction);
    %end
    
    %frameError = zeros(1, curvesNum);
    %AVN = zeros(1, curvesNum);
    frameError = 0;
    AVN = 0;
    frameError2 = 0;
    AVN2 = 0;
    frameError3 = 0;
    AVN3 = 0;
    if isParallel

        for j = 1:maxIterations / clustNum(i)

            parfor k = 1:clustNum(i)
                msg = randi([0 1], 1, KI);
                x = encoder(msg, mask, [1 1 1 0 0 0 1 0 0 0 0 1], convG, encoderType);
                x2 = encoder(msg, mask2, [1 1 0 0 0 0 1], convG, encoderType2);
                x3 = encoder(msg, mask3, [1 1 0 0 0 0 1], convG, encoderType3);
                modulated = 1 - 2 * x;
                modulated2 = 1 - 2 * x2;
                modulated3 = 1 - 2 * x3;
                noise = sqrt(sigmaSquare) * randn(1, N);
                y = modulated + noise;
                y2 = modulated2 + noise;
                y3 = modulated3 + noise;
                llr = 2 / sigmaSquare * y;
                llr2 = 2 / sigmaSquare * y2;
                llr3 = 2 / sigmaSquare * y3;
                [decoded, iterations] = decoder(y, llr, mask, [1 1 1 0 0 0 1 0 0 0 0 1], convG, parameter, dsnrdB, true, snrType, decoderType, msg);
                [decoded2, iterations2] = decoder(y2, llr2, mask2, [1 1 0 0 0 0 1], convG, parameter2, dsnrdB, true, snrType, decoderType2, msg);
                [decoded3, iterations3] = decoder(y3, llr3, mask3, [1 1 0 0 0 0 1], convG, parameter3, dsnrdB, false, snrType, decoderType3, msg);
                %if isCRC
                 %   decoded = decoded(1:KI);
                %end
                decoded = decoded(1:KI);
                decoded2 = decoded2(1:KI);

                frameError = frameError + any(decoded~=msg);
                AVN = AVN + iterations;
                frameError2 = frameError2 + any(decoded2~=msg);
                AVN2 = AVN2 + iterations2;
                frameError3 = frameError3 + any(decoded3~=msg);
                AVN3 = AVN3 + iterations3;
            end
            if min([frameError, frameError2, frameError3]) >= errCnt && (j * clustNum(i)) >= 2 * 1e4
                fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, decoder type = CA-SD\n', snrdB, snrType, frameError / (j * clustNum(i)))
                fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, fer = %.2e, decoder type = CA-SCL\n', snrdB, snrType, frameError2 / (j * clustNum(i)))
                fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, fer = %.2e, decoder type = SCL\n\n', snrdB, snrType, frameError3 / (j * clustNum(i)))
                break;
            end
            fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, decoder type = CA-SD\n', snrdB, snrType, j * clustNum(i), frameError, frameError / (j * clustNum(i)))
            fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, decoder type = CA-SCL\n', snrdB, snrType, j * clustNum(i), frameError2, frameError2 / (j * clustNum(i)))
            fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, decoder type = SCL\n', snrdB, snrType, j * clustNum(i), frameError3, frameError3 / (j * clustNum(i)))
        end
        temp = j * clustNum(i);
        FER(1, i) = frameError / temp;
        Iterations(1, i) = AVN / temp;
        FER(2, i) = frameError2 / temp;
        Iterations(2, i) = AVN2 / temp;
        FER(3, i) = frameError3 / temp;
        Iterations(3, i) = AVN3 / temp;
    else

        for j = 1:maxIterations

            
            msg = randi([0 1], 1, KI);
            noise = sqrt(sigmaSquare) * randn(1, N);
            for ii = 1:curvesNum
                x{ii} = encoder(msg, SimulatorParameter(ii).mask, ...
                                SimulatorParameter(ii).crcPolynomial, convG, ...
                                SimulatorParameter(ii).encoderType);
                modulated{ii} = 1 - 2 * x{ii};
                y{ii} = modulated{ii} + noise;
                llr{ii} = 2 / sigmaSquare * y{ii};
                [decoded{ii}, iterations{ii}] = decoder(y{ii}, llr{ii}, SimulatorParameter(ii).mask, ...
                    SimulatorParameter(ii).crcPolynomial, convG, SimulatorParameter(ii).parameter, ...
                    dsnrdB, SimulatorParameter(ii).isCRC, snrType, SimulatorParameter(ii).decoderType, msg);
                if SimulatorParameter(ii).isCRC
                    decoded{ii} = decoded{ii}(1:KI);
                end
                frameError(ii) = frameError(ii) + any(decoded{ii}~=msg);
                AVN(ii) = AVN(ii) + iterations{ii};
                if any(decoded{ii}~=msg)
                    fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, decoder type = %s\n', ...
                            snrdB, snrType, j, frameError(ii), frameError(ii) / j, upper(SimulatorParameter(ii).decoderType))
                end
            end
            %x = encoder(msg, mask, crcPolynomial, convG, encoderType);
            %modulated = 1 - 2 * x;
            %y = modulated + sqrt(sigmaSquare) * randn(1, N);
            %llr = 2 / sigmaSquare * y;
            %{
            eucDistance = 0;
            for h = 1:N
                eucDistance = eucDistance + (((1 - y(h)) / 2) - x(h)) ^ 2;
            end
            eucDistance;
            %}
            %[decoded, iterations] = decoder(y, llr, mask, crcPolynomial, convG, parameter, dsnrdB, isCRC, snrType, decoderType, msg);
            %{
            if isCRC
                decoded = decoded(1:KI);
                %decoded2 = decoded2(1:KI);
            end
            %}
            %{
            eucDistance = 0;
            eucDistance2 = 0;
            error = encoder(decoded, mask, crcPolynomial, convG, encoderType);
            for h = 1:N
                eucDistance = eucDistance + (((1 - y(h)) / 2) - x(h)) ^ 2;
                eucDistance2 = eucDistance2 + (((1 - y(h)) / 2) - error(h)) ^ 2;
            end
            %disp([eucDistance, eucDistance2])
            if eucDistance < eucDistance2
                disp('error warning12')
                disp(eucDistance)
                disp(eucDistance2)
                pause()
            end
            %}

            %frameError = frameError + any(decoded~=msg);  
            %frameError2 = frameError2 + any(decoded2~=msg);
            if min(frameError) >= errCnt && j >= 2 * 1e4
                %fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, \n\n', snrdB, snrType, frameError / j)
                for ii = 1:curvesNum
                    fprintf('@ SNR = %.2fdB, snr type = %s fer = %.2e, decoder type = %s\n', snrdB, snrType, ...
                            frameError{ii} / j, upper(SimulatorParameter(ii).decoderType))
                end
                break;
            end
            %if any(decoded~=msg)
             %  fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e\n', snrdB, snrType, j, frameError, frameError / j)
               %fprintf('AVN = %.2e\n', AVN / j)
            %end
            %if any(decoded2~=msg)
             %  fprintf('SNR = %.2fdB, snr type = %s, t = %d, error %d, fer = %.2e, pscl\n', snrdB, snrType, j, frameError2, frameError2 / j)
              % fprintf('AVN = %.2e\n', AVN2 / j)
            %end
        end
        temp = j;
        for ii = 1:curvesNum
            FER(ii, i) = frameError(ii) / temp;
            Iterations(ii, i) = AVN(ii) / temp;
        end
    end
end
%}
%%
%{
semilogy(SNRdB, FER(1, :), '-o', 'LineWidth', 1.5, 'Color', '#228B22', 'DisplayName', 'Sphere')
hold on
semilogy(SNRdB, FER(2, :), '-s', 'LineWidth', 1.5, 'Color', '#FF0000', 'DisplayName', 'Fano (\Delta = 1)')
semilogy(SNRdB, FER(3, :), '-v', 'LineWidth', 1.5, 'Color', '#191970', 'DisplayName', 'List (L = 16)')
hold off
xlabel('E_b/N_0(dB)')
ylabel('BLER')
ylim([1e-3 1e-0])
%xlim([-5 0.5])
xticks(1:0.5:3.5)
legend()
grid on
%}
%{
bar(SNRdB, Iterations', 0.75, 'DisplayName', 'Fano (\Delta = 1)');
%hold on
%bar(SNRdB, Iterations(2, :), 0.45, 'FaceColor', '#FF0000', 'DisplayName', 'Fano (\Delta = 1)')
%bar(SNRdB, Iterations(3, :), 0.45, 'FaceColor', '#191970', 'DisplayName', 'List (L = 16)')
%hold off
set(gca, 'YScale', 'log', 'YLim', [1e2 1e6], 'YGrid', 'on', 'XGrid', 'on', ...
    'XColor', [.1 .1 .1])
xlabel('E_b/N_0(dB)')
ylabel('AVN')
legend()
%}


