%% EXPERIMENT 10 - SIMULATION OF THE OFDM TRANSMITTER AND RECEIVER (Simplified)

clc;
clear;

% Define system parameters
nFFT = 64; % FFT size
nDSC = 52; % Data subcarriers
nBitPerSym = nDSC; % BPSK uses one bit per subcarrier
nSym = 1e4; % Number of symbols
EbN0dB = 0:10; % SNR range
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % Convert Eb/N0 to Es/N0

simBer = zeros(size(EbN0dB)); % Store BER results

for ii = 1:length(EbN0dB)
    % Transmitter
    ipBit = randi([0 1], 1, nBitPerSym * nSym); % Generate random bits
    ipMod = 2 * ipBit - 1; % BPSK modulation
    ipMod = reshape(ipMod, nBitPerSym, nSym).';
    
    % OFDM modulation (insert into subcarriers)
    xF = [zeros(nSym, 6), ipMod(:, 1:nBitPerSym/2), zeros(nSym, 1), ipMod(:, nBitPerSym/2+1:nBitPerSym), zeros(nSym, 5)];
    xt = ifft(fftshift(xF, 2), nFFT, 2) * (nFFT / sqrt(nDSC)); % IFFT operation
    xt = [xt(:, 49:64), xt]; % Append cyclic prefix
    xt = reshape(xt.', 1, []); % Serialize data
    
    % Add noise
    noise = (randn(size(xt)) + 1j * randn(size(xt))) / sqrt(2); % AWGN noise
    yt = sqrt(80/64) * xt + 10^(-EsN0dB(ii)/20) * noise;
    
    % Receiver
    yt = reshape(yt, 80, []).';
    yt = yt(:, 17:end); % Remove cyclic prefix
    yF = fftshift(fft(yt, nFFT, 2), 2) * sqrt(nDSC) / nFFT; % FFT operation
    yMod = yF(:, [6+(1:nBitPerSym/2), 7+(nBitPerSym/2+1:nBitPerSym)]);
    
    % BPSK demodulation
    ipBitHat = (real(yMod) > 0); % Convert to bits
    ipBitHat = reshape(ipBitHat.', 1, []);
    
    % Count errors
    simBer(ii) = sum(ipBitHat ~= ipBit) / length(ipBit);
end

% Theoretical BER
theoryBer = (1/2) * erfc(sqrt(10.^(EbN0dB/10)));

% Plot results
figure;
semilogy(EbN0dB, theoryBer, 'bs-', 'LineWidth', 2);
hold on;
semilogy(EbN0dB, simBer, 'mx-', 'LineWidth', 2);
grid on;
legend('Theory', 'Simulation');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate');
title('BER Curve for BPSK using OFDM');
