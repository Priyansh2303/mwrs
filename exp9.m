clc;
close all;

% Number of samples
N = 20000;
SNR_db = -5:0.5:35;
SNR = 10.^(SNR_db/10);

% Generate Rayleigh random variable
r = sqrt(-2 * log(rand(1, N))); % Simplified formula

% Plot histogram
figure;
histogram(r, 100);
title('Rayleigh Random Variable Histogram');
xlabel('Random Variable');
ylabel('Frequency');

% Plot Rayleigh PDF
x = 0:0.01:10;
rayleigh_pdf = (x) .* exp(-x.^2 / 2);
figure;
plot(x, rayleigh_pdf);
title('Rayleigh PDF');
xlabel('Random Variable');
ylabel('Probability');
legend('Variance = 1');

% Error Probability Calculations
Pe_BPSK = 0.5 * (1 - sqrt(SNR ./ (1 + SNR)));
Pe_BFSK = 0.5 * (1 - sqrt(SNR ./ (2 + SNR)));
Pe_DPSK = 0.5 ./ (1 + SNR);

% AWGN Channel Error Probabilities
Pe_BPSK_AWGN = 0.5 * erfc(sqrt(SNR));
Pe_BFSK_AWGN = 0.5 * erfc(sqrt(SNR / 2));
Pe_DPSK_AWGN = 0.5 * exp(-SNR);

% Plot Performance Comparison
figure;
semilogy(SNR_db, Pe_BPSK, 'r.-', SNR_db, Pe_BFSK, 'r*-', SNR_db, Pe_DPSK, 'r--', ...
         SNR_db, Pe_BPSK_AWGN, 'b.-', SNR_db, Pe_BFSK_AWGN, 'b*-', SNR_db, Pe_DPSK_AWGN, 'b--');
title('Performance of BPSK, BFSK, and DPSK in Rayleigh Fading');
xlabel('SNR (dB)');
ylabel('Probability of Error');
legend('BPSK with Fading', 'BFSK with Fading', 'DPSK with Fading', ...
       'BPSK without Fading', 'BFSK without Fading', 'DPSK without Fading');
grid on;

%% Monte Carlo Simulation for BPSK in Rayleigh Fading
BER_BPSK = zeros(1, length(SNR_db));
for i = 1:length(SNR_db)
    errors = 0;
    bits = 0;
    while errors <= 10
        alpha = sqrt(-2 * log(rand)); % Rayleigh fading
        noise = sqrt(1/(2*SNR(i))) * randn; % AWGN noise
        y = alpha + noise;
        errors = errors + (y < 0);
        bits = bits + 1;
    end
    BER_BPSK(i) = errors / bits;
end

figure;
semilogy(SNR_db, BER_BPSK, '-o');
title('Monte Carlo Simulation for BPSK in Rayleigh Fading');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
grid on;

%% Monte Carlo Simulation for BFSK in Rayleigh Fading
BER_BFSK = zeros(1, length(SNR_db));
for i = 1:length(SNR_db)
    errors = 0;
    bits = 0;
    while errors <= 10
        alpha = sqrt(-2 * log(rand));
        noise = sqrt(1/(2*SNR(i))) * randn;
        y = alpha + noise;
        errors = errors + (y < 0);
        bits = bits + 1;
    end
    BER_BFSK(i) = errors / bits;
end

figure;
semilogy(SNR_db, BER_BFSK, '-*');
title('Monte Carlo Simulation for BFSK in Rayleigh Fading');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
grid on;
