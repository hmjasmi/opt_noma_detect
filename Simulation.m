% This code simulates the NOMA system and plots the BER of the far user
% (U1) and near user (U2). The code generates the BER of the simulated conventional SIC and optimal detector (derived by the paper)
function [BER2_conv,BER2_novl]=Simulation(SNR_dB,alpha,Omega_2)
%% Parameters
SNR = 10.^(SNR_dB./10);  % Avgerage SNR - linear
A11 = sqrt(alpha(1)) + sqrt(alpha(2));   % NOMA amplitude A11
A10 = sqrt(alpha(1)) - sqrt(alpha(2));    % NOMA amplitude A10
Iter = 1e5;                    % Number of iterations
Omega_1=2;
%%
for i_snr = 1:length(SNR)
    snr = SNR(i_snr);
    sd = 1./sqrt(2*snr); % Noise standard deviation
    % Initialization
    iter = 0;
    BER_U1_conv = 0; % simulated BER for conventional SIC of U1
    BER_U2_conv = 0; % simulated BER for conventional SIC of U2
    BER_U2_novl = 0; % simulated BER for optimal detector  of U2
    E1_count = 0;
    E0_count = 0;
    while (iter < Iter)
        h_1 = sqrt(Omega_1)./sqrt(2)*(randn(1,1) + 1j*randn(1,1)); % complex channnel coefficient of U1
        h_2 = sqrt(Omega_2)./sqrt(2)*(randn(1,1) + 1j*randn(1,1)); % complex channnel coefficient of U2
        beta_1 = abs(h_1);    % beta_1=|h_1| abs value of the channel
        beta_2 = abs(h_2);    % beta_2=|h_2| abs value of the channel
        w_1 = sd*(randn(1,1) + 1j*randn(1,1));  % AWGN at U1
        w_2 = sd*(randn(1,1) + 1j*randn(1,1));  % AWGN at U2

        data1 = rand(1,1)>0.5;    % Transmitted signal by U1 (bits)
        data2 = rand(1,1)>0.5;    % Transmitted signal by U2 (bits)
        x1 = 2*data1 - 1;         % Transmitted signal by U1 (modulated)
        x2 = 2*data2 - 1;         % Transmitted signal by U2 (modulated)

        x = sqrt(alpha(1)).*x1 + sqrt(alpha(2)).*x2;  % Transmitted NOMA signal
        y_1 = beta_1.*x + w_2;                 % Received signal at U1
        y_2 = beta_2.*x + w_2;                 % Received signal at U2
        x_trials = [-1; +1];             % BPSK constellation
        %% U1 Data Detection
        data1_arg_U1 = abs(y_1 - beta_1.*sqrt(alpha(1)).*x_trials).^2;
        [~,data1_idx_U1] = min(data1_arg_U1,[],1);   % MLD to detect U1's signal at U1
        data1_detected_U1 = data1_idx_U1 - 1; % Demodulation of U1's signal at U1
        BER_U1_conv = BER_U1_conv + sum(data1_detected_U1~=data1);   %  BER of the convetional SIC detector at U1

        %% U2 Data Detection
        data1_arg_U2 = abs(y_2 - beta_2.*sqrt(alpha(1)).*x_trials).^2;
        [~,data1_idx] = min(data1_arg_U2,[],1);   % MLD to detect U1's signal at U2
        data1_detected = data1_idx - 1;             % Demodulation of U1's signal at U2
        x1_hat_U2 = 2*data1_detected - 1;
        %-------------------------------- Conv. SIC ----------------------------
        y_sic = y_2 - beta_2.*sqrt(alpha(1)).*x1_hat_U2;
        data2_arg_conv = abs(y_sic - beta_2.*sqrt(alpha(2)).*x_trials).^2;
        [~,data2_idx_conv] = min(data2_arg_conv, [],1);
        data2_detected_conv = data2_idx_conv - 1;
        BER_U2_conv = BER_U2_conv + sum(data2_detected_conv~=data2);  % BER of the convetional SIC detector at U2
        E1 = y_sic>0;
        %-------------------------------- Optimal Detector ----------------------------

        RHS = sd.^2/4/beta_2/sqrt(alpha(2))*log((1 - exp(-4*abs(y_2)*beta_2*A11/sd.^2))./(1 - exp(-4*abs(y_2)*beta_2*A10/sd.^2)));
        data2_detected_novl = y_sic>-sign(y_2)*RHS;
        BER_U2_novl = BER_U2_novl + sum(data2_detected_novl~=data2);
        iter = iter + 1;
    end
    BER1_conv(i_snr) = BER_U1_conv./iter;
    BER2_conv(i_snr) = BER_U2_conv./iter;
    BER2_novl(i_snr) = BER_U2_novl./iter;

end
