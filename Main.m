clear;
clc;
% close all;
%% Para
SNR_dB = -20:2:20;           % Average SNR - dB
alpha = [0.8 1-0.8];       % Power allocation   
Omega = 1;   % Channel variance

%% 
[BER2_conv,BER2_novl]=Simulation(SNR_dB,alpha,Omega);  % Monte-Carlo Simulations for U_2 using the conv. SIC and optimal detector.
[BER_OPT_theory]=Theory_Optimal_Detector(SNR_dB,alpha,Omega); % Theortical BER for the optimal detector 
[BER_Conv_Theory]=Theory_Conv_SIC_JML(SNR_dB,alpha,Omega);    % Theortical BER for the conv. SIC 
%% Plotting
figure;
semilogy(SNR_dB,BER2_conv,'ro');
hold on;
semilogy(SNR_dB,BER2_novl,'bp'); 
hold on;
semilogy(SNR_dB,BER_OPT_theory,'k--'); 
hold on;
semilogy(SNR_dB,BER_Conv_Theory,'k--')
legend('Conv. SIC','U_2: Opt. Detector','Theory')
xlabel('$\bar{\gamma}$ (dB)',Interpreter='latex')
ylabel('BER')
grid;
ylim([1e-3 1e0])
