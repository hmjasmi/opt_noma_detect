%% This code presents the theortical BER of the convetional SIC/JML detection  at U2.
function [BER_Conv_Theory]=Theory_Conv_SIC_JML(SNR_dB,alpha,Omega)
%% Parameters
SNR = 10.^(SNR_dB./10);
a = sqrt(alpha(1));
b = sqrt(alpha(2));
A11 = a + b;  % NOMA amplitude 
A10 = a - b;    % NOMA amplitude
sigma_n2 = 1 ./ (2*SNR);     % Noise var - real 
gamma_bar = Omega ./ sigma_n2;   % average SNR
%% ============================================================
% Conventional SIC/JML theory
% ============================================================

M_A11 = sqrt((A11.^2 .* gamma_bar) ./ (A11.^2 .* gamma_bar + 2));
M_A10 = sqrt((A10.^2 .* gamma_bar) ./ (A10.^2 .* gamma_bar + 2));

M_b = sqrt((b.^2 .* gamma_bar) ./ (b.^2 .* gamma_bar + 2));

M_2apb = sqrt(((2*a + b).^2 .* gamma_bar) ./ (((2*a + b).^2 .* gamma_bar) + 2));
M_2amb = sqrt(((2*a - b).^2 .* gamma_bar) ./ (((2*a - b).^2 .* gamma_bar) + 2));

BER_Conv_Theory = ...
    (2 ...
    + M_A11 ...
    - M_A10 ...
    + M_2amb ...
    - M_2apb ...
    - 2*M_b) ./ 4;
