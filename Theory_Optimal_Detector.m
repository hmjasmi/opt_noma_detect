%% This code presents the theortical BER of the optimal detector at U2.
function [BER_OPT_theory]=Theory_Optimal_Detector(SNR_dB,alpha,Omega)
SNR = 10.^(SNR_dB./10);  % Avgerage SNR - linear
a = sqrt(alpha(1));
b = sqrt(alpha(2));
A11 = a + b;
A10 = a - b;
sigma_n2 = 1 ./ (2*SNR);   % Noise variance
gamma_bar = Omega ./ sigma_n2;   % Average SNR

%% Allocate arrays for the BERs
BER_OPT_theory = zeros(size(SNR));
PeS_OPT = zeros(size(SNR));
PeF_OPT = zeros(size(SNR));
PrS_vec = zeros(size(SNR));
PrF_vec = zeros(size(SNR));
%% Integration grid
uMax = 45;   % To replace the integration upper limit (exp(-45)\approx 2.86e-20.
Nu = 7000;   % Integration points
u = linspace(1e-10, uMax, Nu);  % to create u - grid
rayleigh_weight = exp(-u);  % transformed Rayleigh density
betaS_grid = sqrt(Omega .* u);   % Success u - grid
betaF_grid = sqrt(Omega .* u);  % Failure u - grid

Delta = A11^2 - A10^2;     % Delta = 4*sqrt(alpha(1)*alpha(2))

for i_snr = 1:length(SNR)
    sig2 = sigma_n2(i_snr);   % sigma_n^2
    sig = sqrt(sig2);       %  sigma_n (standard deviation)
    gamma_i = gamma_bar(i_snr);

    % --------------------------------------------------------
    % Success/failure probabilities
    % --------------------------------------------------------

    M11 = sqrt((A11^2 * gamma_i) / (A11^2 * gamma_i + 2));   % Eq. (32)and M(d)
    M10 = sqrt((A10^2 * gamma_i) / (A10^2 * gamma_i + 2));   % Eq. (32)and M(d)

    p11 = (1 + M11) / 2;      % Eq. (33)
    p10 = (1 + M10) / 2;       % Eq. (34)

    q11 = 1 - p11;    
    q10 = 1 - p10;

    PrS = (p11 + p10) / 2;  % Eq. (35) 
    PrF = (q11 + q10) / 2;   % Eq. (35) 

    PrS_vec(i_snr) = PrS;
    PrF_vec(i_snr) = PrF;

    %% Compute tau_beta for every fading amplitude
 
    tauS_grid = zeros(size(betaS_grid));
    tauF_grid = zeros(size(betaF_grid));

    for k = 1:length(u)

        % ====================================================
        % Success branch: tau(beta_S)
        % ====================================================

        rho = betaS_grid(k);    % Temp variable for beta_S

        if rho < 1e-14

            tauS_grid(k) = 0;    %  To avoid the division by a value numerically equal to zero

        else
            % F_beta(t)=log sinh(2*t*beta A11/\sigma_n^2)-log sinh(2*t*beta
            % A10/\sigma_n^2)-beta^2*Delta/sigma_n^2 >> % Equivalent to Eqs. (26),(28) and (29)
            F0 = log(A11/A10) - rho^2 * Delta / sig2;    % F_beta(t)=0 

            if F0 >= 0     % To satisify the boundary limit, i.e., if Delta_beta <= A11/A10, then tau_beta=0 

                tauS_grid(k) = 0;

            else

                t_low = 1e-6 * sig2 / (2*rho*A11);  % t_low is set to a small positive value to avoid undefined numerical problem
                t_high = rho * A11;

                F = @(t) ...
                    (2*t*rho*A11/sig2 ...
                    - log(2) ...
                    + log1p(-exp(-4*t*rho*A11/sig2))) ...
                    - ...
                    (2*t*rho*A10/sig2 ...
                    - log(2) ...
                    + log1p(-exp(-4*t*rho*A10/sig2))) ...
                    - rho^2 * Delta / sig2;

                while F(t_high) < 0
                    t_high = 2*t_high;

                    if t_high > 1e6
                        error('Could not bracket tau(beta_S).');
                    end
                end

                tauS_grid(k) = fzero(F, [t_low, t_high]);

            end
        end

        % ====================================================
        % Failure branch: tau(beta_F)
        % ====================================================


        rho = betaF_grid(k);

        if rho < 1e-14

            tauF_grid(k) = 0;

        else

            F0 = log(A11/A10) - rho^2 * Delta / sig2;

            if F0 >= 0

                tauF_grid(k) = 0;

            else

                t_low = 1e-6 * sig2 / (2*rho*A11);
                t_high = rho * A11;

                F = @(t) ...
                    (2*t*rho*A11/sig2 ...
                    - log(2) ...
                    + log1p(-exp(-4*t*rho*A11/sig2))) ...
                    - ...
                    (2*t*rho*A10/sig2 ...
                    - log(2) ...
                    + log1p(-exp(-4*t*rho*A10/sig2))) ...
                    - rho^2 * Delta / sig2;

                while F(t_high) < 0
                    t_high = 2*t_high;

                    if t_high > 1e6
                        error('Could not bracket tau(beta_F).');
                    end
                end

                tauF_grid(k) = fzero(F, [t_low, t_high]);

            end
        end

    end

    % --------------------------------------------------------
    % Gaussian CDF terms
    % Phi_G(x/sigma_n) = 0.5*erfc(-x/(sqrt(2)*sigma_n))
    % --------------------------------------------------------

    % Successful SIC branch: beta_S, W

    Phi_S_A11 = ...
        0.5 .* erfc(-(betaS_grid .* A11) ./ (sqrt(2)*sig));

    Phi_S_A11_minus_tau = ...
        0.5 .* erfc(-(betaS_grid .* A11 - tauS_grid) ./ (sqrt(2)*sig));

    Phi_S_A10_minus_tau = ...
        0.5 .* erfc(-(betaS_grid .* A10 - tauS_grid) ./ (sqrt(2)*sig));

    % Failed SIC branch: beta_F, Z

    Phi_F_A11_plus_tau = ...
        0.5 .* erfc(-(betaF_grid .* A11 + tauF_grid) ./ (sqrt(2)*sig));

    Phi_F_A10_plus_tau = ...
        0.5 .* erfc(-(betaF_grid .* A10 + tauF_grid) ./ (sqrt(2)*sig));

    Phi_F_A10 = ...
        0.5 .* erfc(-(betaF_grid .* A10) ./ (sqrt(2)*sig));

   % I_11^S:
    % A00/A11 contribution under successful SIC

    I11S_integrand = ...
        Phi_S_A11 ...
        - Phi_S_A11_minus_tau;

    % I_10^S:
    % A01/A10 contribution under successful SIC

    I10S_integrand = ...
        Phi_S_A10_minus_tau;

    % I_11^F:
    % A00/A11 contribution under failed SIC

    I11F_integrand = ...
        1 - Phi_F_A11_plus_tau;

    % I_10^F:
    % A01/A10 contribution under failed SIC

    I10F_integrand = ...
        Phi_F_A10_plus_tau ...
        - Phi_F_A10;

    % --------------------------------------------------------
    % Integrals using f_beta(beta_S)d_beta_S and
    % f_beta(beta_F)d_beta_F through the u transformation
    % --------------------------------------------------------

    I11S = trapz(u, I11S_integrand .* rayleigh_weight);
    I10S = trapz(u, I10S_integrand .* rayleigh_weight);

    I11F = trapz(u, I11F_integrand .* rayleigh_weight);
    I10F = trapz(u, I10F_integrand .* rayleigh_weight);

    % --------------------------------------------------------
    % Conditional BERs
    % --------------------------------------------------------

    PeS_OPT(i_snr) = (I11S + I10S) / (p11 + p10);
    PeF_OPT(i_snr) = (I11F + I10F) / (q11 + q10);

    % --------------------------------------------------------
    % Total BER
    % --------------------------------------------------------

    BER_OPT_theory(i_snr) = ...
        PeS_OPT(i_snr)*PrS ...
        + PeF_OPT(i_snr)*PrF;
 

end

BER_OPT_theory = max(BER_OPT_theory, 0);

