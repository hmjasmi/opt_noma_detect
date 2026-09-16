# On the Detection Optimality and Exact BER Analysis of NOMA
# MATLAB Code for the BER Evaluation of Two-User NOMA Detectors

## Overview

This package evaluates the bit-error rate (BER) at the near user, (U\_2), in a two-user downlink NOMA system with BPSK signaling. It compares:

* Monte Carlo simulation of conventional successive interference cancellation (SIC);
* Monte Carlo simulation of the optimal detector;
* The theoretical BER of conventional SIC/joint maximum-likelihood detection; and
* The theoretical BER of the optimal detector, including numerical computation of the decision threshold (\\tau\_\\beta) obtained from \[1]. 

## Files

|File|Purpose|
|-|-|
|`Main.m`|Defines the system parameters, calls the simulation and theoretical functions, and plots the BER curves.|
|`Simulation.m`|Simulates transmission over independent Rayleigh-fading channels and returns the BER of conventional SIC and the optimal detector at U\_2.|
|`Theory\_Conv\_SIC\_JML.m`|Evaluates the closed-form theoretical BER of the conventional SIC/JML detector.|
|`Theory\_Optimal\_Detector.m`|Evaluates the theoretical BER of the optimal detector by numerically finding (\\tau(\\beta)) and integrating over Rayleigh fading.|

Keep all four MATLAB files in the same folder.

## Requirements

* MATLAB; R2021a or later is recommended because `Main.m` uses the `Name=Value` plotting syntax.
* No specialized toolbox is required; the implementation uses standard MATLAB functions such as `randn`, `fzero`, `erfc`, `trapz`, and `semilogy`.

For older MATLAB releases, replace expressions such as

```matlab
xlabel('$\\overline{\\gamma}$ (dB)',Interpreter='latex')
```

with

```matlab
xlabel('$\\overline{\\gamma}$ (dB)','Interpreter','latex')
```

## Running the code

1. Place the four `.m` files in one MATLAB working directory.
2. Open MATLAB and change the current folder to that directory.
3. Run:

```matlab
Main
```

The script calculates the four BER curves and displays them on a semilogarithmic plot.

## Main parameters

The principal parameters are defined at the beginning of `Main.m`:

```matlab
SNR\_dB = -20:2:20;
alpha  = \[0.8 1-0.8];
Omega  = 1;
```

* `SNR\_dB` is the average signal-to-noise-ratio range in decibels.
* `alpha(1)` and `alpha(2)` are the NOMA power-allocation coefficients. They must be positive and satisfy `sum(alpha)=1` and 'alpha(1)>alpha(2)'.
* `Omega` is the average channel power used for (U\_2), so that (\\mathbb{E}\[|h\_2|^2]=\\Omega).

The composite BPSK amplitudes are

\[
A\_{11}=\\sqrt{\\alpha\_1}+\\sqrt{\\alpha\_2},\\qquad
A\_{10}=\\sqrt{\\alpha\_1}-\\sqrt{\\alpha\_2}.
]



## Program flow

`Main.m` calls the three computational functions as follows:

```matlab
\[BER2\_conv,BER2\_novl] = Simulation(SNR\_dB,alpha,Omega);
BER\_OPT\_theory         = Theory\_Optimal\_Detector(SNR\_dB,alpha,Omega);
BER\_Conv\_Theory        = Theory\_Conv\_SIC\_JML(SNR\_dB,alpha,Omega);
```

The returned quantities are:

|Variable|Meaning|
|-|-|
|`BER2\_conv`|Simulated BER at (U\_2) using conventional SIC.|
|`BER2\_novl`|Simulated BER at (U\_2) using the optimal detector.|
|`BER\_OPT\_theory`|Theoretical BER of the optimal detector.|
|`BER\_Conv\_Theory`|Theoretical BER of conventional SIC/JML detection.|

## Monte Carlo simulation

`Simulation.m` performs `Iter = 1e5` independent trials at each SNR value. In every trial it:

1. generates independent Rayleigh-fading coefficients for (U\_1) and (U\_2);
2. generates two equiprobable BPSK symbols;
3. forms the NOMA signal

   \[
x=\\sqrt{\\alpha\_1}x\_1+\\sqrt{\\alpha\_2}x\_2;
]

4. adds complex Gaussian noise;
5. detects the first user's symbol at (U\_2);
6. applies conventional SIC and detects the second user's symbol;
7. applies the optimal decision rule and detects the second user's symbol; and
8. estimates each BER as the number of errors divided by `Iter`.

The real and imaginary noise components have variance

\[
\\sigma\_{\\mathrm{real}}^2=\\frac{1}{2\*\\mathrm{SNR}}.
]

Because no random-number seed is set, the simulated curves vary slightly between runs. For exactly repeatable results, insert the following line in `Main.m` before calling `Simulation`:

```matlab
rng(1)
```

## Conventional theoretical detector

`Theory\_Conv\_SIC\_JML.m` calculates the closed-form BER using terms of the form

\[
M(d)=\\sqrt{\\frac{d^2\\bar\\gamma}{d^2\\bar\\gamma+2}},
\\qquad
\\bar\\gamma=\\frac{\\Omega}{\\sigma\_{\\mathrm{real}}^2}.
]

The final expression combines the terms associated with (A\_{11}), (A\_{10}), (2\\sqrt{\\alpha\_2}), and the remaining composite amplitudes. The final BER expression is provided in \[1, Eq. (45)]. 

## Optimal theoretical detector and numerical threshold

`Theory\_Optimal\_Detector.m` averages the conditional error probabilities over Rayleigh fading. It uses the transformation

\[
u=\\frac{\\beta^2}{\\Omega},\\qquad
f\_\\beta(\\beta),d\\beta=e^{-u},du,
]

and evaluates the resulting integrals numerically over `u = 0` to `uMax = 45` using `Nu = 7000` samples and `trapz`.

For each fading amplitude (\\rho=\\beta), the code obtains the nonnegative threshold (\\tau(\\rho)) in \[1, Eqs. (28),(29)]. Defining

\[
\\Delta=A\_{11}^2-A\_{10}^2,
]

the root function implemented in the code is

\[
F(t)=\\log!\\left\[\\sinh!\\left(\\frac{2t\\rho A\_{11}}{\\sigma\_n^2}\\right)\\right]
-\\log!\\left\[\\sinh!\\left(\\frac{2t\\rho A\_{10}}{\\sigma\_n^2}\\right)\\right]
-\\frac{\\rho^2\\Delta}{\\sigma\_n^2}.
]

Instead of evaluating `log(sinh(x))` directly, the program uses the numerically stable identity

\[
\\log(\\sinh x)=x-\\log 2+\\log(1-e^{-2x}),\\qquad x>0.
]

This is represented in MATLAB by `log1p(-exp(...))`, which reduces loss of numerical accuracy.

The threshold is selected as follows:

* At the boundary (t\\rightarrow0^+), the code evaluates

  \[
F(0^+)=\\log!\\left(\\frac{A\_{11}}{A\_{10}}\\right)
-\\frac{\\rho^2\\Delta}{\\sigma\_n^2}.
]

* If `F0 >= 0`, the constrained nonnegative threshold is set to `tau = 0`.
* If `F0 < 0`, the code brackets the sign-changing positive root. It repeatedly doubles the upper endpoint until `F(t\_high) >= 0`, then calls

```matlab
  tau = fzero(F,\[t\_low t\_high]);
  ```

  to obtain the desired positive root.

The resulting threshold is used in the four conditional Gaussian-CDF terms corresponding to successful and failed SIC. These quantities are integrated over the Rayleigh distribution and combined using the probabilities of SIC success and failure. The BER of the optimal detector at U\_2 is provided in (44). 

## Numerical settings

|Setting|Location|Effect|
|-|-|-|
|`Iter = 1e5`|`Simulation.m`|Larger values reduce Monte Carlo variation but increase runtime.|
|`uMax = 45`|`Theory\_Optimal\_Detector.m`|Truncates the transformed Rayleigh integral; the omitted exponential tail is very small.|
|`Nu = 7000`|`Theory\_Optimal\_Detector.m`|Controls the integration-grid resolution and runtime.|
|`t\_high` doubling|`Theory\_Optimal\_Detector.m`|Ensures that the positive threshold root is bracketed before `fzero` is called.|

## Changing an experiment

* Change the SNR range by editing `SNR\_dB` in `Main.m`.
* Change the power allocation by editing `alpha`, while maintaining positive coefficients that sum to one.
* Improve Monte Carlo accuracy by increasing `Iter` in `Simulation.m`.
* Improve numerical-integration resolution by increasing `Nu` in `Theory\_Optimal\_Detector.m`.

After changing a numerical parameter, rerun `Main.m` and compare the new results with the previous curves to confirm convergence.

## Citation

\[1] The optimal detector derivation and BER are obtained from T. Assaf, H. Yahya, and A. Al-Dweik, “On the detection optimality and exact BER analysis of NOMA,” arXiv preprint arXiv:2607.17755, Jul. 2026, doi: 10.48550/arXiv.2607.17755. \[Online]. Available: https://arxiv.org/abs/2607.17755.
