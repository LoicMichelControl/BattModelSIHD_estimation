# Semi-implicit estimation of a battery model using the SIHD method — Preliminary results
---

## Overview

This repository presents a proof-of-concept application of the  
Semi-Implicit Homogeneous Differentiation (SIHD) technique to the estimation of a  
Li-ion battery equivalent circuit model, composed of four parallel first-order ODEs.

The simulated battery follows a classical 2-RC equivalent circuit with thermal dynamics.

---

## Battery model

### **State equations**

$$
\begin{aligned}
\dot{\mathrm{SoC}} &= -\frac{I}{Q_{\text{nom}}}, \\[6pt]
\dot{V}_{RC1} &= -\frac{V_{RC1}}{R_1 C_1} + \frac{I}{C_1}, \\[6pt]
\dot{V}_{RC2} &= -\frac{V_{RC2}}{R_2 C_2} + \frac{I}{C_2}, \\[6pt]
\dot{T} &= \frac{I^2\,(R_s + R_1 + R_2) - (T - T_{\text{amb}})/R_{\text{th}}}{C_{\text{th}}}.
\end{aligned}
$$

These equations describe the evolution of:
- the State of Charge (SoC),  
- the RC branch voltage \(V_{RC1}\),  
- the RC branch voltage \(V_{RC2}\),  
- the temperature \(T\).

The objective is to analyse the behavior of the SIHD observer when applied to 
a multi-state, decoupled, current-driven system.

---

## SIHD-based observer principle

To obtain a prediction of each state, SIHD relies on a zero-initialized duplicate model, acting as a *minimal* Luenberger-type observer, which:
- does not attempt to match the real initial conditions;
- only provides a prediction error, used by SIHD to update its internal states;
Each ODE branch (SoC, RC1, RC2, T) is treated independently.

In this preliminary work, the terminal voltage \(V_t\) is not injected into the SIHD loop.  
The primary objective is to investigate the intrinsic convergence properties of the discrete SIHD algorithm before introducing measurement-fed synchronization in future developments.

---

## Results

<p align="center">
  <img width="90%" src="/Figures/Fig_1.png" alt="Estimated states">
</p>
<p align="center"><em>Fig. 1 – Estimated states. Good tracking performance is observed despite noisy current excitation.</em></p>

<p align="center">
  <img width="90%" src="/Figures/Fig_2.png" alt="Estimated derivatives">
</p>
<p align="center"><em>Fig. 2 – Estimated time derivatives of the states.</em></p>

---

## References

[1] V. Acary, B. Brogliato, and Y. Orlov.  
*Chattering-free digital sliding-mode control with state observer and disturbance rejection.*  
IEEE Trans. on Automatic Control, 57(5):1087–1101, 2012.

[2] L. Michel, M. Ghanes, Y. Aoustin, and J.-P. Barbot.  
*An interconnected discrete-time cascaded semi-implicit differentiation.*  
17th International Workshop on Variable Structure Systems (VSS 2024), Accepted.  
<https://hal.science/hal-04564290/>

[3] L. Michel, M. Ghanes, Y. Aoustin, J.-P. Barbot.  
*Semi-implicit discrete-time cascaded observer: Chua circuit case study.*  
In: Sliding mode, variable-structure discontinuous control systems. Springer, In press.  
<https://hal.science/hal-04993406>

---

## Credits and License

(c) 2025 — Nantes Université / Centrale Nantes  
Laboratoire des Sciences du Numérique de Nantes — LS2N UMR 6004, France

Author: Loïc Michel

Distributed under the Creative Commons Attribution–NonCommercial–ShareAlike 4.0 International license.
