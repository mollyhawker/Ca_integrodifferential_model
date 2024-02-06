# A calcium puff model based on integrodifferential equations (Hawker et al., 2024)

Ca_integrodifferential_model implements three hybrid stochastic systems used to simulate the correct ion channel and calcium signalling dynamics. We ask that reference be made to our manuscript **A calcium puff model based on integrodifferential equations** (<https://doi.org/10.48550/arXiv.2401.17326>) when using this code. All code is contained within `ca_puff_model.zip`.

## 6-state model

The 6-state model is based on the Markov model by Siekmann et al., (2012) and Cao et al., (2013). We replace the four ODEs used to model the gating variables in the Cao et al., (2013) model with integrodifferential equations, based on the method by Brady (1972).

The 6-state model is simulated in `Hawker_et_al_ca_puff_model_6_state.m`

## Reduced 6-state model

Our 6-state model is reduced to a 6-state model with one gating variable modelled using the integral terms. This is achieved by using quasi-steady-state approximation and setting three of the gating variables (m24, h24, m42) immediately to their steady state. This is possible as the rate at which the three gating variables reach equilibrium is extremely fast.

The reduced 6-state model is simulated in `Hawker_et_al_ca_puff_model_reduced_6_state.m`

## Reduced 2-state model

We further simplify the reduced 6-state model by using quasi-steady-state approximation and ignoring low dwell times.

The reduced 2-state model is simulated in `Hawker_et_al_ca_puff_model_2_state.m`

## The integral term 

Our integral term is a finite version of the integrodifferential equation by Brady (1972). We interpret the integral as a temporal weighted average over past calcium concentrations. Using our model we can analyse the effect different integral lengths have on the calcium puff dynamics and relate this to the "memory" of the ion channel. Within the equation $G$ represents the gating variables, $\lambda_{G}$ represents the rate the gating variable reaches equilibrium, $\alpha_{G}=\lambda_{G}G_{\infty}$, $c$ represents calcium, $\tau$ represents how far into the past the ion channel can *recognise* calcium concentrations.

$$
\Phi_{G}(t,c)=G(0)exp\left[-\int_{t-\tau}^{t}(\lambda_{G}\circ c)(x) dx\right]    - exp\left[- \int_{t-\tau}^{t}  (\lambda_{G}\circ c)(x) dx\right] \int_{t-\tau}^{t}(-\alpha_{G} \circ c)(s) \times \hspace{0.5cm} exp\left[\int_{s-\tau}^{s}(\lambda_{G}\circ c)(x) dx\right]ds
$$

The integral is implemented in `gatingSolutionMH.m`. `Ca_model.m` and `RK4.m` are required for the numerical solutions within the puff models.
