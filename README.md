# A calcium puff model based on integrodifferential equations (Hawker et al., 2024)

Ca_integrodifferential_model implements three versions of a model for simulating stochastic intracellular calcium release based on hybrid stochastic systems. The model is developed in the manuscript "A calcium puff model based on integrodifferential equations” by Molly Hawker, Pengxing Cao, James Sneyd and Ivo Siekmann” available on [arXiv](https://doi.org/10.48550/arXiv.2401.17326). We ask that reference be made to our manuscript when using this code: Molly Hawker, Pengxing Cao, James Sneyd and Ivo Siekmann (2024). A calcium puff model based on integrodifferential equations, arXiv, arXiv:2401.17326, [https://doi.org/10.48550/arXiv.24](https://doi.org/10.48550/arXiv.2401.17326){.uri}

## 6-state model

The 6-state model is based on the Markov model by Siekmann et al., (2012) and Cao et al., (2013). We replace the four ODEs used to model the gating variables in the Cao et al., (2013) model with integrodifferential equations, based on the method by Brady (1972).

The 6-state model is simulated in `Hawker_et_al_ca_puff_model_6_state.m`. This model can be used to reproduce the calcium trace shown in Fig 4a (Hawker et al., 2024).

## Reduced 6-state model

Our 6-state model is reduced to a 6-state model with one gating variable modelled using the integral terms. This is achieved by using quasi-steady-state approximation and setting three of the gating variables (m24, h24, m42) immediately to their steady state. This is possible as the rate at which the three gating variables reach equilibrium is extremely fast.

The reduced 6-state model is simulated in `Hawker_et_al_ca_puff_model_reduced_6_state.m` . Using this model, one can reproduce the calcium trace shown in Fig 4b (Hawker et al., 2024).

## Reduced 2-state model

We further simplify the reduced 6-state model by using quasi-steady-state approximation and ignoring low dwell times.

The reduced 2-state model is simulated in `Hawker_et_al_ca_puff_model_2_state.m`. This model can be used to reproduce the calcium trace shown in Fig 4c (Hawker et al., 2024).

## The integral term

Our integral term is a finite version of the integrodifferential equation by Brady (1972). We interpret the integral as a temporal weighted average over past calcium concentrations. Using our model we can analyse the effect different integral lengths have on the calcium puff dynamics and relate this to the "memory" of the ion channel. Within the equation $G$ represents the gating variables, $\lambda_{G}$ represents the rate the gating variable reaches equilibrium, $\alpha_{G}=\lambda_{G}G_{\infty}$, $c$ represents calcium, $\tau$ represents how far into the past the ion channel can *recognise* calcium concentrations.

$$
\Phi_{G}(t,c)=G(0)exp\left[-\int_{t-\tau}^{t}(\lambda_{G}\circ c)(x) dx\right]    - exp\left[- \int_{t-\tau}^{t}  (\lambda_{G}\circ c)(x) dx\right] \int_{t-\tau}^{t}(-\alpha_{G} \circ c)(s) \times \hspace{0.5cm} exp\left[\int_{s-\tau}^{s}(\lambda_{G}\circ c)(x) dx\right]ds
$$

The integral is implemented in `gatingSolutionMH.m`. `Ca_model.m` and `RK4.m` are required for the numerical solutions within the puff models.

## Instructions for running code

The calcium puff model based on integrodifferential equations, described in Hawker et al., (2024), is implemented in Matlab. Each model depends on the functions `gatingSolutionMH.m`, `Ca_model.m` and `RK4.m`. By changing the `NumTimes` parameter, one can choose how long they want their calcium puff trace to be. Running the calcium puff model produces `NumTimes` .`mat` files containing data on the calcium concentration, buffer concentration and gating variables for each second. The `.mat` files are compiled in the code to produce a single `.mat` file for the entire run time.

### Parameters

Our model is based on the hybrid stochastic system by Cao et al., (2013). Transition rates for the Markov model are the same as those estimated by Siekmann (2012) and Cao et al., (2013) using single channel stationary data (Wagner and Yule, 2012) and kinetic single channel data (Mak et al., 2007). We choose a ion channel cluster size of 10, this can be changed by altering the `Num_IPR` parameter.

$a_{h42}$ can be considered a tuning parameter within the model and represents the basal concentration of $\lambda_{h42}$. In our manuscript, we keep $a_{h42}$ constant at 0.5s. However, increasing $a_{h42}$ to 5s increases the frequency of calcium puffs (Cao et al., 2013,2014,2017).

$\tau$ represents how much of past calcium concentrations the ion channel has *knowledge* of. We show in Hawker et al., (2024) that if $\tau$ is set to 0.1s, calcium puffs are not produced. A longer integral length (e.g. 3s) is sufficiently long enough to produce calcium puffs. The reduced 2-state model is used to implement the results presented in the section *The effect of* $\tau$ *on calcium dynamics* in Hawker et al., (2024).
