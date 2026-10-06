# Sterile Neutrino Constraints

This repository lists the current terrestrial, astrophysical and cosmological constraints and future sensitivities on the **active-sterile mixing** between the three active neutrino flavours of the Standard Model ($\nu_e$, $\nu_\mu$, and $\nu_\tau$) and a single sterile neutrino flavour ($\nu_s$). The active and sterile states mix to yield the mass eigenstates $\nu_i$ as

$$ \nu_\ell = \sum_i^{3+1} U_{\ell i} \nu_i ~,$$

where $\ell = e,\mu,\tau,s$. The state $\nu_4$ is often called a **sterile neutrino** (despite having a small active component) or **heavy neutral lepton** ($N$) for $m_4 = m_N \gtrsim 1~\text{MeV}$. Here, limits on $|U_{eN}|^2$, $|U_{\mu N}|^2$ and $|U_{\tau N}|^2$ are given as a function of the sterile neutrino mass in the range $10~\mu \text{eV} < m_N < 10~\text{TeV}$. These constraints were originally collected for the paper [1912.03058](https://arxiv.org/abs/1912.03058) and hosted at this [webpage](https://www.hep.ucl.ac.uk/~pbolton/). This repository is an updated database and makes corrections where errors were spotted.

**Model assumptions**: We consider the phenomenological scenario of a single Majorana $\nu_4/N$ with a single flavour of active-sterile mixing. This setup cannot generate masses for two of the light neutrinos, and therefore is incompatible with the neutrino oscillation data. It nevertheless provides a useful benchmark for experimental searches, and the bounds can be scaled appropriately in the context of a specific neutrino mass model.

For the bounds from active-sterile neutrino oscillation searches, i.e $\nu_\alpha \to \nu_s$, we make use of the results assuming a normal ordering of neutrino masses, with the lightest neutrino being massless, $m_1 = 0$. Thus, the sterile neutrino mass is determined as $m_N = \sqrt{\Delta m_{41}^2}$ and the active-sterile mixing strength as

$$|U_{\alpha N}|^2 \approx \sin^2\theta_{\alpha 4} = \frac{1- \sqrt{1-\sin^22\theta_{\alpha 4}}}{2} \approx \frac{\sin^22\theta_{\alpha 4}}{4}$$

where $\theta_{e 4} = \theta_{14}$, $\theta_{\mu 4} = \theta_{24}$ and $\theta_{\tau 4} = \theta_{34}$.

For any queries, comments or additional constraints that should be included here, contact me at patrick.d.bolton@gmail.com.

# Mixing Landscape

<p align="center"><img src="plots/summary_constraints.png" width="700"></p>

# Electron Mixing

<!-- ### [**Electron mixing**](docs/electron_mixing) -->

Current constraints ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_constraints.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_constraints.png))\
Future sensitivities ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_future.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_future.png))

[<p align="center"><img src="plots/UeNsq_constraints.png" width="700"></p>](docs/electron_mixing.md)

[<p align="center"><img align="center" src="plots/UeNsq_future.png" width="700"></p>](docs/electron_mixing.md)

# Muon Mixing

<!-- ### [**Muon mixing**](docs/muon_mixing) -->

Current constraints ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_constraints.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_constraints.png))\
Future sensitivities ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_future.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_future.png))

[<p align="center"><img src="plots/UmuNsq_constraints.png" width="700"></p>](docs/muon_mixing.md)

[<p align="center"><img align="center" src="plots/UmuNsq_future.png" width="700"></p>](docs/muon_mixing.md)

# Tau Mixing

<!-- ### [**Tau mixing**](docs/tau_mixing) -->

Current constraints ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_constraints.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_constraints.png))\
Future sensitivities ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_future.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_future.png))

[<p align="center"><img src="plots/UtauNsq_constraints.png" width="700"></p>](docs/tau_mixing.md)

[<p align="center"><img align="center" src="plots/UtauNsq_future.png" width="700"></p>](docs/tau_mixing.md)