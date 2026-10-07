# Sterile Neutrino Constraints

This repository lists the current terrestrial, astrophysical and cosmological constraints and future sensitivities on the **active-sterile mixing** between the three active neutrino flavours of the Standard Model ($\nu_e$, $\nu_\mu$, and $\nu_\tau$) and a single sterile neutrino. The active neutrinos are composed of mass eigenstates as

$$ \nu_\alpha = \sum_i^3 U_{\alpha i} \nu_i + U_{\alpha N} N,$$

where $\alpha = e,\mu,\tau$. The state $N$ is often called a **sterile neutrino** (despite having a small active component) or **heavy neutral lepton** for $m_N \gtrsim 1~\text{MeV}$. Here, limits on $|U_{eN}|^2$, $|U_{\mu N}|^2$ and $|U_{\tau N}|^2$ are given as a function of the sterile neutrino mass in the range $10~\mu \text{eV} < m_N < 10~\text{PeV}$. These constraints were originally collected for the paper [1912.03058](https://arxiv.org/abs/1912.03058) and hosted at this [webpage](https://www.hep.ucl.ac.uk/~pbolton/). This repository is an updated database and corrects previous errors.

For any queries, comments or additional constraints that you think should be included here, contact me at patrick.d.bolton@gmail.com.

## Model Assumptions

We consider the phenomenological scenario of a massive **Majorana** fermion $N$ with a **single flavour** active-sterile mixing. This setup cannot generate masses for two of the light neutrinos and therefore is incompatible with neutrino oscillation data. Nevertheless, it provides a useful benchmark for experimental searches, and the bounds can be rescaled appropriately in the context of a specific neutrino mass model. 

To plot some of the constraints below, we make the additional assumptions:

**Seesaw line**:

$$|U_{\alpha N}|^2 = \frac{m_3}{m_3 + m_N},$$

with $m_3 = 0.05~\text{eV}$.

**Active-sterile neutrino oscillations**: Bounds from active-sterile oscillation searches are often given in terms of $\Delta m_{41}^2 = m_4^2 - m_1^2$ and $\sin^2 2\theta_{i4}$ for $i = 1,2,3$. To plot these bounds, we assume the normal ordering of neutrino masses ($m_1 < m_2 < m_3$) and take the lightest neutrino to be massless, $m_1 = 0$. Thus, the sterile neutrino mass and active-sterile mixing are determined as

$$m_N = \sqrt{\Delta m_{41}^2} , \quad |U_{\alpha N}|^2 = \sin^2\theta_{\alpha 4} = \frac{1}{2}- \frac{1}{2}\sqrt{1-\sin^22\theta_{\alpha 4}} \approx \frac{\sin^22\theta_{\alpha 4}}{4} ,$$

respectively, where $\theta_{e 4} = \theta_{14}$, $\theta_{\mu 4} = \theta_{24}$ and $\theta_{\tau 4} = \theta_{34}$. Note that $m_N > m_{2,3}$ is not strictly required, and we consider $m_N < m_{2,3}$ in the plots below.

**Neutrinoless double beta decay**:

$$(T_{1/2}^{0\nu})^{-1} = \frac{G_{01} |V_{ud}|^4 g_A^4 }{m_e^2} \Big|\sum_i^3 U_{ei}^2 m_i \mathcal{A}_\nu + U_{eN}^2 m_N \mathcal{A}_N(m_N)\Big|^2 ,$$

$$ \mathcal{A}_N(m_N) = \begin{cases}
\mathcal{A}(m_N) + m_N \frac{d}{dm_N}\mathcal{M}(m_N) + \mathcal{A}_\nu^{\text{usoft}}(m_N) & m_N < 0.1~\text{GeV} \\
\mathcal{A}(m_N) & 0.1~\text{GeV} \leq m_N < 2~\text{GeV}\\
\mathcal{A}^{(9)}(m_N) & m_N \geq 2~\text{GeV}
\end{cases},$$

## Mixing Landscape

Below, the general landscape of current bounds is shown for different active flavours.

<p align="center"><img src="plots/summary_constraints.png" width="700"></p>

## Electron Mixing

<!-- ### [**Electron mixing**](docs/electron_mixing) -->

Current constraints ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_constraints.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_constraints.png))\
Future sensitivities ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_future.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UeNsq_future.png))

[<p align="center"><img src="plots/UeNsq_constraints.png" width="700"></p>](docs/electron_mixing.md)

[<p align="center"><img align="center" src="plots/UeNsq_future.png" width="700"></p>](docs/electron_mixing.md)

## Muon Mixing

<!-- ### [**Muon mixing**](docs/muon_mixing) -->

Current constraints ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_constraints.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_constraints.png))\
Future sensitivities ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_future.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UmuNsq_future.png))

[<p align="center"><img src="plots/UmuNsq_constraints.png" width="700"></p>](docs/muon_mixing.md)

[<p align="center"><img align="center" src="plots/UmuNsq_future.png" width="700"></p>](docs/muon_mixing.md)

## Tau Mixing

<!-- ### [**Tau mixing**](docs/tau_mixing) -->

Current constraints ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_constraints.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_constraints.png))\
Future sensitivities ([pdf](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_future.pdf), [png](https://raw.githubusercontent.com/PatrickBolton/SterileNeutrinoConstraints/main/plots/UtauNsq_future.png))

[<p align="center"><img src="plots/UtauNsq_constraints.png" width="700"></p>](docs/tau_mixing.md)

[<p align="center"><img align="center" src="plots/UtauNsq_future.png" width="700"></p>](docs/tau_mixing.md)