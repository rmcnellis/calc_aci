# Calculate vcmax and vqmax from A/Ci curve
Vcmax is Rubisco-limited photosynthesis
Vqmax is Cytochrome b6f-limited photosnthesis as descibed by Johnson & Berry 2022

## Inputs
- Dataframe including:
  - Photo = Photosynthesis, umol m-2 s-1
  - Ci = Intercellular CO2, umol mol-1
  - Rd = Dark respiration, umol m-2 s-1
  - PAR = Photosynthetically active radiation, umol m-2 s-1
  - Tleaf = Leaf temperature, C
- ci_trans = Ci transition point, estimated from data
- nl = ATP per e- in linear flow, ATP/e- (defaults to 0.75)
- nc = ATP per e- in cyclic flow, ATP/e- (defaults to 1)
- Abs = Total leaf absorptance to PAR, mol PPFD absorbed mol-1 PPFD incident (defaults to 0.85)
- beta = PSII fraction of total leaf absorptance, mol PPFD absorbed by PSII mol-1 PPFD absorbed (defaults to 0.52)
- Kf = Rate constant for fluoresence at PSII and PSI, s-1 (defaults to 0.05e09)
- Kd = Rate constant for constitutive heat loss at PSII and PSI, s-1 (defaults to 0.55e09)
- Kp1 = Rate constant for photochemistry at PSI, s-1 (defaults to 14.5e09)
- start_vmax = a named list or named numeric vector of starting estimates (defaults to 30)
- control_maxiter = positive integer, termination occurs when the number of iterations reaches maxiter (defaults to 500)
- control_minFactor = positive numeric value specifying the minimum step-size factor allowed on any step in the iteration (defaults to 1/10000)

## Citation

Johnson, J. E. and J. A. Berry. 2021. The role of Cytochrome b<sub>6</sub>f in the 
control of steady-state photosynthesis: a conceptual and quantitative model.
 *Photosynthesis Research*, DOI: 10.1007/s11120-021-00840-4
