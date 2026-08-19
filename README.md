SSICE — Integrated Super Script for Enzyme Kinetics
MATLAB scripts for automated enzyme kinetics analysis, multi-model discrimination
(AIC-based), and biocatalytic reactor simulation. This repository accompanies the
manuscript describing the SSICE framework (submitted to Biotechnology and
Bioengineering).
Contents
Script	Description	Manuscript figure
`michaelis_menten_analysis.m`	Interactive Michaelis–Menten analysis: nonlinear fit plus Lineweaver–Burk, Eadie–Hofstee and Hanes–Woolf linearizations, with residual diagnostics.	Figure 1
`ssice_integrated_kinetics.m`	Core SSICE script. Module I extracts initial velocities (v0) from raw time–signal curves with a moving-window linearity algorithm. Module II performs multi-model fitting (Michaelis–Menten, Hill, Haldane, MWC, substrate activation; or competitive / non-competitive / uncompetitive / mixed inhibition) and selects the best model by AIC.	Figure 2
`batch_reactor_product_inhibition.m`	Batch-reactor progress-curve simulation with uncompetitive product inhibition (Euler integration), based on the maltose hydrolysis model of Baş et al. (2007).	Figure 3
`immobilized_catalase_reusability.m`	Cycle-by-cycle reusability simulation of immobilized catalase (chitosan vs. chitosan–bentonite) with first-order activity decay, calibrated with the retentions reported by Kaushal et al. (2018).	Figure 4
Python validation scripts (`python/`)
Reference kinetic models validated against two experimental systems published
by the authors:
Script	Description	Reference system
`python/lipase_orlistat_kinetics.py`	Michaelian enzyme case: Michaelis–Menten fit of porcine pancreatic lipase, global competitive-inhibition fit, Lineweaver–Burk and Dixon plots (Km, Vmax, Ki).	Lipase–orlistat (Candela et al., J. Chem. Educ. 2021)
`python/catalase_resorcinarene_kinetics.py`	Non-Michaelian enzyme case: first-order kinetics of bovine liver catalase (no saturation), non-competitive inhibition model, double-reciprocal diagnosis and IC50 ranking of five sulfonated resorcinarenes.	Catalase–resorcinarenes (Collazos et al., Int. J. Biol. Macromol. 2019)
Run them with:
```bash
cd python
pip install -r requirements.txt
python3 lipase_orlistat_kinetics.py
python3 catalase_resorcinarene_kinetics.py
```
Figures are saved automatically to the `figures/` folder.
Requirements
MATLAB R2018b or later recommended.
`michaelis_menten_analysis.m` and `ssice_integrated_kinetics.m` require the
Optimization Toolbox (`lsqcurvefit`).
Scripts 3 and 4 run on a base MATLAB installation
(`immobilized_catalase_reusability.m` uses `yline`, R2018b+).
Usage
Clone the repository and run any script from the MATLAB command window:
```matlab
>> michaelis_menten_analysis   % interactive menu (example/manual/file/simulation)
>> ssice_integrated_kinetics   % runs with the built-in example dataset
>> batch_reactor_product_inhibition
>> immobilized_catalase_reusability
```
To use your own data with the SSICE core script, replace the `macro_data` matrix
(`[S], [I], v` columns) in `ssice_integrated_kinetics.m`, or feed it with the v0
values produced by Module I. In `michaelis_menten_analysis.m`, option 3 of the
menu loads any plain-text file with two columns: `[S]` and `v`.
References
Baş, D., Duduk, F. C., & Boyacı, İ. H. (2007). Modeling and optimization IV:
Investigation of reaction kinetics and kinetic parameters using artificial
neural networks. International Journal of Food Engineering, 3(1).
Kaushal, J., Singh, G., & Arya, S. K. (2018). Immobilization of catalase onto
chitosan and chitosan–bentonite complex: A comparative study. Biotechnology
Reports, 18, e00258. https://doi.org/10.1016/j.btre.2018.e00258
Candela, M. F., Arenas, N. E., Caicedo, O., & Malagón, A. (2021). In silico
inhibition of lipase by orlistat: Kinetics combined with approaches to
visualize interactions. Journal of Chemical Education.
https://doi.org/10.1021/acs.jchemed.0c01184
Collazos, N., García, G., Malagón, A., Caicedo, O., & Vargas, E. F. (2019).
Binding interactions of a series of sulfonated water-soluble resorcinarenes
with bovine liver catalase. International Journal of Biological
Macromolecules, 139, 75–84. https://doi.org/10.1016/j.ijbiomac.2019.07.197
Citation
If you use these scripts, please cite the associated manuscript
(citation details will be added upon publication).
License
<!-- Choose a license before making the repository public, e.g. MIT:
