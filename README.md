# Delay_differential_equation_HostParasitoid_models_LTCC
Matlab code to simulate delay differential equation, stage-structured, Host-Parasitoid models under long-term environmental change.

Temperature drives changes in Host & Parasitoid larval development time,
reduce Parasitoid attack rate and increase Parasitoid larval mortality 
rate, all based on empirical measurements.

Code originally written by Christophe Coste <christophe.coste@nina.no>, 
<c.f.d.coste.ecology@gmail.com>, funded by NERC grant number
NE/W006731/1 awarded to Mike S Fowler, Steve M Sait and Chenggui Yuan.
Developed from a previous stochastic HP model version by David Gilljam
<david.gilljam@slu.se>, funded by NERC grant number NE/N002849/1 awarded 
to Mike S Fowler, Steve M Sait and Chenggui Yuan.
This version cleaned up by by Mike S Fowler <m.s.fowler@swansea.ac.uk> 
03 Nov 2026.
Runs successfully on Matlab R2023b (v23.2) and R2025a on 64-bit Mac Apple 
M3 Pro. Doesn't require additional Toolboxes.

N.B. Fix required: Parasitoid Adult densities should be reset to 0 after
declining below ymin (e.g., ymin = 1e-6), but remain at positive value
just below this pseudo-extinction threshold.
