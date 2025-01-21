# AFM-corneocyte
This project is used to measure the elastoviscoplastic properties of corneocytes.
Scripts should be run in the following order:

1. "Tip geometry_PDMS_power law"
Nanoindentation force curves performed on PDMS should be used to characterise AFM tip geometry before performing nanoindentation on corneocytes.

This algorithm considers PDMS AFM data to obtain the parameters 'c' and 'n' to describe the geometry of the AFM probe. Young's modulus (E*) was % obtained using Micromanipulation technique.
Note power law F=bh^m where 'b' and 'm' are the coefficients found with this script and used to characterise probe geometry.
These were used to calculate AFM tip geometry parameters as explained in the manuscript (please see main text and SI for detailed derivation aof equations).
The load index, m, of the power law was calculated at different contact depths, hc, based on a polynomial fit of PDMS loading curves. The plane strain elastic modulus, E^*, of the PDMS specimen required to derive the tip geometry was independently calculated from micromanipulation experiments based on the indentation of 3 different regions. 
The derived parameters c and n define the local geometry of the AFM tip at different selected contact depths.

2. "Oliver-Pharr method"
This algorithm uses the Oliver-Pharr method to analyse the retract curves and calculate the stiffness of the sample as the slope of the tangent at maxium applied force.
This script only serves to calculate Young's modulus of cells on tape and is based on the Oliver-Pharr method.

3. "Viscoelastic properties"
This algorithm uses the Oliver-Pharr method to analyse the retract curves and calculate the stiffness of the sample as the slope of the tangent at maximum applied force. Viscoelastic properties are also measured using force-relaxation curves. Parameters mesured with this script are Young's modulus (MPa), Relaxation times	(τ1 and τ2) Hardness (H0 (MPa)	H∞ (MPa)).


4. "Create_vpmodel" and "Viscoplastic model_total"
