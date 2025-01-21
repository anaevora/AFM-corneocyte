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
