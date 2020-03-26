<b>Author:</b> Payman Tohidifar  
<b>Corresponding author:</b> Christopher V. Rao  
Department of Chemical & Biomolecular Engineering  
University of Illinois at Urbana-Champaugn  

<b>Motivation:</b>  
Motivation was to quantify true concentration of ethanol that cells experience in the pond near the mouth of capillary.  
See fem_cap3d.pdf for formulation of fem problem.


Required tools and python libraries:  
---
>Gmsh (v-4.4.2)  
fenics (v-2019.1.0) (prebuilt Anaconda Python packages were installed in fenicsproject environment)  
Python (v-3.6)  
numpy (v-1.16.4)  
matplotlib (v-3.1.1)  


Usage:  
---
Mesh geometry was created as a cap3d.geo file.  
Gmsh was used to create cap3d.msh file:  
`gmsh -3 cap3d.geo`  
fenics built-in meshio was used to generate cap3d.xml file from cap3d.msh file:  
`meshio-convert cap3d.msh cap3d.xml`  
Run `conda activate fenicsproject` to be able use fenics libraries.  
Run `python cap_assay.py` to solve finite-element problem.



