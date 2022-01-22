# Bone-Remodelling
Mathematical modelling of bone remodelling cycles including the NFkB signalling pathway


Bone diseases become heavily prevalent in older communities (after 50 years of age). Modelling and simulation of bone remodelling processes substantiates our understanding of said diseases. Bone remodelling involves the RANKL-RANK protein binding in osteoclasts and its subsequent transmission to the nucleous, where osteoclastogenesis is promoted.

Ji et al. (https://doi.org/10.1016/j.compbiomed.2019.03.003) presented a quantitative bone remodelling model in normal bone microenvironment from which bone volume and cell variation can be obtained, as well as biochemical factors involved in the aforementioned transmission. The work presented in this repository replicates this model using MATLAB, and presents satisfactory results which enable the analysis of RANKL impact on bone remodelling.

NFK.m -> function for the NFkB mechanism.

NFK_bone.m -> NFK function with bone remodelling additions.

NFK_bone_RANKL.m -> function enabling user-input simulation of RANKL protein inhibition/injection.

parametros.m -> Script with parameters and constants.

Euler.m -> Script that "solves" the model. 


Run "Euler.m" to obtain the desired reuslts/plots. Includes an user-input binary-coded RANKL inhibition/injection.

Contact me for the full work report (only in portuguese) and any additional questions.
