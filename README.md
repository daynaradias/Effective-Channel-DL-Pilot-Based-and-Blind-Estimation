# Effective-Channel-DL-Pilot-Based-and-Blind-stimation

This code simulates the results of the published papers: 

D. D. Souza, M. M. M. Freitas, G. S. Borges, A. M. Cavalcante, D. B. da Costa, J. C. W. A. Costa,
"Effective Channel Blind Estimation in Cell-Free Massive MIMO Networks," in IEEE Wireless
Communications Letters, vol. 11, no. 3, pp. 468-472, March 2022. [Doi:10.1109/LWC.2021.3132418](https://doi.org/10.1109/LWC.2021.3132418)

D. D. Souza, M. M. M. Freitas, D. B. da Costa, G. S. Borges, A. M. Cavalcante, L. Valcarenghi, and J. C.
W. A. Costa, “Effective channel DL pilot-based estimation in user-centric cell-free massive MIMO
networks,” in Proc. IEEE Global Commun. Conf., pp. 705-710, Dec. 2022. [Doi:10.1109/GLOBECOM48099.2022.10001150](https://doi.org/10.1109/GLOBECOM48099.2022.10001150)

## Instructions

run the **'main_vartiations.m'** code to simulate

This code defines the main parameters of the simulation, as the number of UEs, number of APs, pilot length, clustering method, etc.

The function **'mainFunction.m'** receives these parameters and saves the results in .mat when finished

The function **'generateSetup.m'** generate the channels and AP clustering, change the channel parameters like carrier frequency and bandwidth inside it.

The function **'functionComputeSE_downlink_sCSI_pCSI_BE_DLPE.m'** is responsible for performing the estimation of the effective channel, using **_DL pilot-based_** method and **_Blind estimation_**. Different methods for DL pilot assignment is also performed. Different approachs for computing the espectral effeciency can also be selected.
