# Modeling of modulated exosome release from differentiated induced neural stem cells for the treatment of brain cancer
This repository provides code for the analyses presented in the manuscript submitted to IEEE Transactions on NanoBioscience. 

The main m.file for executing simulations on modulated exocytosis in neurons (with the focus on released exosomes) is controlledExocytosisNeurons.m. The execution requires a definition of the so-called control signal – induced current into a cell. Here, we consider a sequence of rectangular pulses. The parameters that should be set comprise the number of pulses, time range (ms), offset (ms) and amplitude (\muA/cm2). 

The main m.file for executing simulations on modulated exocytosis in astrocytes (with the focus on released exosomes) is controlledExocytosisAstrocytes.m. The execution requires a definition of the so-called control signal – induced voltage into a cell. Here, we consider a sequence of rectangular pulses. The parameters that should be set comprise the number of pulses, time range (ms), offset (ms) and amplitude (mV). The user is also instructed to set an initial value as of IP3 concentration. 
