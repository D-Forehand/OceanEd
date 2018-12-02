## OceanEd

The OceanEd software is a suite of MATLAB and Simulink software codes for creating and testing linear hydrodynamic time-domain models of arrays of rigid-body Wave Energy Converters (WECs). The resulting time-domain models are fast because they use state-space techniques to approximate the radiation forces. The WECs in the array can have any hull shape, can be arranged in any layout and each WEC can move in up to all six of its degrees of freedom. The resulting time-domain models take into account all the hydrodynamic interactions (i.e. all the radiated and diffracted waves) between all the converters. At the moment, the code uses frequency-domain input data from WAMIT but it is hoped that in the future this will be extended to Nemoh too.

The advantages of a time-domain model (over the frequency-domain model, from which it is derived) are that the time-domain model can simulate transient behaviour and it can also incorporate nonlinear external forces, such as nonlinear Power Take-Off (PTO), nonlinear hydrostatic and nonlinear mooring forces.

The heart of the main "create_td_WEC_array_model.m" code is the derivation of the single radiation state-space model which approximates the complete set of radiation convolution integrals. It is difficult to derive a numerically stable radiation state-space model and this code is designed to attempt to do this. More details on the theory behind, and implementation of, this code can be found in [1].

Copyright (C) 2018, David Forehand, The University of Edinburgh.

Contact: D.Forehand@ed.ac.uk

Repository: https://github.com/D-Forehand/OceanEd

README last modified: 1st December 2018

## MATLAB Toolbox Requirements

The OceanEd software requires the following MATLAB toolboxes:
- MATLAB
- Simulink
- Single Processing Toolbox
- Control System Toolbox

## Repository Structure

The OceanEd repository contains three main folders:
- The "Program_Files" folder: This contains all the source codes.
- The "Examples" folder: This contains all the example folders.
- The "Documentation" folder: This contains documentation on how to run the examples.

The two main codes in the "Program_Files" folder are "create_td_WEC_array_model.m" and "td_v_WAMIT_fd_comparison.m":
- "create_td_WEC_array_model.m" creates the hydrodynamic time-domain model from the WAMIT frequency-domain data.
- "td_v_WAMIT_fd_comparison.m" tests to see if the created time-domain model is working correctly by comparing its output for steady-state oscillatory conditions against the frequency-domain data.

## About

The OceanEd software was created by Dr David Forehand from the School of Engineering in The University of Edinburgh.  The work was supported by the UK's Engineering and Physical Sciences Research Council (EPSRC) funded SuperGen Marine programme.

The author would like to gratefully acknowledge the help from the late Dr Andy McCabe from Lancaster University in creating this code. In particular, the author would like to thank Dr McCabe for the "get_TF.m" function. Andy is very sadly missed.

The author would also like to thank the following colleagues from The University of Edinburgh for their help with this code and the wider wave-to-wire (i.e. from the waves all the way to the electricity network) numerical model:
- Dr Anup Nambiar
- Dr Aristides Kiprakis
- Prof Robin Wallace

## References

[1] Forehand, D.I.M., Kiprakis, A.E., Nambiar, A.J. and Wallace, A.R. (2016) “A fully coupled wave-to-wire model of an array of wave energy converters”. IEEE Transactions on Sustainable Energy. Vol. 7, No. 1, pp. 118-128.




