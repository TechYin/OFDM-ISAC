# OFDM-ISAC
This repository is pubic for the optimized pilot-patterns of our work on low-PAPR underwater acoustic OFDM-ISAC system. 

We offered 5 predefined sets of pilot-patterns (also called pilot-pattern dictionaries in the code) after optimization in it, and offered two scripts to test the SVRs and the correlation functions.

The public data is in `PilotPatterns.mat`. There are 5 different cells in this data file, and each of them includes a predefined set of 100 different pilot-patterns.

We also provide a script to test the SVRs of pilot observation matrices (run `SVR_Analysis.m`) and a script to test the P-CCFs between pilot-pattern detection (PD) signals (run `Correlation_Analysis.m`).

We are sorry for some explanatory notes still writing in Chinese, and we will modify them soon in the future version.


Best regards,
The authors.
