# Thermal Perfromance Curve Modeling
Files associated with the Thermal Performance Curve model analysis of Fusarium fungal growth data. This is part of a larger Masters project.


Fusarium fungi were exposed to temperatures ranging from 10 to 40°C to examine their fungal growth rate performance under different thermal conditions. Ten functions aimed at estimating thermal optima were fitted to radial growth rate data using the “rTPC” package (Padfield et al., 2025). 

Models were selected specifically for their suitability in modeling fungal growth but this might be different for other projects and depends on the question being answered. 

More infromation can be found at https://padpadpadpad.github.io/rTPC/articles/rTPC.html 

The best fit thermal performance curve model was chosen by considering the small sample size Akaike 
Information Criterion (AICc), Bayesian Information Criterion (BIC) scores and adjusted R-squared (goodness of fit). Through non-parametric bootstrapping (1000 resamples with replacements) was used to generate 95% confidence intervals for thermal performance parameters, specifically for thermal optimum (Topt) and maximum performance (Rmax) for each timepoint. Thermal parameters of Topt and Rmax were considered significant when 
their 95% confidence intervals did not overlap across temperatures or days considered.

