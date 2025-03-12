## Version 1.4.7
  - added CITATION file
  - cleaned up C++ code
## Version 1.4
  - Cleaner README and Vignettes
  - Extend support to M1 processors where sizeof(long double) < 16
  - Comply with _R_CHECK_LENGTH_0_LOGIC2_ 
## Version 1.3
  - GUILDS is now on GitHub: https://github.com/thijsjanzen/GUILDS
  - Wrote code tests to check code integrity, code coverage is >95
  - Modified maximum likelihood functions to take into account theta_x = theta_y = theta / 2
  - added a function to plot preston-style plots
## Version 1.2.1
  - Updated the User manual
## Version 1.2
  - fixed memory leak issues by adding extra vector access checks
  - fixed memory leak issues by introducing vectors in KDA code
  - renamed logLik to avoid shadowing of the function logLik in the package stats
## Version 1.1
  - removed malloc header from KDA code
