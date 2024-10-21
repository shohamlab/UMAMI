The UMAMI demdulation procedure works on single frames at a time. Breifly, first thing done,
and optionally, is motion corretions with NoRMCorre[1]. Second is the filtering procedure
which is included as a function in demodulationProceudre.m.

For custom use and depending on scan and iamge parameters, the values for the filter location and width may need
to be adjusted.




[1] Pnevmatikakis, E.A. & Giovannucci, A. NoRMCorre: An online algorithm for piecewise rigid
motion correction of calcium imaging data. Journal of Neuroscience Methods 291, 83-94 (2017).
