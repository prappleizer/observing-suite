# observing-suite
A suite of tools for planning and carrying out astronomical observations, building on astropy and inspired by astroplan.

This package is built around three core classes: `Target()`, where one specifies astronomical targets and various configurations of desired target, `ObservingPlan`, which ingests a target list, sets up the observatory and nights of observation, and creates beautiful html observing plans (with airmass plots, finder charts, etc) and exports observatory-ready target lists for all unique configurations, and `ObservingLog`, which spawns an easy to use, but highly flexible, html entry tool that saves an observing log during the night. 

Numerous convenience functions exist at each stage: my goal is to make it as easy as possible for me (and others) to set up complex observing runs yet keep track of everything. Some elements of this package are made better with instrument specific and observatory specific information. For now, Keck Observatory and Palomar Observatory are the two "most supported". I'll happily add more if people are interested and can provide, e.g., the template for that observatory's targetlist. 
