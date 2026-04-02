# timmins_project

This project is used to determine how vapor pressure deficit impacts fungal loads,

tissue invasion, and survival of little brown bats with white-nose syndrome over time.

I am making changes to test commit and push

# week 2

upload paragraph from first assignment:

Maybe I can use the time this semester to explore VPD in other species?

This dataset is a collection of infection data and climate data accumulated over the past

decade or so. My main variables of interest are lgdL or fungal loads, UV scores or

tissue invasion, logVPD or log transformed vapor pressure deficit, and avgTEMP or

hierarchical temperature averages. I would like to explore the relationship between

vapor pressure deficit, temperature, and different measures of disease.

# week 3

I made three figures:

Average VPD from early hiber effects on late hiber UV scores

Average VPD effects (from heirarichal model) on late hiber UV scores (as a comparison)

Average VPD from early hiber effects on late hiber fungal loads

Data: mylu.working.extrapolated

# week 5

code: week5_tests

data: mylu.working.extrapolated

I completed a permutation test to assess if avglogVPD is different between WI and NY.
I completed a shapiro-wilks test to assess if avglogVPD is normally distributed.


# week 7 & 8

Week 7

code: week7_lm

data: mylu.working.extrapolated

I created a linear model with my predictor variable being VPD values from 
early hibernation, and my y variable being fungal loads in late hibernation. I ran
model diagnostics on this and created a figure. I also looked at this in an 
early hibernation model. Conclusions seem to be that nothing is normal and
I may need to transform VPD differently.

Week 8

code: week7_lm

data: mylu.working.extrapolated

Additive model interpretation:
In both the established and invasion phases, early hibernation VPD has a statistically significant 
negative impact on fungal loads in late hibernation. In other words, drier early hibernation
conditions are correlated with lower fungal loads in late hibernation.

Fungal loads are not statistically different in the invasion versus established phase,
although interestingly the trend is that fungal loads are higher in the established phase.
This could be due to survival bias.

Interactive model interpretation:
In both the established and invasion phases, early hibernation VPD has a statistically 
significant negative impact on fungal loads in late hibernation. 
This trend is much steeper in the invasion phase, and the difference between 
the trend in established versus invasion phase is statistically significant. 
There is significantly higher fungal loads in invasion phase 
than in the established phase at low values of VPD.

# Week 10

code: 

data: mylu.working.imputated

Results section:

In the epidemic phase of the pathogen, we found that fungal loads are 
significantly higher in warmer environmental temperatures in late hibernation 
(environmental temperature coefficient β = 0.16429 ± 0.05384 SE, p = 0.00233).

In the epidemic phase of the pathogen, we found that fungal loads are 
significantly lower in dry environmental temperatures in late hibernation 
(environmental vapor pressure deficit coefficient β = -0.40213 ± 0.14640 SE, p = 0.00612).

The interaction between environmental temperature and vapor pressure deficit has a
statistically unclear effect on fungal loads, although the trend is positive
(interaction coefficient β = 0.04504  ± 0.02924 SE, p = 0.12380). 
