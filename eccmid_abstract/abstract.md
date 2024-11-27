Background

The bystander effect of antimicrobials can be defined as the process by which
they exert selection pressure on microbes that are not the intended treatment
target. Current antimicrobial stewardship protocols aim to minimise it by
restricting antimicrobial spectrum or duration, but the optimal approach to
measure the bystander effect of different stewardship strategies is unknown. We
present here a novel sequencing/modelling approach to address this data gap to
understand and quantify the bystander effect of different antimicrobials on the
human gut resistome and microbiome in a low-resource setting.

Methods

We longitudinally sampled up to 4 stool or rectal swabs from each of 163
Malawian adults in Blantye, Malawi, over 6 months; antibiotic exposed/unexposed
hospital inpatients (n=109, n=29, respectively) and community members (n=25),
followed by shotgun metagenomic sequencing. Resistome (defined using
AMRfinderplus) and microbiome composition (using Kraken2) were related to
antibiotic exposure using Bayesian logistic or negative binomial regression with
a correlation structure to account for repeated sampling, and a time-dependent
antibiotic effect that exponentially decayed over time, fit with Stan.

Results 

Ceftriaxone (received by 87% of the antibiotic exposed), cotrimoxazole (54%)
ciprofloxacin (28%), and amoxicillin (17%) were common exposures. Carriage of
cephalosporin resistance determinants was almost universal, and further
increased folloing ceftriaxone exposure. Ceftriaxone was also associated with a
dramatic increase in aminoglycoside and clindamycin/erthromycin resistance
genes, changes not seen with cotrimoxazole or amoxicillin (Figure 1).
Ceftriaxone was also associated with profound changes in microbiome composition
(Figure 2); an increase in Proteobacteria, particularly Enterobacterales, and
Escherichia within this order. The effect of amoxcillin was less marked, and
Cotrimoxazole did not demonstrate this effect.

Conclusions

Ceftriaxone acts to profoundly alter microbiome composition and exert a stromg
bystander effect, enriching for aminoglycoside and clindamycin/erythromycin
resistance genes in this setting. Our flexible modelling approach allows us to
quantify this effect. Simulation from fitted models can rapidly test *in silico*
different stewardship interventions to identify optimal ways to minimise
selection pressures, which can be taken forward to inteventional trials.
