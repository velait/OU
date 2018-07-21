## Hierarchical Ornstein-Uhlenbeck Process in RStan

This repository contains material for the [Helsinki StanCon 2018](http://mc-stan.org/events/stancon2018Helsinki/) presentation by Ville Laitinen & Leo Lahti, [Open Research Labs](http://openresearchlabs.github.io), University of Turku, Finland. The work provides an implementation and initial analysis of a hierarchical Ornstein-Uhlenbeck model by extending the earlier work by [Aaron Goodman, 2018](https://github.com/stan-dev/stancon_talks/tree/master/2018/Contributed-Talks/05_goodman). For further details, see the reproducible [markdown document](https://github.com/velait/OU/blob/master/hierarchical_OU.main.md).


### Replicating the analysis

The material is available under the CC-4.0 license.

Running the entire script from scratch takes a long time (hours on a basic laptop). For this reason the most time consuming parts in OU.source.R have been
commented. Some necessary source files that take long times to generate can be readily [downloaded](https://drive.google.com/drive/folders/15kd6Y2CgoXEH6y0mqEU4YRLR4i_fTTrk?usp=sharing). If you wish to generate these files yourself
simply uncomment rows 38-43, 64-65, 106-115 in [source.R].

To generate the report in html format, unzip the downloaded files and run the script hierarchical_OU.main.Rmd as follows:

```
library(rmarkdown)
render("hierarchical_OU.main.Rmd")
```


### Citing the work

Kindly cite the work as follows: Ville Laitinen and Leo Lahti. A
Hierarchical extension to Ornstein-Uhlenbeck-type Student’s
t-processes. StanCon 2018, Helsinki. URL:
https://github.com/velait/OU/