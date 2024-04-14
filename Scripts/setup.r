#libraries

#stan
require(rstan)
require(rstanarm)

#dirichlet
require(gtools)

#skew-normal
require(sn)

#tidyverse
require(tidyverse)
require(bayesplot)
require(reshape2)

#editing, auxiliary libraries
require(flextable)
require(equatags) #for math in flextable
require(beepr)
require(stringr)

#images and plots
require(png)
require(grid)
require(gridExtra)
require(ggrgl)
require(plotly)
require(latex2exp)

#time series tests
require(tseries)
require(strucchange)

#set global knitr options
knitr::opts_chunk$set(
	echo = FALSE,
	message = FALSE,
	warning = FALSE,
	digits = 0
)

#global parameters
min.age=0
max.age=100
sex="Maschi"

set.seed(13)