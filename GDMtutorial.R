GDM Workshop
========================================================

```{r}
gitpath<-"C:/Users/Jorge/Documents/GDMWorkshop/"
setwd(paste(gitpath,"Code_and_data",sep=""))

require(raster)
require(ecodist)
```

generate GDM input.r
---------------------
To add the environment values for each site pair

```{r fig.width=7, fig.height=6}
source("generate_site_pairs.r")
```

use the script generate GDM input.r to add the environment values for each site pair
--------------------------------

setwd(paste(gitpath,"Code_and_data",sep=""))
```{r}
source("generate_GDM_input.r")
```

2 Model fitting and testing
2.1 Fit a GDM model
====================
  
```{r}
source(paste(gitpath,"gdm_fit_and_test.r"))
```

