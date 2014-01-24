GDM Workshop
========================================================


```r
gitpath <- "C:/Users/Jorge/Documents/GDMWorkshop/"
setwd(paste(gitpath, "Code_and_data", sep = ""))

require(raster)
```

```
## Loading required package: raster Loading required package: sp
```

```r
require(ecodist)
```

```
## Loading required package: ecodist
```

```
## Warning: package 'ecodist' was built under R version 3.0.2
```

```
## Attaching package: 'ecodist'
## 
## The following object is masked from 'package:raster':
## 
## crosstab, distance
```


generate GDM input.r
---------------------
To add the environment values for each site pair


```r
source("generate_site_pairs.r")
```

```
## Warning: cannot open file 'generate_site_pairs.r': No such file or
## directory
```

```
## Error: cannot open the connection
```


use the script generate GDM input.r to add the environment values for each site pair
--------------------------------

setwd(paste(gitpath,"Code_and_data",sep=""))

```r
source("generate_GDM_input.r")
```

```
## Warning: cannot open file 'generate_GDM_input.r': No such file or
## directory
```

```
## Error: cannot open the connection
```


2 Model fitting and testing
2.1 Fit a GDM model
====================
  

```r
source("gdm_fit_and_test.r")
```

```
## Warning: cannot open file 'gdm_fit_and_test.r': No such file or directory
```

```
## Error: cannot open the connection
```


