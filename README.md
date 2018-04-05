## BCUH

BCUH provides statistical methods to post-process climate model simulations for impact study purposes. A tutorial is available [here](https://RatyO.github.io/BCUH)

## Features

 - 10 univariate methods designed for (daily mean) temperature
 - 10 univariate methods designed for (daily) precipitation
 - copula-based bias correction method for the joint distribution of temperature and precipitation
 
## Installation

```s
devtools::install_github("Ratyo/BCUH")
```

The package depends on the following R packages:

moments,
MASS,
fitdistrplus,
copula,
methods,
Rdpack,
lubridate

## Usage

There are two main wrapper functions available. 

To apply univariate bias correction use

```r
  BCUH::biasco()
```

For copula-based joint bias correction method, use

```r
  BCUH::biasco2D()
```


A more detailed description is availbale the package's [documentation](https://RatyO.github.io/BCUH)

