## BCUH

BCUH provides statistical methods to post-process climate model simulations for impact study purposes. A tutorial is available [here](http://ratyo.github.io/BCUH)

## Features

 - 10 univariate methods designed for (daily mean) temperature
 - 10 univariate methods designed for (daily) precipitation
 - copula-based bias correction method for the joint distribution of temperature and precipitation
 
## Installation

```s
devtools:install_github("Ratyo/BCUH")
```

## Usage

There are two main wrapper functions. 

To apply univariate bias correction use

```r
  BCUH::biasco()
```

For using the joint bias correction method use

```r
  BCUH::biasco2D()
```


A more detailed is given in the package's [documentation](http://ratyo.github.io/BCUH)

