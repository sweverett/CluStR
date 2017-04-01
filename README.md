# CluStR

This package calculates various scaling relations from cluster catalogs. While only appropriate for fitting power laws, it uses a methodology developed by Brandon Kelly for a fully Bayesian approach that incorporates correlated and heteroscedastic measurement errors, intrinsic scatter, selection effects, censored data, and gaussian mixture modeling for the covariates. Alternatively, a different methodology developed by Adam Mantz can be used that handles multivariate regression and replaces the gaussian mixture model with a Dirichlet process - at the expense of no censored data handling. If desired, both methods can be used and compared.

For more information on the Kelly method, read his [2007 paper](https://arxiv.org/pdf/0705.2774.pdf), check out his original IDL implementation [linmix_err](https://idlastro.gsfc.nasa.gov/ftp/pro/math/linmix_err.pro), or see the python port [linmix](https://github.com/jmeyers314/linmix) by Josh Meyers, which was used in this work.

For more information on the Mantz method, read his [2016 paper](https://arxiv.org/pdf/1509.00908.pdf) or check out his R implementation [lrgs](https://github.com/abmantz/lrgs).

## Getting Started

First you will need to clone the repo by moving to the desired local directory and using

```
git clone https://github.com/sweverett/CluStR
```


### Dependencies

Besides some standard packages like numpy, scipy, and matplotlib that can be aquired through common distributions or pip, CluStR requires the following python packages be installed:

* [astropy](http://www.astropy.org/)
* [corner]()
* [PyPDF2]()
* [rpy2]()
* [linmix]()


Note that astropy is now included in Anaconda.


You can look at these links for details, or simply use the following:

```
pip install astropy
pip install corner
```

CluStR also requires R and the following R package:

* [lrgs](https://github.com/abmantz/lrgs)

This is needed as CluStR interfaces with the original lrgs code through rpy2. The easiest way to do this is to use the Anaconda distribution and then use the following:

```
conda install -c r r-essentials
conda update -c r r-essentials
```

This installs R and some essential packages that works well with Anaconda and Jupyter. Once installed, open R in the terminal and input the following commands:

```
install.packages("devtools")
devtools::install_github("hoxo-m/githubinstall")
library(devtools)
install_github("abmantz/lrgs")
```

## Config File

## Example Use

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Add Tesla, Devon, etc.
