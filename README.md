# CluStR

This package calculates various scaling relations from cluster catalogs. While only appropriate for fitting power laws, it uses a methodology developed by Brandon Kelly for a fully Bayesian approach that incorporates correlated and heteroscedastic measurement errors, intrinsic scatter, selection effects, censored data, and gaussian mixture modeling for the covariates. Alternatively, a different methodology developed by Adam Mantz can be used that handles multivariate regression and replaces the gaussian mixture model with a Dirichlet process - at the expense of no censored data handling. If desired, both methods can be used and compared.

For more information on the Kelly method, read his [2007 paper](https://arxiv.org/pdf/0705.2774.pdf), check out his original IDL implementation [linmix_err](https://idlastro.gsfc.nasa.gov/ftp/pro/math/linmix_err.pro), or see the python port [linmix](https://github.com/jmeyers314/linmix) by Josh Meyers, which was used in this work.

For more information on the Mantz method, read his [2016 paper](https://arxiv.org/pdf/1509.00908.pdf) or check out his R implementation [lrgs](https://github.com/abmantz/lrgs).

## Getting Started (using `conda env`)

First, you will need [conda](https://conda.io/docs/intro.html).

Once you have that, clone this repository into a local directory, and switch to that directory:

```bash
git clone https://github.com/sweverett/CluStR
cd CluStR
```

Next, create a `conda` environment using the supplied `environment.yml`:

```bash
conda env create -f environment.yml
```

Now activate the newly-created `conda` environment:
```bash
source activate clustr
```

Next, we need to install the necessary `R` packages. Run `R`, and then enter:
```R
install.packages("devtools")
install_github("abmantz/lrgs")
```

**Note on Mac OS**: On OSX Mavericks or more recent, you need to run `R` as:
```bash
TAR=“/usr/bin/tar” R
```

Now you should be ready to use `CluStR`! Whenever you want to run `CluStR`, activate the `conda` environment with:
```bash
source activate clustr
```

Whenever you are finished, deactivate the `conda` environment with
```bash
source deactivate clustr
```

## Getting Started (using `pip` and `conda`)

First you will need to clone the repo by moving to the desired local directory and using

```
git clone https://github.com/sweverett/CluStR
```

### Dependencies

Besides some standard packages like numpy, scipy, and matplotlib that can be aquired through common distributions or pip, CluStR requires the following python packages be installed:

* [astropy](http://www.astropy.org/)
* [corner](http://corner.readthedocs.io/en/latest/)
* [PyPDF2](http://pythonhosted.org/PyPDF2/)
* [linmix](https://github.com/jmeyers314/linmix)
* [rpy2]()


Note that astropy is now included in Anaconda.


You can look at these links for details, or simply paste the following into terminal (don't install rpy2 if you don't have R yet):

```
pip install astropy
pip install corner
pip install pypdf2
```

The simplest way to get `linmix` is to clone the repo and install using the given setup file:

```
git clone https://github.com/jmeyers314/linmix.git
cd linmix
python setup.py install
```

CluStR also requires R and the following R package:

* [lrgs](https://github.com/abmantz/lrgs)

This is needed as CluStR interfaces with the original lrgs code through rpy2. The easiest way to do this is to use the Anaconda distribution and then use the following:

```
conda install -c r r-essentials
conda update -c r r-essentials
```

This installs R and some essential packages that works well with Anaconda and Jupyter. Now that R is installed, grab rpy2.

```
pip install rpy2
```

Once everything is installed, open R in the terminal and input the following commands:

```
install.packages("devtools")
devtools::install_github("hoxo-m/githubinstall")
library(devtools)
install_github("abmantz/lrgs")
```

Now you should be ready to use `CluStR`!

## Config File <a name="config"></a>

Most parameters are set in the `param.config` file. Here you can set the cosmology, default regression method, plotting options, and most importantly any desired flags. There are three possible flag types: bool, cutoff, and range. For each, you must specify the exact catalog column name you want to make cuts along with the flag type and, if a cutoff or range, the corresponding cut values. All name:value pairs must be separated by a colon.

There are two important things to note that might be unclear:

* Setting a flag type and cut value **does not mean the cut will be used!** A flag is set to be used in the actual method call - see [Example Use](#exuse) below. This allows you to set many flag parameters without having to change the config file everytime you want to use a different combination of flags.

* While it may seem counter-intuitive at first, the flag parameters are set to what data you want to *keep*, not remove - i.e. set what redshift range you want your clusters to be found in rather than what ranges you want to remove. I found this to eliminate some mental gynmastics when setting cuts but may feel awkward for bools.

Here is an example for each flag type:

### Bool: 

*<column_name>_bool_type: <True/False>* 

To only include clusters that are not within r500 of a chip edge, use

```
edge_exclude_r500_bool_type: False
```

In other words - only use data that is *not* flagged with `edge_exclude_r500`.

### Cutoff:

*<column_name>_cut_type: <above/below>*

*<column_name>_cut: <value>*

To analyze clusters whose redshift is below 1.0, use

```
redshift_cut_type: below
redshift_cut: 1.0
```

### Range:

*<column_name>_range_type: <inside/outside>*

*<column_name>_range_min: <value>*

*<column_name>_range_max: <value>*

To only use clusters with redshift between 0.3 and 0.5, use

```
redshift_range_type: inside
redshift_range_min: 0.3
redshift_range_max: 0.5
```

This system may seem inefficient, but allows for quite a bit of flexibility in selecting interesting data subsets.

## Example Use <a name="exuse"></a>

*python clustr.py <catalog.fits> <response> <covariate>*

CluStR has three mandatory inputs: An appropriate cluster FITS file, the response variable (y-axis), and the covariate variable (x-axis). The available columns for axis variables (on the right) and their corresponding input labels (on the left) are:

- lambda : lambda (richness)
- l500kpc : 500_kiloparsecs_band_lumin
- lr2500 : r2500_band_lumin
- lr500 : r500_band_lumin
- lr500cc : r500_core_cropped_band_lumin
- t500kpc : 500_kiloparsecs_temperature
- tr2500 : r2500_temperature
- tr500 : r500_temperature
- tr500cc : r500_core_cropped_temperature

To plot the scaling relation between r2500 temperature and richness using the default regression method (set in `param.config`), use

```
python clustr.py <catalog.fits> tr2500 lambda
```

The output file will be named `<default_prefix>r2500_temperature-lambda.pdf`, where you can set the default prefix in `param.config`.

Additionally there are three optional arguments: A filename prefix (`-p`), regression method (`-m`), and flag options (`-f`). The order of the first two is arbitrary, but all desired flags must be specified last. As described in the [Config File](#config) section, flag paramters are set in `param.config` but are only used if present in the method call following `-f`.

To plot the scaling relation between r2500 temperature and richness using both regression methods but only including clusters with redshift between 0.3 and 0.5, use

```
python clustr.py <catalog.fits> tr2500 lambda -p SDSS_redshifts_ -m both -f redshift
```

The output file will be named `SDSS_redshifts_r2500_band_lumin-lambda.pdf`

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
