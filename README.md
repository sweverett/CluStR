# CluStR

This package calculates various scaling relations from cluster catalogs. While only appropriate for fitting power laws, it uses a methodology developed by Brandon Kelly for a fully Bayesian approach that incorporates correlated and heteroscedastic measurement errors, intrinsic scatter, selection effects, censored data, and gaussian mixture modeling for the covariates. Alternatively, a different methodology developed by Adam Mantz can be used that handles multivariate regression and replaces the gaussian mixture model with a Dirichlet process - at the expense of no censored data handling. If desired, both methods can be used and compared.

For more information on the Kelly method, read his [2007 paper](https://arxiv.org/pdf/0705.2774.pdf), check out his original IDL implementation [linmix_err](https://idlastro.gsfc.nasa.gov/ftp/pro/math/linmix_err.pro), or see the python port [linmix](https://github.com/jmeyers314/linmix) by Josh Meyers, which was used in this work.

For more information on the Mantz method, read his [2016 paper](https://arxiv.org/pdf/1509.00908.pdf) or check out his R implementation [lrgs](https://github.com/abmantz/lrgs).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them



### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be


And repeat


End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
