# MCnebula

MCnebula algorithm integration in R.

## Installation 

The following installation has been tested in Pop!_OS 20.04 (Ubuntu 20.04).

### Install dependencies

Install R. We have tested MCnebula in R version 4.2 and 4.1. Maybe R version
$\geq$ 4.0 is feasible.
Herein, the codes are giving an example with installing R 4.2.

```{bash}
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt-get update
sudo apt install --no-install-recommends r-base
```

