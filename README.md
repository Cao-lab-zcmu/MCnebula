# MCnebula

&ensp;&ensp; 
MCnebula algorithm integration in R.

## Installation 

&ensp;&ensp; 
The following installation has been tested in Pop!_OS 20.04 (Ubuntu 20.04).

### Install R

&ensp;&ensp; 
Install R. We have tested MCnebula in R version 4.2 and 4.1. Maybe R version ≥ 4.0 is feasible.
Herein, the codes are giving an example with installing R 4.2.
In bash:

```{bash}
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt-get update
sudo apt install --no-install-recommends r-base
```

### Install dependencies

```
## Libraries for installing 'usethis' and 'devtools'.
sudo apt install libssl-dev libcurl4-openssl-dev
## Libraries for installing 'BiocManager' and its some packages.
sudo apt install libnetcdf-dev libopenbabel-dev libeigen3-dev
## Libraries For installing other graphic packages.
sudo apt install libfontconfig1-dev librsvg2-dev libmagick++-dev
```

### Install dependent R packages for MCnebula

&ensp;&ensp; 
Some data processing tools:

```
install.packages(c("data.table", "dplyr", "ggplot2", "stringr", "tidyr", "reshape2", "pbapply", "ggsci"))
```

&ensp;&ensp; 
Installing R packages from Bioconductor:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("MSnbase", "ChemmineOB"))
```

&ensp;&ensp; 
Of note, for current, 'ChemmineOB' is only available in Linux.
Fortunately, 'ChemmineOB' is **non-essential** for MCnebula.
However, the chemical structure mapping in child-nebula is depending on 'ChemmineOB'.  

&ensp;&ensp; 
Installing graphic tools:

```
install.packages("igraph", "ggraph", "svglite", "tidygraph", "rsvg", "grImport2", "ggimage", "ggtext", "ggsci")
```

&ensp;&ensp; 
For installing packages in github:

```
install.packages(c("usethis", "devtools"))
```

### Install MCnebula

MCnebula is available in github:

```
devtools::install_github("Cao-lab-zcmu/MCnebula")
```

