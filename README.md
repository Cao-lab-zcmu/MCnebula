# MCnebula

MCnebula algorithm integration in R.  
MCnebula has been published at <https://cao-lab-zcmu.github.io/MCnebula/>.
Guidance for MCnebula application:
[MCnebula_workflow](http://htmlpreview.github.io/?https://github.com/Cao-lab-zcmu/VIgnette/blob/master/MCnebula_workflow.html).

## Installation 

The following installation has been tested in Pop!_OS 20.04 (Ubuntu 20.04).

### Install R

We have tested MCnebula in R version 4.2 and 4.1. Maybe R version ≥ 4.0 is feasible.
Herein, the codes are giving an example with installing R 4.2.
In bash:

```{bash}
sudo apt-key adv \
--keyserver keyserver.ubuntu.com \
--recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository \
"deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt-get update
sudo apt install --no-install-recommends r-base
```

### Install dependencies

In bash:
```
## Libraries for installing 'usethis' and 'devtools'.
sudo apt install libssl-dev libcurl4-openssl-dev libblas-dev liblapack-dev libgfortran-11-dev
## Libraries for installing 'BiocManager' and its some packages.
sudo apt install libnetcdf-dev libopenbabel-dev libeigen3-dev
## Libraries For installing other graphic packages.
sudo apt install libfontconfig1-dev librsvg2-dev libmagick++-dev
```

### Install dependent R packages for MCnebula

Some data processing tools:

```
install.packages(c("data.table", "dplyr", "ggplot2",
                   "stringr", "tidyr", "reshape2",
                   "pbapply", "ggsci"))
```

Installing R packages from Bioconductor:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("MSnbase", "ChemmineOB"))
```

Of note, for current, 'ChemmineOB' is only available in Linux.
Fortunately, 'ChemmineOB' is **non-essential** for MCnebula.
However, the chemical structure mapping in child-nebula is depending on 'ChemmineOB'.  

Installing graphic tools:

```
install.packages(c("igraph", "ggraph", "svglite",
                 "tidygraph", "rsvg", "grImport2",
                 "ggimage", "ggtext", "ggsci"))
```

For installing packages in github:

```
install.packages(c("usethis", "devtools"))
```

### Install MCnebula

MCnebula is available in github:

```
devtools::install_github("Cao-lab-zcmu/MCnebula")
```

## Usage

### Main functions

Loading library.

```
library(ggplot2)
library(ggraph)
library(grid)
library(MCnebula)
```

MCnebula perform data collating and integration in SIRIUS project.  

First, user should initialize MCnebula at the directory of SIRIUS project.

```
path <- "my.sirius.project"
initialize_mcnebula(path)
```

The above will set some global var:

```
> .MCn.
.MCn.output         .MCn.palette_label  .MCn.palette_stat   .MCn.sirius
.MCn.palette        .MCn.palette_ppcp   .MCn.results
```

Then, for data collating:

```
collate_structure(
  exclude_element = c("Cl", "S", "P"),
  ppm_error = 20
)
build_classes_tree_list()
collate_ppcp(
  nebula_class = T,
  ## ------------------------------------- 
  ## detail control for nebula-class
  ## level 5, subclass, class, superclass
  hierarchy_priority = c(6, 5, 4, 3),
  max_classes = NA,
  ## ------------------------------------- 
  ## summarise nebula-index
  nebula_index = T,
  ## ------------------------------------- 
  ## detail control for nebula-index
  min_possess = 30,
  max_possess_pct = 0.07,
  filter_identical = c("top_hierarchy" = 4),
  identical_factor = 0.7,
  ## following maybe remove lots of classes with many low score structures,
  ## if not, set NA
  ## ------------------ 
  filter_via_struc_score = "tanimotoSimilarity", 
  struc_score_cutoff = 0.3,
  min_reached_pct = 0.6,
  ## ------------------ 
  rm_position_describe_class = T
  ## ------------------------------------- 
)
```

`collate_ppcp` is flexible with multiple augments. In addition, ignore the complexity, user can jsut run likewise:

```
collate_ppcp(
  min_possess = 30,
  max_possess_pct = 0.2
)
```

As well, the above commands will done with some global vars (data.frame or list) of which the name
begin with `.MCn.`.  

Third, for generate chemical nebulae (network):

```
generate_parent_nebula(
  rm_parent_isolate_nodes = T
)
generate_child_nebulae(
  max_edges = 5,
  ## this will write a output of .graphml, which support by Cytoscape.
  output_format = "graphml"
)
```

Fourth, visualization for nebulae:

```
visualize_parent_nebula(
  layout = "mds",
  ## map nodes color with superclass
  nodes_color = c("hierarchy" = 3)
)
visualize_child_nebulae(
  layout = "fr",
  ## herein, you can mark some features with specific color.
  ## if not, just ignore the augment.
  nodes_mark = data.frame(
    .id = c("1", "300", "500"),
    mark = c("biomarker1", "biomarker2", "biomarker3")
  ),
  palette = c(
    biomarker1 = 'pink',
    biomarker2 = 'lightblue',
    biomarker3 = 'white',
    Others = '#B8B8B8'
  ),
  nodes_size_range = c(2, 4),
  ## adjust the legend
  remove_legend_lab = T,
  legend_fill = T,
  legend_position = "bottom",
  nodes_stroke = 0,
  edges_color = "#D9D9D9"
)
```

Last, users are encouraged for in-depth visualization of child-nebula:

```
annotate_child_nebulae(
  ## string, i.e. class name in nebula-index
  nebula_name,
  layout = "fr",
  ## a table to mark color of nodes
  nodes_mark = data.frame(
    .id = c("1", "300", "500"),
    mark = c("biomarker1", "biomarker2", "biomarker3")
  ),
  plot_nodes_id = T,
  plot_structure = T,
  plot_ppcp = T,
  ## manually define the color of nodes
  palette = c(
    biomarker1 = 'pink',
    biomarker2 = 'lightblue',
    biomarker3 = 'white',
    Others = '#B8B8B8'
  ),
  ## feature quantification table
  ratio_df = data.frame(
    .id = 1:1000,
    sample1 = rnorm(1000, 1E5, 1E4),
    sample2 = rnorm(1000, 1E6, 5E4),
    sample3 = rnorm(1000, 1E4, 5E3),
  ),
  ## A vector of the hex color with names or not
  palette_stat = c(
    sample1 = 'blue',
    sample2 = 'yellow',
    sample3 = 'red'
  ),
  ## control nodes size in child-nebula, zoom in or zoom out globally.
  global.node.size = 0.8,
  ## the args of `ggplot::theme`
  theme_args = list(
    panel.background = element_rect(),
    panel.grid = element_line()
  ),
  return_plot = F
)
```

### Other functions

...

### Standalone tools

**Use InChiKey planar to get all possible InChIKey.**
`pubchem_curl_inchikey` function use `RCurl` package and through PubChem API to get InChIKey.
We design this function with the stability of the network in consideration.
Therefore, the process of getting data will not fail due to `Status: 503`.

```
pubchem_curl_inchikey(
  ## a vector of InChIKey planar
  inchikey2d,
  ## output path
  dir,
  ## number of threads
  curl_cl = NULL,
  ## each InChiKey planar will get a .TSV file. Here, this setting could gather them
  ## as 'list' into .rdata. Indeed, this can save storage space.
  gather_as_rdata = T
)
## if 'gather_as_rdata' == T, user can use following function to get results.
inchikey_set <- extract_rdata_list(
  paste0(dir, "/", "inchikey.rdata")
)
```

**Use PubChem CID to get synonyms of compounds**.
Note that `pubchem_curl_inchikey` not only get InChIKey, but also PubChem CID.
Herein, `pubchem_get_synonyms` is used to get synonyms via PubChem CID.

```
pubchem_get_synonyms(
  ## a vector of PubChem CID
  cid,
  ## output path
  dir,
  ## number of threads
  curl_cl = NULL,
  gather_as_rdata = T
)
syno_set <- extract_rdata_list(
  paste0(dir, "/", "cid.rdata")
)
```

