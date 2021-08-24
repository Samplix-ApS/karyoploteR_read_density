# karyoploteR_read_density
## Table of Contents
- [Description](#descrip)
- [Getting Started](#get_started)
  - [Dependencies](#dependencies)
  - [Installation](#install_)
- [Executing program](#execute_program)
- [Help](#help_)
- [Authors](#authors_)

## <a name="descrip"></a> Description
Creating read density plots with karyoploteR to give an overview of where reads have aligned to the genome. 

<p align="center">
  <image src="https://user-images.githubusercontent.com/60882704/130640561-8f0b6f73-d4c9-4a2a-8626-4655cf6fe15b.jpg">
    </p>
  
## <a name="get_started"></a> Getting Started
Custom genome file must be generated using prep reference before creating karyoplots.
  
### <a name="dependencies"></a>Dependencies
* R 4.0.3 (later versions likely also work)
  
* R library optparse v 1.6.6

* R library karyoploteR 1.16.0

### <a name="install_"></a> Installation

Download the R script

## <a name="execute_program"></a> Executing program

To creat karyoplots with read densities:
```
Rscript karyoploteR_read_density.R -i BAM -g GENOMEFILE <optional>
```
  
The following parameters can be used:
  
```
-i	input text file containg read name and read lengths. Required
-g	Path for genome file (karyoploteR). Required.
-t	Title convention for plots, e.g. 'Title (1MB window-size)'. Default is 'Read density for reads mapping to genome (1XX window-size)'. 
-o	output file name. Output convention is 'OUTPUT_read_density_1XX.jpg'. Default is 'Figure_read_density_1XX.jpg'.
-r	Resolution of output images. Default is 300. Maximum value is 900.
```


## <a name="help_"></a>Help

```
Rscript karyoploteR_read_density.R -h
```
## <a name="refs_"></a>References
R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
  
Trevor L Davis (2020). optparse: Command Line Option Parser. R package version 1.6.6. https://CRAN.R-project.org/package=optparse
  
Jakson Aquino and Dominique-Laurent Couturier (2019). colorout:   Colorize R Output on Terminal Emulators. R package version 1.2-2.   https://github.com/jalvesaq/colorout
  

Gel B, Serra E (2017). "karyoploteR : an R / Bioconductor package to plot customizable genomes displaying arbitrary data." _Bioinformatics_,*33*(19), 3088-3090. doi: 10.1093/bioinformatics/btx346 https://doi.org/10.1093/bioinformatics/btx346.

  
## <a name="authors_"></a>Authors

Camille Johnston 


