defaultW <- getOption("warn")
writeLines(paste('\nLoading libraries This may take a while.\n'))
options(warn = -1)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(colorout))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(karyoploteR))
options(warn = defaultW)

option_list = list(
  make_option(c("-i", "--ibam"), help="Input bam file for read density. Default is 'output.primary.mapped.sort.bam'", type="character", metavar="character"),
    make_option(c("-g","--genomefile"), help="Path for genome file (karyoploteR). This must be provided.", type="character", metavar="character"),
    make_option(c("-o","--output"), help="Output file name.", type="character", metavar="character"),
    make_option(c("-t","--title"), help="Title convention for plots, e.g. 'Title (1MB window-size). Default is 'Read density for reads mapping to genome (1MB window-size)', etc.", type="character", metavar="character"),
    make_option(c("-r","--res"), help="Resolution of output images. Default is 300. Maximum value is 900.", type="integer", metavar="integer")
    # make_option(c("-s","--svg"), help="Output in .svg. Default is jpg.", type="character", metavar="character")

);

###################### Input files and parameters ######################
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genomefile)){
  warning('Path to genome file must be provided.')
  q()
} else {
  genome_file <- opt$genomefile
  custom.genome <- toGRanges(genome_file)
}

# if (is.null(opt$svg)){
# image.end <- '.jpg'
# } else {
# image.end <- '.svg'
# }

if (is.null(opt$res)) {
  w <- 2400
  h <- 1800
  r <- 300
} else if (opt$res <= 900 ){
  r <- opt$res
  w <- 8*r
  h <- 6*r
} else {
  warning('Max resolution is 900')
  q()
}


if (is.null(opt$output)) {
  output_1mb <- 'Figure_read_density_1MB.jpg'
  output_100kb <- 'Figure_read_density_100KB.jpg'
  output_10kb <- 'Figure_read_density_10KB.jpg'
  output_01kb <- 'Figure_read_density_1KB.jpg'
} else {
  output_1mb <- paste0(opt$output,'_read_density_1MB.jpg')
  output_100kb <- paste0(opt$output,'_read_density_100KB.jpg')
  output_10kb <- paste0(opt$output,'_read_density_10KB.jpg')
  output_01kb <- paste0(opt$output,'_read_density_1KB.jpg')
}

if (is.null(opt$title)) {
  title_1mb <- 'Read density for reads mapping to genome (1 mb window-size)'
  title_100kb <- 'Read density for reads mapping to genome (100 kb window-size)'
  title_10kb <- 'Read density for reads mapping to genome (10 kb window-size)'
  title_01kb <- 'Read density for reads mapping to genome (1 kb window-size)'
} else {
  title_1mb <- paste(opt$title,'(1MB window-size)')
  title_100kb <- paste(opt$title,'(100 kb window-size)')
  title_10kb <- paste(opt$title,'(10 kb window-size)')
  title_01kb <- paste(opt$title,'(1 kb window-size)')
}

if (is.null(opt$ibam)) {
  warning('Input bam file must be provided')
} else {
  density_bam <- opt$ibam
}

writeLines(paste('Input file: ',density_bam))
writeLines(paste('Path for genome file', genome_file))
writeLines(paste('Image resolution:', r))


###################### Colors ######################
samplix_green <- "#96bf21"
samplix_pastel <- "#8bcab5"
samplix_blue <- "#006c9a"
samplix_red <- "#ce1d31"
samplix.grey <- "#c2c2c1"
samplix.60.green <- "#C0D97A"
samplix.60.grey <- "#DADADA"
samplix.60.red <- "#E27783"
samplix.60.pastel <- "#B9DFD3"
samplix.60.blue <- "#66A7C2"



### Bai index size for maximum chromosome length that can be handled by standard karyoploteR
BAI_INDEX_SIZE <- 536870912

###################### Functions ######################
### Custom BAM density plot for long chromosomes
kpPlotBAMDensity_custom <- function(karyoplot, data=NULL, window.size=1e6, normalize=FALSE, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, col="gray80", border=NA, clipping=TRUE,...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotBAMDensity: 'karyoplot' must be a valid 'KaryoPlot' object"))
  if(is.character("data")) {
  data <- Rsamtools::BamFile(file = data)
  }
  if(!methods::is(data, "BamFile")) stop(paste0("In kpPlotBAMDensity: 'data' must be a character or a 'BamFile' object."))

  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())

  #use tileGenome to create windows only on the visible part of the genome, that
  #is in the karyoplot$plot.region
  plot.region.lengths <- setNames(width(karyoplot$plot.region), as.character(seqnames(karyoplot$plot.region)))
  windows <- tileGenome(plot.region.lengths, tilewidth = window.size, cut.last.tile.in.chrom = TRUE)
  seqinfo(windows) <- GenomeInfoDb::Seqinfo(seqnames=seqlevels(windows)) #remove the seqlength info from seqinfo to avoid a potential out-of-bounds warning when shifting the windows

  #Now, move the windows start(karyoplot$plor.region) bases to the right.
  #It's only necessary when zoomed or with chromosomes not starting at position 1
  windows <- shift(windows, shift=start(karyoplot$plot.region[seqnames(windows)])-1)
  windows_data <- data.frame(windows)
  windows_bed <- windows_data %>% select(seqnames,start,end)
  write.table(windows_bed, file='windows.bed', sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE)

  handle <- '/usr/local/bedtools/bin/bedtools'
  args <- c('intersect', '-a','windows.bed', '-b', 'output.sort.bam', '-C', '-bed')
  records <- read.table(text = system2(handle, args, stdout=TRUE), sep = "\t", header=FALSE)
  dens <- as.numeric(unlist(records$V4))
  file.remove('windows.bed')

  if(normalize==TRUE) {
  total.reads <- sum(Rsamtools::idxstatsBam(file=data)$mapped)
  dens <- dens/total.reads
  }

  if(is.null(ymax)) {
  ymax <- max(dens)
  }

  #Specify the missing colors if possible
  if(is.null(border) & !is.null(col) & !is.na(col)) {
  border=darker(col, amount = 100)
  }
  if(is.null(col) & !is.null(border) & !is.na(col)) {
  col=lighter(border)
  }
  if(is.na(col) & is.null(border)) {
  border <- "black"
  }

  karyoplot <- kpBars(karyoplot, data=windows, y0=0, y1=dens,ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, border=border, col=col, clipping=clipping, ...)

  karyoplot$latest.plot <- list(funct="kpPlotBAMDensity", computed.values=list(density=dens, windows=windows, max.density=max(dens)))

  invisible(karyoplot)
}

density_plot <- function(poutput,wSize,pTitle){
  writeLines(paste('Processing read density graph',wSize,'bp. This may take a while'))
  jpeg(poutput, width=w,height=h, res=r)
  pp <- getDefaultPlotParams(plot.type=1)
  pp$ideogramheight <- 0.1
  kp <- plotKaryotype(genome = custom.genome, cex=0.4, plot.params = pp)
  if (max.chromosome > BAI_INDEX_SIZE){
  bam.read <- kpPlotBAMDensity_custom(kp, data=density_bam, col=samplix_blue, r0=0, r1=0.5, window.size=wSize)
  } else {
  bam.read <- kpPlotBAMDensity(kp, data=density_bam, col=samplix_blue, r0=0, r1=0.5, window.size=wSize)
  }
  ymax.read <- bam.read$latest.plot$computed.values$max.density
  kpAxis(kp, ymax=ymax.read, numticks=2, cex=0.25, r0=0, r1=0.5, col=samplix.grey, side=2)
  kpAddBaseNumbers(kp, tick.dist=tickDist, add.units=TRUE, minor.ticks=TRUE, minor.tick.dist = minorDist, cex=0.25, tick.len=10, minor.tick.len=4)
  kpAddMainTitle(kp, main=pTitle, cex=0.5)
  dev.off()
  writeLines(paste('Output:', poutput))
}


###################### Tick distances for plot ######################

k.genome <- read.table(genome_file, header=FALSE)
max.chromosome <- max(k.genome$V3)
if (max.chromosome > 250000000){
  tickDist <- 20000000
  minorDist <- 5000000
  } else {
  tickDist <- 10000000
  minorDist <- 2000000
}


###################### Output files ######################
if (max.chromosome < 1e5){
  density_plot(output_10kb,1e4,title_10kb)
  density_plot(output_01kb,1e3,title_1kb)
  } else if (max.chromosome < 1e6){
  density_plot(output_100kb,1e5,title_100kb)
  density_plot(output_10kb,1e4,title_10kb)
  } else {
  density_plot(output_1mb,1e6,title_1mb)
  density_plot(output_100kb,1e5,title_100kb)
}


q()
