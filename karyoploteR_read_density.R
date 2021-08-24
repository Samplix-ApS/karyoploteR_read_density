defaultW <- getOption("warn")
options(warn = -1)
library(optparse)
library(colorout)
options(warn = defaultW)

option_list = list(
  make_option(c("-i", "--ibam"), help="Input bam file for read density. Default is 'output.primary.mapped.sort.bam'", type="character", metavar="character"),
    make_option(c("-g","--genomefile"), help="Path for genome file (karyoploteR). This must be provided.", type="character", metavar="character"),
    make_option(c("-o","--output"), help="Output file name.", type="character", metavar="character"),
    make_option(c("-t","--title"), help="Title convention for plots, e.g. 'Title (1MB window-size). Default is 'Read density for reads mapping to genome (1MB window-size)', etc.", type="character", metavar="character"),
    make_option(c("-r","--res"), help="Resolution of output images. Default is 300. Maximum value is 900.", type="integer", metavar="integer")
    # make_option(c("-s","--svg"), help="Output in .svg. Default is jpg.", type="character", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genomefile)){
  warning('Path to genome file must be provided.')
  q()
} else {
  genome_file <- opt$genomefile
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

options(warn = -1)
writeLines(paste('\nLoading karyoploteR. This may take a while.\n'))
library(karyoploteR)
options(warn = defaultW)
writeLines(paste('\nkaryoploteR loaded.\n'))



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


###input KaryoploteR figures
custom.genome <- toGRanges(genome_file)


density_plot <- function(poutput,wSize,pTitle){
  jpeg(poutput, width=w,height=h, res=r)
  pp <- getDefaultPlotParams(plot.type=1)
  pp$ideogramheight <- 0.1
  kp <- plotKaryotype(genome = custom.genome, cex=0.4, plot.params = pp)
  bam.read <- kpPlotBAMDensity(kp, data=density_bam, col=samplix_blue, r0=0, r1=0.5, window.size=wSize)
  ymax.read <- bam.read$latest.plot$computed.values$max.density
  kpAxis(kp, ymax=ymax.read, numticks=2, cex=0.25, r0=0, r1=0.5, col=samplix.grey, side=2)
  kpAddBaseNumbers(kp, tick.dist=tickDist, add.units=TRUE, minor.ticks=TRUE, minor.tick.dist = minorDist, cex=0.25, tick.len=10, minor.tick.len=4)
  kpAddMainTitle(kp, main=pTitle, cex=0.5)
  dev.off()
}

k.genome <- read.table(genome_file, header=FALSE)
max.chromosome <- max(k.genome$V3)
if (max.chromosome > 250000000){
  tickDist <- 20000000
  minorDist <- 5000000
} else {
  tickDist <- 10000000
  minorDist <- 2000000
}


### 1 region
if (max.chromosome < 1e5){
  writeLines(paste('Processing read density graph 10 kb window. This may take a while'))
  density_plot(output_10kb,1e4,title_10kb)
  writeLines(paste('Processing read density graph 1 kb window. This may take a while'))
  density_plot(output_01kb,1e3,title_1kb)
} else if (max.chromosome < 1e6){
  writeLines(paste('Processing read density graph 100 kb window. This may take a while'))
  density_plot(output_100kb,1e5,title_100kb)
  writeLines(paste('Processing read density graph 10 kb window. This may take a while'))
  density_plot(output_10kb,1e4,title_10kb)
} else {
  writeLines(paste('Processing read density graph 1 mb window. This may take a while'))
  density_plot(output_1mb,1e6,title_1mb)
  writeLines(paste('Output', output_1mb))
  writeLines(paste('Processing read density graph 100 kb window. This may take a while'))
  density_plot(output_100kb,1e5,title_100kb)
  writeLines(paste('Output', output_100kb))
}


q()
