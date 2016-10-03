# gfunc_liftover()
#
# A function to perform LiftOVer of coordinates.
# Users Bioconductor's LiftOver function (rtracklayer), which is a
# native implementation of UCSC's LiftOver command line tool.
# 
# Reference: https://www.biostars.org/p/65558/  
# Nick Burns
# Oct, 2016

# Parameters:
#   datafile: full path to input datafile
#   chainfile: full path to appropriate chain file
#   bed: vector ("CHR", "POS", "SNP") matching the corresponding column names
#        in the input datafile.
#
# Output:
#   a merged data.table containing only those positions which were mapped over successfully.
gfunc_liftover <- function (datafile, 
                            chainfile,
                            bed = c("CHR", "POS", "SNP")) {
    
    data <- data.table::fread(datafile)
    data.bed <- data[, bed, with = FALSE]
    
    # reformat the CHR column (prefix: chr)
    if (any(grepl("chr", data[[bed[1]]]))) {
        data.bed[, CHR := data[[bed[1]]]]
    } else {
        data.bed[, CHR := paste0("chr", data[[bed[1]]], sep = "")]
    }
    
    # wrangle into GRanges object, then perform liftover
    data_ <- GenomicRanges::GRanges(data.bed[, .(CHR, 
                                                 Start = data.bed[[bed[2]]],
                                                 End = data.bed[[bed[2]]] + 1,
                                                 data.bed[[bed[3]]])])
    names(GenomicRanges::mcols(data_)) <- bed[3]
    data_ <- data.table::as.data.table(
                rtracklayer::liftOver(data_,
                                      rtracklayer::import.chain(chainfile))
                )
    
    # merge in new coordinates
    data.table::setkeyv(data_, bed[3])
    data.table::setkeyv(data, bed[3])
    data <- data[data_[, .(data_[[bed[3]]], start, end)]]
  
    return (data)
}