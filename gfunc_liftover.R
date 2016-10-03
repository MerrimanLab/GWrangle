# gfunc_liftover()
#
# A function to perform LiftOVer of coordinates
# from one version of the human genome to another.
# Users Bioconductor's LiftOver function (rtracklayer), which is a
# native implementation of UCSC's LiftOver command line tool.
# 
# Reference: https://www.biostars.org/p/65558/  
# Nick Burns
# Oct, 2016

gfunc_liftover <- function (datafile, 
                            chainfile,
                            bed_columns = c("CHR", "POS", "SNP")) {
    
    data <- data.table::fread(datafile)
    data.bed <- data[, bed_columns, with = FALSE]
    
    # reformat the CHR column (prefix: chr)
    if (any(grepl("chr", data[[bed_columns[1]]]))) {
        data.bed[, CHR := data[[bed_columns[1]]]]
    } else {
        data.bed[, CHR := paste0("chr", data[[bed_columns[1]]], sep = "")]
    }
    
    data_ <- GenomicRanges::GRanges(data.bed[, .(CHR = CHR, 
                                                 Start = bed_columns[2],
                                                 End = bed_columns[2] + 1,
                                                 SNP = bed_columns[3])])
    data_ <- unlist(rtracklayer::liftOver(data_, 
                                          rtracklayer::import.chain(chainfile)))
    
    data <- merge(data, data_, by = bed_columns[3])
    
    return (data)
}