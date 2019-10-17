#' Read allele data from a file
#'
#' @param fileName
#'
#' @return
#' @export
#'
#' @examples
readData = function(fileName, locusNames = FALSE, delim = "[,[:space:]]+"){
  Lines = readLines(fileName)
  numProfiles = length(Lines)

  if(locusNames){
    loci = strsplit(Lines[1], delim)[[1]]
    numLoci = length(loci)
    Lines = Lines[-1]
  }else{
    loci = NULL
  }

  Profiles = strsplit(Lines, delim)
  numAllelesPerProfile = sapply(Profiles, length)

  if(length(unique(numAllelesPerProfile)) != 1){
    stop("Each profile must have the same number of alleles")
  }

  if(is.null(loci)){
    numAllelesPerProfile = unique(sapply(Profiles, length))
    if(numAllelesPerProfile %% 2 != 0){
      stop("There must be an even number of alleles per profile (2 per locus)")
    }

    numLoci = numAllelesPerProfile / 2L
    loci = paste0("Locus.", 1:numLoci)
  }

  Profiles = do.call(rbind, Profiles)

  Counts = Freqs = list(length = numLoci)
  numAllelesPerLocus = rep(0, numAlleles)


  return(Profiles)

}
