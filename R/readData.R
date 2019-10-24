#' Read allele data from a file
#'
#' @param fileName
#'
#' @return
#' @export
#'
#' @examples
#' inFile = system.file("extdata", "infile.txt", package = "ConStruct", mustWork = TRUE)
#' readData(inFile, missing = 0)
readData = function(fileName, locusNames = FALSE, delim = "[,[:space:]]+", missing = NA){
  Lines = readLines(fileName)
  numProfiles = length(Lines)

  if(locusNames){
    Loci = strsplit(Lines[1], delim)[[1]]
    numLoci = length(Loci)
    Lines = Lines[-1]
  }else{
    Loci = NULL
  }

  Profiles = strsplit(Lines, delim)
  numAllelesPerProfile = sapply(Profiles, length)

  if(length(unique(numAllelesPerProfile)) != 1){
    stop("Each profile must have the same number of alleles")
  }

  if(is.null(Loci)){
    numAllelesPerProfile = unique(sapply(Profiles, length))
    if(numAllelesPerProfile %% 2 != 0){
      stop("There must be an even number of alleles per profile (2 per locus)")
    }

    numLoci = numAllelesPerProfile / 2L
    Loci = paste0("Locus.", 1:numLoci)
  }

  Profiles = do.call(rbind, Profiles)

  if(!is.na(missing)){
      Profiles = gsub(paste0("^", missing, "$"), NA, Profiles)
  }

  Counts = Freqs = Alleles = list(length = numLoci)
  numAllelesPerLocus = rep(0, numLoci)
  numMissing = numObs = numHoms = rep(0, numLoci)

  Recoded.Profiles = matrix(0, nrow = numProfiles, ncol = 2 * numLoci)

  for(loc in 1:numLoci){
    i1 = 2 * loc -1
    i2 = i1 + 1

    Counts[[loc]] = table(c(Profiles[,i1:i2]), useNA = "no")
    Freqs[[loc]] = Counts[[loc]] / sum(Counts[[loc]])
    Alleles[[loc]] = names(Counts[[loc]])
    numMissing[loc] = sum(is.na(Profiles[, i1:i2]))
    numObs[loc] = 2 * numProfiles - numMissing[loc]
    numAllelesPerLocus[loc] = length(Counts[[loc]])

    ## Recode Profiles

    Recoded.Profiles[, i1] = match(Profiles[, i1], Alleles[[loc]])
    Recoded.Profiles[, i2] = match(Profiles[, i2], Alleles[[loc]])
  }

  return(list(Profiles = Recoded.Profiles,
              Profiles.raw = Profiles,
              Counts = Counts,
              Freqs = Freqs,
              Loci = Loci,
              Alleles = Alleles,
              numProfiles = numProfiles,
              numLoci = numLoci,
              numAllelesPerLocus = numAllelesPerLocus,
              numObs = numObs,
              numMissing = numMissing)
  )
}
