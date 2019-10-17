#' maxLikelihood
#'
#' @param data the input file
#' @param max.alleles max.alleles places an uppermost limit on the number of alleles considered
#' @param resolution resolution is the resolution of the F parameter
#'
#' @return
#' @export
#'
#' @examples
#' maxLikelihood(data = "infile.txt", max.alleles = 1000, resolution = 100)
maxLikelihood = function(data, max.alleles, resolution) {
  # Input file defined and NAs converted to 0
  infile = read.table(data)
  infile [is.na(infile)] <- 0
  # Population size is defined by the number of lines
  N = nrow(infile)
  # Number of loci is defined by number of columns
  num.loc = ncol(infile) / 2
  # Initialise variables
  all.count = matrix(0, num.loc, max.alleles)
  all.freq = matrix(0, num.loc, max.alleles)
  homozygous.frequency = matrix(0, num.loc, max.alleles)
  count.alleles = matrix(0, max.alleles)
  likelihood.individual = matrix(0, N)
  ln.likelihood = matrix(0, resolution)
  e.likelihood = matrix(0, resolution)
  hom.total = matrix(0, num.loc)
  #
  #set up input data
  #
  # Read diploid genotypes as sets of two alleles
  allele.1 = matrix(0, N, (2 * num.loc))
  allele.2 = matrix(0, N, (2 * num.loc))
  for (i in 1:N) {
    for (j in seq(1, by = 2, (2 * num.loc - 1))) {
      allele.1[i, j] = infile[i, j]
    }
  }
  for (i in 1:N) {
    for (j in seq(2, by = 2, (2 * num.loc))) {
      allele.2[i, j] = infile[i, j]
    }
  }
  allele.1.seq = seq(1, 2 * num.loc, by = 2)
  allele.2.seq = seq(2, 2 * num.loc, by = 2)
  for (i in 1:N) {
    for (j in 1:num.loc) {
      allele.1[i, j] = allele.1[i, allele.1.seq[j]]
    }
  }
  for (i in 1:N) {
    for (j in 1:num.loc) {
      allele.2[i, j] = allele.2[i, allele.2.seq[j]]
    }
  }

  # Calculate allele frequencies
  for (i in 1:N) {
    for (j in 1:num.loc) {
      if (allele.1[i, j] > 0) {
        all.count[j, allele.1[i, j]] = all.count[j, allele.1[i, j]] + 1
      }
    }
  }
  for (i in 1:N) {
    for (j in 1:num.loc) {
      if (allele.2[i, j] > 0) {
        all.count[j, allele.2[i, j]] = all.count[j, allele.2[i, j]] + 1
      }
    }
  }
  #
  locus.n = matrix(0, num.loc)
  for (i in 1:num.loc) {
    for (j in 1:max.alleles) {
      locus.n[i] = sum(all.count[i,])
    }
  }
  for (i in 1:num.loc) {
    for (j in 1:max.alleles) {
      all.freq[i, j] = all.count[i, j] / locus.n[i]
    }
  }
  #
  # Calculate maximum Fst using Hedrick’s equation (1-Hs)/Hs, where
  # Hs is the expected heterozygosity
  homozygous.frequency = matrix(0, num.loc, max.alleles)
  for (i in 1:num.loc) {
    for (j in 1:max.alleles) {
      if (all.freq[i, j] > 0) {
        homozygous.frequency[i, j] = homozygous.frequency[i, j] + (all.freq[i, j]) ^ 2
      }
    }
  }
  fst.max = (sum(homozygous.frequency / num.loc) / (1 - sum(homozygous.frequency) / num.loc))
  if (fst.max >= 1) {
    fst.max = 1
  }
  #
  # Calculate probability of multilocus genotypes for range of f values
  #
  # substructure.homozygous calculates the probability of homozygous genotypes
  # given a value of f
  substructure.homozygous = function(f, all.freq.1) {
    all.freq.1 * (f + (1 - f) * all.freq.1)
  }
  # substructure.heterozygous calculates the probability of homozygous genotypes
  # given a value of f
  substructure.heterozygous = function(f, all.freq.1, all.freq.2) {
    2 * all.freq.1 * all.freq.2 * (1 - f)
  }
  # Initialise f value
  f = 0
  # res keeps track of number of number of iterations
  res = 1
  # Initialise ln.likelihood values (the ln of the probabilities of the
  # individual geneotypes
  ln.likelihood = matrix(0, resolution)
  # likelihood.locus = the probabilities of the individual genotypes, for
  # each individual and each locus
  likelihood.locus = matrix(0, N, num.loc)
  # The probability of the multicolous dataset is calculated for a range of
  # values of f (0 - fst.max)
  repeat {
    likelihood.locus = matrix(0, N, num.loc)
    # Take missing genotypes and make likelihood = 1
    for (i in 1:N) {
      for (j in 1:num.loc) {
        if (allele.1[i, j] == 0) {
          likelihood.locus[i, j] = 1
        }
      }
    }
    # Identify homozygous genotypes and calculate probability, given f
    for (i in 1:N) {
      for (j in 1:num.loc) {
        if (allele.1[i, j] > 0) {
          if (allele.1[i, j] == allele.2[i, j]) {
            likelihood.locus[i, j] = substructure.homozygous(f, all.freq[j, allele.1[i, j]])
          }
        }
      }
    }
    # Identify heterozygous genotypes and calculate probability, given f
    for (i in 1:N) {
      for (j in 1:num.loc) {
        if (allele.1[i, j] > 0) {
          if (allele.1[i, j] != allele.2[i, j]) {
            likelihood.locus[i, j] = substructure.heterozygous(f, all.freq[j, allele.1[i, j]], all.freq[j, allele.2[i, j]])
          }
        }
      }
    }
    # Take log of the product of each individual’s multilocus genotype
    for (i in 1:N) {
      likelihood.individual[i] = log(prod(likelihood.locus[i,]))
    }
    # ln.likelihood[res] is the log.probability of the entire dataset, given a value of f
    ln.likelihood[res] = sum(likelihood.individual[, 1])
    # To identify the maximum likelihood, the maximum value is
    # subtracted from each of the ln.likelihood values
    max.value = max(ln.likelihood)
    temp = matrix(0, resolution)
    for (i in 1:resolution) {
      temp[i] = ln.likelihood[i] - max.value
    }
    for (i in 1:resolution) {
      ln.likelihood[i] = temp[i]
    }
    res = res + 1
    f.res = fst.max / resolution
    f = f + f.res
    if (res == (resolution + 1)) {
      break
    }
  }
  #
  # The maximum likelihood is then the ln.likelihood value that = 0
  # The corresponding value of F is found:
  for (i in 1:resolution) {
    if (ln.likelihood[i] == 0) {
      ML = ((i - 1) / resolution) * fst.max
    }
  }
  # A plot of the distribution of F; the f-axis goes from 0 - fst.max
  f = seq(0, by = f.res, (fst.max - f.res))
  f.axis <<- f
  #
  #
  # Likelihood values are converted back from logs
  for (i in 1:resolution) {
    e.likelihood[i] = exp(ln.likelihood[i])
  }
  # Area of distribution sums to 1
  total = sum(e.likelihood)
  e.likelihood = e.likelihood / total
  probability <<- e.likelihood
  plot(f,
       e.likelihood,
       type = "l",
       xlab = "F",
       ylab = "Likelihood")
  #
  #
  #
  writeLines("Maximum value of Fst = ")
  print(fst.max)
  writeLines("Maximum Likelihood value of Fst = ")
  print(ML)
  writeLines("G = ")
  G = exp(log(max(e.likelihood)) - 2)
  print(G)
}
# end of function
