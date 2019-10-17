#' Simulate a population with desired levels of substructure and inbreeding
#'
#' @param N the total sample size
#' @param num.loc the number of loci
#' @param fst the value of Fst that is to be simulated between two populations
#' @param r.actual the inbreeding coefficient of the inbred individuals
#' @param c  the proportion of the population inbred to degree r.actual
#' @param r.consider the value of the inbreeding coefficient being considered for the
#     analysis of the simulated dataset
#' @param max.alleles places an uppermost limit on the number of alleles considered
#' @param f.resolution the resolution on the Fst parameter
#' @param c.resolution the resolution on the c parameter
#' @param iteration the number of iterations of the simulation run through in order
#     to arrive at the specified, simulated Fst
#'
#' @return
#' @export
#'
#' @examples
#' simulate(200, 12, 0.05, 0.05, 0.5, 0.05, 100, 100, 100, 10000)
simulate = function(N,
                    num.loc,
                    fst,
                    r.actual,
                    c,
                    r.consider,
                    max.alleles,
                    f.resolution,
                    c.resolution,
                    iteration) {
  # It is necessary to specify number of alleles per locus
  # ***********************************
  # specify number of alleles per locus
  # ***********************************
  num.alleles = c(4, 4, 4, 4, 4)
  #
  # Check to see if specified number of alleles corresponds with the specified number of
  # loci (num.loc):
  if (length(num.alleles) != num.loc) {
    stop("number of alleles per locus does not match number of loci specified.")
  }
  
  # Generate random allele frequencies for each locus for two pops
  # x.old.1, for example are the initial allele frequencies in population 1
  # x.new.1, for example are the updated allele frequencies in population 1
  # Once the required Fst value has been found, these frequencies are referred to
  # as allele.freq.1 and allele.freq.2
  # Initialise variables
  
  limit = array(0, dim = c(num.loc))
  x.old.1 = matrix(0, num.loc, max.alleles)
  
  x.old.2 = matrix(0, num.loc, max.alleles)
  
  x.new.1 = matrix(0, num.loc, max.alleles)
  
  x.new.2 = matrix(0, num.loc, max.alleles)
  allele.freq.1 = matrix(0, num.loc, max.alleles)
  allele.freq.2 = matrix(0, num.loc, max.alleles)
  homo.1.F = matrix(0, num.loc, max.alleles)
  homo.1.O = matrix(0, num.loc, max.alleles)
  homo.cumulative.1.F = matrix(0, num.loc, max.alleles)
  homo.cumulative.1.O = matrix(0, num.loc, max.alleles)
  het.1.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
  het.1.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
  het.cumulative.1.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
  het.cumulative.1.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
  homo.2.F = matrix(0, num.loc, max.alleles)
  homo.2.O = matrix(0, num.loc, max.alleles)
  homo.cumulative.2.F = matrix(0, num.loc, max.alleles)
  homo.cumulative.2.O = matrix(0, num.loc, max.alleles)
  het.2.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
  het.2.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
  het.cumulative.2.F = array(0, dim = c(num.loc, max.alleles, max.alleles))
  het.cumulative.2.O = array(0, dim = c(num.loc, max.alleles, max.alleles))
  genotype.1.allele.1.F = matrix(0, N, num.loc)
  genotype.1.allele.2.F = matrix(0, N, num.loc)
  genotype.2.allele.1.F = matrix(0, N, num.loc)
  genotype.2.allele.2.F = matrix(0, N, num.loc)
  genotype.1.allele.1.O = matrix(0, N, num.loc)
  genotype.1.allele.2.O = matrix(0, N, num.loc)
  genotype.2.allele.1.O = matrix(0, N, num.loc)
  genotype.2.allele.2.O = matrix(0, N, num.loc)
  all.count.1 = matrix(0, num.loc, max.alleles)
  all.count.2 = matrix(0, num.loc, max.alleles)
  all.freq = matrix(0, num.loc, max.alleles)
  x.average = array(0, dim = c(max.alleles))
  
  x.average.new = array(0, dim = c(max.alleles))
  numerator = array(0, dim = c(max.alleles))
  
  denominator = array(0, dim = c(max.alleles))
  f.limit = f.resolution + 1
  c.limit = c.resolution + 1
  sum.allele = array(0, dim = c(num.loc))
  fst.old = array(0, dim = c(num.loc))
  fst.new = array(0, dim = c(num.loc))
  likelihood.individual.c = array(0, dim = c(N))
  likelihood.individual.s = array(0, dim = c(N))
  ln.likelihood.individual = array(0, dim = c(N))
  ln.likelihood.cs = matrix(0, f.resolution, c.resolution)
  e.likelihood.cs = matrix(0, f.resolution, c.resolution)
  
  num = 0
  #
  #
  # Initial allele frequencies are uniformly distributed
  
  for (i in 1:num.loc) {
    for (j in 1:num.alleles[i]) {
      x.old.1[i, j] = 1 / num.alleles[i]
      
      x.old.2[i, j] = 1 / num.alleles[i]
      
    }
    
  }
  #
  #
  # Specifies maximum number of alleles given by num.alleles
  for (i in 1:num.loc) {
    limit[i] = num.alleles[i]
    
    #
    
    # Calculates initial Fst between two simulated populations
    for (j in 1:num.alleles[i]) {
      x.average[j] = (x.old.1[i, j] + x.old.2[i, j]) / 2
      numerator[j] = (x.old.1[i, j] - x.average[j]) ^ 2
      denominator[j] = x.average[j] * (1 - x.average[j])
    }
    # Randomly explores similar allele frequencies in search of required Fst
    fst.old[i] = sum(numerator) / sum(denominator)
    for (k in 1:iteration) {
      #  Number of iterations to find desired Fst
      for (j in 1:(limit[i] - 1)) {
        #  For each allele
        sign = runif(1, 0, 1)           #  Random number between 0 - 1
        if (sign > 0.5)
          x = 1
        else
          x = -1
        if (x.old.1[i, j] + (x * 0.01) > 0 &
            x.old.1[i, j] + (x * 0.01) < 1)
          x.new.1[i, j] = x.old.1[i, j] + (x * 0.01)
        sign = runif(1, 0, 1)
        if (sign > 0.5)
          x = 1
        else
          x = -1
        if (x.old.2[i, j] + (x * 0.01) > 0 &
            x.old.2[i, j] + (x * 0.01) < 1)
          x.new.2[i, j] = x.old.2[i, j] + (x * 0.01)
      }
      if (sum(x.new.1[i, 1:(limit[i] - 1)]) <= 1) {
        x.new.1[i, limit[i]] = 1 - sum(x.new.1[i, 1:limit[i] - 1])
      }
      if (sum(x.new.2[i, 1:(limit[i] - 1)]) <= 1) {
        x.new.2[i, limit[i]] = 1 - sum(x.new.2[i, 1:limit[i] - 1])
      }
      if (sum(x.new.1[i, ]) == 1 & sum(x.new.2[i, ]) == 1) {
        for (j in 1:limit[i]) {
          x.average.new[j] = (x.new.1[i, j] + x.new.2[i, j]) / 2
          numerator[j] = (x.new.1[i, j] - x.average.new[j]) ^ 2
          denominator[j] = x.average.new[j] * (1 - x.average.new[j])
        }
        fst.new[i] = sum(numerator) / sum(denominator)    #fst across all alleles
        if ((fst - fst.new[i]) ^ 2 < (fst - fst.old[i]) ^ 2) {
          for (j in 1:limit[i]) {
            x.old.1[i, j] = x.new.1[i, j]
            x.old.2[i, j] = x.new.2[i, j]
            fst.old[i] = fst.new[i]
          }
        }
      }
    }
    # Required Fst is found and corresponding allele frequencies are specified
    for (j in 1:limit[i]) {
      allele.freq.1[i, j] = x.old.1[i, j]
      allele.freq.2[i, j] = x.old.2[i, j]
    }
  }
  print("Simulated Fst values for each locus")
  print(data.frame(fst.new))
  #
  # Simulated populations of multicolour genotypes are generated according to
  # these allele frequencies.
  N.O = round((1 - c) * N / 2)   #  Proportion of simulated population non-inbred
  N.F = round(c * N / 2)         #  Proportion of simulated population inbred
  #
  #
  #
  # Generate population number 1
  # Generate homozygotes in inbred population(homo.F) and non-inbred(homo.O)
  for (i in 1:num.loc) {
    for (j in 1:limit[i]) {
      homo.1.F[i, j] = allele.freq.1[i, j] * (r.actual + (1 - r.actual) * allele.freq.1[i, j])
      homo.1.O[i, j] = (allele.freq.1[i, j]) ^ 2
    }
    # Generate heterozygous genotypes
    for (j in 1:(limit[i] - 1)) {
      for (k in (j + 1):limit[i]) {
        het.1.F[i, j, k] = 2 * allele.freq.1[i, j] * allele.freq.1[i, k] * (1 - r.actual)
        het.1.O[i, j, k] = 2 * allele.freq.1[i, j] * allele.freq.1[i, k]
      }
    }
    # Generate cumulative sum. This allows a random number to “pick” either a homozygous genotype
    # or a heterozygous genotype according to the relative frequency of the genotype
    for (j in 1:limit[i]) {
      homo.cumulative.1.F[i, j] = sum(homo.1.F[i, 1:j])
      homo.cumulative.1.O[i, j] = sum(homo.1.O[i, 1:j])
    }
    het.cumulative.1.F.temp = homo.cumulative.1.F[i, limit[i]]
    het.cumulative.1.O.temp = homo.cumulative.1.O[i, limit[i]]
    for (j in 1:(limit[i] - 1)) {
      for (k in (j + 1):limit[i]) {
        het.cumulative.1.F[i, j, k] = het.1.F[i, j, k] + het.cumulative.1.F.temp
        het.cumulative.1.O[i, j, k] = het.1.O[i, j, k] + het.cumulative.1.O.temp
        het.cumulative.1.F.temp = het.cumulative.1.F[i, j, k]
        het.cumulative.1.O.temp = het.cumulative.1.O[i, j, k]
      }
    }
  }
  if (N.F > 0) {
    for (i in 1:N.F) {
      #inbred population
      for (j in 1:num.loc) {
        random.number = runif(1, 0, 1)
        check = FALSE
        for (k in 1:limit[j]) {
          if (random.number <= homo.cumulative.1.F[j, k] & check == FALSE) {
            genotype.1.allele.1.F[i, j] = k
            genotype.1.allele.2.F[i, j] = k
            all.count.1[j, k] = all.count.1[j, k] + 2
            check = TRUE
          }
        }
        for (l in 1:(limit[j] - 1)) {
          for (m in (l + 1):limit[j]) {
            if (random.number <= het.cumulative.1.F[j, l, m] & check == FALSE) {
              genotype.1.allele.1.F[i, j] = l
              genotype.1.allele.2.F[i, j] = m
              all.count.1[j, l] = all.count.1[j, l] + 1
              all.count.1[j, m] = all.count.1[j, m] + 1
              check = TRUE
            }
          }
        }
      }
    }
  }
  if (N.O > 0) {
    for (i in (2 * N.F + 1):(2 * N.F + N.O)) {
      #non-inbred population
      for (j in 1:num.loc) {
        random.number = runif(1, 0, 1)
        check = FALSE
        for (k in 1:limit[j]) {
          if (random.number <= homo.cumulative.1.O[j, k] & check == FALSE) {
            genotype.1.allele.1.O[i, j] = k
            genotype.1.allele.2.O[i, j] = k
            all.count.1[j, k] = all.count.1[j, k] + 2
            check = TRUE
          }
        }
        for (l in 1:(limit[j] - 1)) {
          for (m in (l + 1):limit[j]) {
            if (random.number <= het.cumulative.1.O[j, l, m] & check == FALSE) {
              genotype.1.allele.1.O[i, j] = l
              genotype.1.allele.2.O[i, j] = m
              all.count.1[j, l] = all.count.1[j, l] + 1
              all.count.1[j, m] = all.count.1[j, m] + 1
              check = TRUE
            }
          }
        }
      }
    }
  }
  #
  # Generate population number 2
  # Generate homozygotes in inbred population(homo.F) and non-inbred(homo.O)
  for (i in 1:num.loc) {
    for (j in 1:limit[i]) {
      homo.2.F[i, j] = allele.freq.2[i, j] * (r.actual + (1 - r.actual) * allele.freq.2[i, j])
      homo.2.O[i, j] = (allele.freq.2[i, j]) ^ 2
    }
    # Generate heterozygous genotypes
    for (j in 1:(limit[i] - 1)) {
      for (k in (j + 1):limit[i]) {
        het.2.F[i, j, k] = 2 * allele.freq.2[i, j] * allele.freq.2[i, k] * (1 - r.actual)
        het.2.O[i, j, k] = 2 * allele.freq.2[i, j] * allele.freq.2[i, k]
      }
    }
    # Generate cumulative sum
    for (j in 1:limit[i]) {
      homo.cumulative.2.F[i, j] = sum(homo.2.F[i, 1:j])
      homo.cumulative.2.O[i, j] = sum(homo.2.O[i, 1:j])
    }
    het.cumulative.2.F.temp = homo.cumulative.2.F[i, limit[i]]
    het.cumulative.2.O.temp = homo.cumulative.2.O[i, limit[i]]
    for (j in 1:(limit[i] - 1)) {
      for (k in (j + 1):limit[i]) {
        het.cumulative.2.F[i, j, k] = het.2.F[i, j, k] + het.cumulative.2.F.temp
        het.cumulative.2.O[i, j, k] = het.2.O[i, j, k] + het.cumulative.2.O.temp
        het.cumulative.2.F.temp = het.cumulative.2.F[i, j, k]
        het.cumulative.2.O.temp = het.cumulative.2.O[i, j, k]
      }
    }
  }
  if (N.F > 0) {
    for (i in (N.F + 1):(2 * N.F)) {
      #inbred population
      for (j in 1:num.loc) {
        random.number = runif(1, 0, 1)
        check = FALSE
        for (k in 1:limit[j]) {
          if (random.number <= homo.cumulative.2.F[j, k] & check == FALSE) {
            genotype.2.allele.1.F[i, j] = k
            genotype.2.allele.2.F[i, j] = k
            all.count.2[j, k] = all.count.2[j, k] + 2
            check = TRUE
          }
        }
        for (l in 1:(limit[j] - 1)) {
          for (m in (l + 1):limit[j]) {
            if (random.number <= het.cumulative.2.F[j, l, m] & check == FALSE) {
              genotype.2.allele.1.F[i, j] = l
              genotype.2.allele.2.F[i, j] = m
              all.count.2[j, l] = all.count.2[j, l] + 1
              all.count.2[j, m] = all.count.2[j, m] + 1
              check = TRUE
            }
          }
        }
      }
    }
  }
  if (N.O > 0) {
    for (i in (2 * N.F + N.O + 1):(2 * N.F + 2 * N.O)) {
      #non-inbred population
      for (j in 1:num.loc) {
        random.number = runif(1, 0, 1)
        check = FALSE
        for (k in 1:limit[j]) {
          if (random.number <= homo.cumulative.2.O[j, k] & check == FALSE) {
            genotype.2.allele.1.O[i, j] = k
            genotype.2.allele.2.O[i, j] = k
            all.count.2[j, k] = all.count.2[j, k] + 2
            check = TRUE
          }
        }
        for (l in 1:(limit[j] - 1)) {
          for (m in (l + 1):limit[j]) {
            if (random.number <= het.cumulative.2.O[j, l, m] & check == FALSE) {
              genotype.2.allele.1.O[i, j] = l
              genotype.2.allele.2.O[i, j] = m
              all.count.2[j, l] = all.count.2[j, l] + 1
              all.count.2[j, m] = all.count.2[j, m] + 1
              check = TRUE
            }
          }
        }
      }
    }
  }
  #
  #
  # Calculate allele frequencies from total, combined population.
  # These are the average allele frequencies
  #
  for (i in 1:num.loc) {
    for (j in 1:max.alleles) {
      all.freq[i, j] = (all.count.1[i, j] + all.count.2[i, j]) / (2 * N)
    }
  }
  #
  # Calculate maximum Fst
  homozygous.frequency = 0
  for (i in 1:num.loc) {
    for (j in 1:max.alleles) {
      if (all.freq[i, j] > 0) {
        homozygous.frequency = homozygous.frequency + (all.freq[i, j]) ^ 2
      }
    }
  }
  homozygous.frequency = homozygous.frequency / num.loc
  fst.max = homozygous.frequency / (1 - homozygous.frequency)
  #
  if (fst.max > 1) {
    fst.max = 1
  }
  writeLines("Maximum value of Fst = ")
  print(fst.max)
  #
  f.values = fst.max / f.resolution
  #
  #
  # ANALYSIS, similar to construct function
  #
  # Calculate probability of multilocus genotypes for range of f values
  # Calculates probability of genotype, if homozygous, given value of f
  substructure.homozygous = function(f, all.freq.1) {
    all.freq.1 * (f + (1 - f) * all.freq.1)
  }
  # Calculates probability of genotype, if heterozygous, given value of f
  substructure.heterozygous = function(f, all.freq.1, all.freq.2) {
    2 * all.freq.1 * all.freq.2 * (1 - f)
  }
  # Calculates probability of genotype, if homozygous, given value of f and r
  consanguinity.homozygous = function(f, r.consider, all.freq.1) {
    all.freq.1 * (r.consider + (1 - r.consider) * (f + (1 - f) * all.freq.1))
  }
  # Calculates probability of genotype, if heterozygous, given value of f and r
  consanguinity.heterozygous = function(f, r.consider, all.freq.1, all.freq.2) {
    2 * all.freq.1 * all.freq.2 * (1 - r.consider) * (1 - f)
  }
  #
  # Iterate through values of f and c to estimate maximum likelihood values given simulated genotypes
  f.count = 1
  f = 0
  repeat {
    c.count = 1
    c = 0
    repeat {
      likelihood.locus.s = matrix(0, N, num.loc)
      likelihood.locus.c = matrix(0, N, num.loc)
      if (N.F > 0) {
        # Calculate probability of genotypes for inbred proportion (c) and non-inbred proportion (1-c)
        # for both simulated populations
        # Inbred population number 1, if genotype is homozygous
        for (i in 1:N.F) {
          for (j in 1:num.loc) {
            if (genotype.1.allele.1.F[i, j] > 0) {
              if (genotype.1.allele.1.F[i, j] == genotype.1.allele.2.F[i, j]) {
                likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.1.allele.1.F[i, j]])
                likelihood.locus.c[i, j] = consanguinity.homozygous(f, r.consider, all.freq[j, genotype.1.allele.1.F[i, j]])
              }
            }
          }
        }
        # Inbred population number 2, if genotype is homozygous
        for (i in (N.F + 1):(2 * N.F)) {
          for (j in 1:num.loc) {
            if (genotype.2.allele.1.F[i, j] > 0) {
              if (genotype.2.allele.1.F[i, j] == genotype.2.allele.2.F[i, j]) {
                likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.2.allele.1.F[i, j]])
                likelihood.locus.c[i, j] = consanguinity.homozygous(f, r.consider, all.freq[j, genotype.2.allele.1.F[i, j]])
              }
            }
          }
        }
        # Inbred population number 1, if heterozygous
        for (i in 1:N.F) {
          for (j in 1:num.loc) {
            if (genotype.1.allele.1.F[i, j] > 0) {
              if (genotype.1.allele.1.F[i, j] != genotype.1.allele.2.F[i, j]) {
                likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.1.allele.1.F[i, j]], all.freq[j, genotype.1.allele.2.F[i, j]])
                likelihood.locus.c[i, j] = consanguinity.heterozygous(f, r.consider, all.freq[j, genotype.1.allele.1.F[i, j]], all.freq[j, genotype.1.allele.2.F[i, j]])
              }
            }
          }
        }
        # Inbred population number 2, if heterozygous
        for (i in (N.F + 1):(2 * N.F)) {
          for (j in 1:num.loc) {
            if (genotype.2.allele.1.F[i, j] > 0) {
              if (genotype.2.allele.1.F[i, j] != genotype.2.allele.2.F[i, j]) {
                likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.2.allele.1.F[i, j]], all.freq[j, genotype.2.allele.2.F[i, j]])
                likelihood.locus.c[i, j] = consanguinity.heterozygous(f, r.consider, all.freq[j, genotype.2.allele.1.F[i, j]], all.freq[j, genotype.2.allele.2.F[i, j]])
              }
            }
          }
        }
      }
      # Non-inbred population number 1, if homozygous
      if (N.O > 0) {
        for (i in (2 * N.F + 1):(2 * N.F + N.O)) {
          for (j in 1:num.loc) {
            if (genotype.1.allele.1.O[i, j] > 0) {
              if (genotype.1.allele.1.O[i, j] == genotype.1.allele.2.O[i, j]) {
                likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.1.allele.1.O[i, j]])
                likelihood.locus.c[i, j] = consanguinity.homozygous(f, r.consider, all.freq[j, genotype.1.allele.1.O[i, j]])
              }
            }
          }
        }
        # Non-inbred population number 2, if homozygous
        for (i in (2 * N.F + N.O + 1):(2 * N.F + 2 * N.O)) {
          for (j in 1:num.loc) {
            if (genotype.2.allele.1.O[i, j] > 0) {
              if (genotype.2.allele.1.O[i, j] == genotype.2.allele.2.O[i, j]) {
                likelihood.locus.s[i, j] = substructure.homozygous(f, all.freq[j, genotype.2.allele.1.O[i, j]])
                likelihood.locus.c[i, j] = consanguinity.homozygous(f, r.consider, all.freq[j, genotype.2.allele.1.O[i, j]])
              }
            }
          }
        }
        # Non-inbred population number 1, if heterozygous
        for (i in (2 * N.F + 1):(2 * N.F + N.O)) {
          for (j in 1:num.loc) {
            if (genotype.1.allele.1.O[i, j] > 0) {
              if (genotype.1.allele.1.O[i, j] != genotype.1.allele.2.O[i, j]) {
                likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.1.allele.1.O[i, j]], all.freq[j, genotype.1.allele.2.O[i, j]])
                likelihood.locus.c[i, j] = consanguinity.heterozygous(f, r.consider, all.freq[j, genotype.1.allele.1.O[i, j]], all.freq[j, genotype.1.allele.2.O[i, j]])
              }
            }
          }
        }
        # Non-inbred population number 2, if heterozygous
        for (i in (2 * N.F + N.O + 1):(2 * N.F + 2 * N.O)) {
          for (j in 1:num.loc) {
            if (genotype.2.allele.1.O[i, j] > 0) {
              if (genotype.2.allele.1.O[i, j] != genotype.2.allele.2.O[i, j]) {
                likelihood.locus.s[i, j] = substructure.heterozygous(f, all.freq[j, genotype.2.allele.1.O[i, j]], all.freq[j, genotype.2.allele.2.O[i, j]])
                likelihood.locus.c[i, j] = consanguinity.heterozygous(f, r.consider, all.freq[j, genotype.2.allele.1.O[i, j]], all.freq[j, genotype.2.allele.2.O[i, j]])
              }
            }
          }
        }
      }
      # Take product of individual genotype probabilities and set limit to smallest probability (here 10^-100)
      for (i in 1:N) {
        if (prod(likelihood.locus.s[i,] > 0)) {
          likelihood.individual.s[i] = prod(likelihood.locus.s[i,])
        } else {
          likelihood.individual.s[i] = 10 ^ (-100)
        }
        if (prod(likelihood.locus.c[i,] > 0)) {
          likelihood.individual.c[i] = prod(likelihood.locus.c[i,])
        } else{
          likelihood.individual.c[i] = 10 ^ (-100)
        }
        # Take log of probability of individual multilocus genotypes according to value of c
        ln.likelihood.individual[i] = log(c * likelihood.individual.c[i] + (1 - c) * likelihood.individual.s[i])
      }
      # Take sum of probabilities of individual multilocus genotypes
      ln.likelihood.cs[f.count, c.count] = sum(ln.likelihood.individual)
      c.count = c.count + 1
      c = c + 1 / c.resolution
      if (c.count == c.limit) {
        break
      }
    }
    f.count = f.count + 1
    
    f = f + f.values
    if (f.count == f.limit) {
      break
    }
  }
  # Find maximum likelihood by subtracting maximum likelihood value from each likelihood value
  max.value = max(ln.likelihood.cs)
  temp = matrix(0, f.resolution, c.resolution)
  for (i in 1:f.resolution) {
    for (j in 1:c.resolution) {
      temp[i, j] = ln.likelihood.cs[i, j] - max.value
    }
  }
  for (i in 1:f.resolution) {
    for (j in 1:c.resolution) {
      ln.likelihood.cs[i, j] = temp[i, j]
    }
  }
  for (i in 1:f.resolution) {
    for (j in 1:c.resolution) {
      if (ln.likelihood.cs[i, j] == 0) {
        ML.F = ((i - 1) / f.resolution) * fst.max
        ML.C = (j / c.resolution) - 0.01
      }
    }
  }
  writeLines("Maximum Likelihood value: ")
  writeLines("Fst=")
  print(ML.F)
  writeLines("C = ")
  print(ML.C)
  #
  #
  # Convert log values
  for (i in 1:f.resolution) {
    for (j in 1:c.resolution) {
      e.likelihood.cs[i, j] = exp(ln.likelihood.cs[i, j])
    }
  }
  total = sum(e.likelihood.cs)
  for (i in 1:f.resolution) {
    for (j in 1:c.resolution) {
      e.likelihood.cs[i, j] = e.likelihood.cs[i, j] / total
    }
  }
  write.table(e.likelihood.cs, file = "ConStruct.Sim.Outfile.txt")
  #
  x = seq(0, fst.max - f.values, by = f.values)
  y = seq(0, 0.99, by = 1 / c.resolution)
  contour(
    x,
    y,
    e.likelihood.cs,
    xlab = expression("F"[ST]),
    ylab = expression("c"[g]),
    nlevels = 4
  )
  writeLines("G = ")
  G = exp(log(max(e.likelihood.cs)) - 2)
  print(G)
  # Global variables FST(f.axis) cg(c.axis) and e.likelihood(probability)
  f.axis <<- x
  c.axis <<- y
  probability <<- e.likelihood.cs
}
# end of function