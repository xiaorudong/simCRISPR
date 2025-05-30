#' Simulate Amplification and Sequencing Process
#'
#' This function simulates the sequencing output by applying amplification and sequencing steps to sgRNA count data.
#' It models biases and variability introduced by PCR amplification, capture efficiency, and sequencing depth.
#'
#' @param wide_data A data frame of sgRNA counts before amplification, typically from a prior simulation step.
#' @param rounds Number of PCR amplification rounds to simulate. Default is 5.
#' @param totalDepth Total sequencing depth to distribute across sgRNAs. Default is \code{1e9}.
#' @param cm_mu Mean of the capture efficiency \eqn{\mu_{\text{cm}}} used in generating \eqn{\mu_{\text{PCR}}}.
#' @param cm_sd Standard deviation \eqn{\sigma_{\text{cm}}} of the capture efficiency.
#' @param ep_mu Mean of the error propagation rate during PCR.
#' @param ep_sd Standard deviation of the error propagation rate.
#' @param pcr_sd Standard deviation \eqn{\sigma_{\text{PCR}}} used in the diagonal covariance matrix \eqn{\Sigma_{\text{PCR}}}.
#' @param sf_sd Standard deviation for sampling fluctuation during the sequencing process.
#'
#' @return A data frame of simulated read counts after sequencing, incorporating amplification and sequencing biases.
#' @export
#'

seq_add <- function(wide_data, rounds=5, totalDepth = 1e9,
                    cm_mu=0.99, cm_sd=0.01,
                    ep_mu=0.01, ep_sd=0.001,
                    pcr_sd=0.02, sf_sd=0.03) {

  afteramp <- simCRISPR::amplifyStep(wide_data, rounds=rounds,
                          cm_mu=cm_mu, cm_sd=cm_sd,
                          ep_mu=ep_mu, ep_sd=ep_sd,
                          pcr_sd=pcr_sd)

  combined_data <- as.data.frame(simCRISPR::sequenceStep(amp_frags=afteramp$amp_frags, totalDepth = totalDepth, sf_sd=sf_sd))

  return(combined_data)
}


#' Simulate PCR Amplification of sgRNA Fragments
#'
#' This function simulates the PCR amplification process applied to sgRNAs during sequencing library preparation.
#'
#' @param capturedMolecules A data frame or matrix of initial sgRNA molecule counts before amplification.
#' @param rounds Number of PCR amplification rounds to simulate. Default is 5.
#' @param cm_mu Mean of the capture efficiency \eqn{\mu_{\text{cm}}} used in generating \eqn{\mu_{\text{PCR}}}.
#' @param cm_sd Standard deviation \eqn{\sigma_{\text{cm}}} of the capture efficiency.
#' @param ep_mu Mean of the error propagation rate during PCR.
#' @param ep_sd Standard deviation of the error propagation rate.
#' @param pcr_sd Standard deviation \eqn{\sigma_{\text{PCR}}} used in the diagonal covariance matrix \eqn{\Sigma_{\text{PCR}}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{amp_frags}}{A matrix of PCR-amplified fragment counts for each sgRNA and sample.}
#'   \item{\code{eff_matrix}}{The amplification efficiency matrix used in the simulation.}
#' }
#' @importFrom stats rnorm
#' @importFrom MASS mvrnorm
#'
amplifyStep <- function(capturedMolecules, rounds,
                        cm_mu=0.99, cm_sd=0.01,
                        ep_mu=0.01, ep_sd=0.001,
                        pcr_sd=0.02){

  X <- as.matrix(capturedMolecules)

  # simulate the PCR efficiency from a multivariate normal distribution
  sigma2 <- matrix(0, nrow = nrow(capturedMolecules), ncol = nrow(capturedMolecules), byrow = TRUE)
  diag(sigma2) <- rep(pcr_sd^2, nrow(capturedMolecules))
  mu <- stats::rnorm(nrow(capturedMolecules), cm_mu, cm_sd)

  efficiencyPCR_mat <- t(MASS::mvrnorm(n=ncol(capturedMolecules), mu = mu, Sigma = sigma2))
  # head(efficiencyPCR_mat)
  # plot(efficiencyPCR_mat[,1:2])
  #

  # considering preferential amplification
  efficiencyPrefersg <- stats::rnorm(round(nrow(efficiencyPCR_mat)/10), ep_mu, ep_sd)
  pickPrefersg <- rep(0, nrow(efficiencyPCR_mat))
  pickPrefersg[sample(1:nrow(efficiencyPCR_mat), length(efficiencyPrefersg))] <- efficiencyPrefersg
  names(pickPrefersg) <- rownames(capturedMolecules)

  efficiencyPCR_mat <- pickPrefersg+efficiencyPCR_mat
  efficiencyPCR_mat <- pmin(efficiencyPCR_mat, 1 - 1e-7)


  amp_mat <- (1 + efficiencyPCR_mat)^rounds


  A <- X*amp_mat
  A <- round(A, digits = 0)

  return(list(amp_frags=A, prefersg=pickPrefersg[pickPrefersg!=0]))
}


#' Simulate Sequencing from Amplified Fragments
#'
#' This function simulates the sequencing output from PCR-amplified sgRNA fragments.
#' After amplification, sequencing introduces additional variability due to technical noise.
#'
#' @param amp_frags A matrix or data frame of PCR-amplified fragment counts for each sgRNA (e.g., the output from \code{amplifyStep}).
#' @param totalDepth Total sequencing depth to allocate across all sgRNAs. Default is \code{1e9}.
#' @param sf_sd Standard deviation of the sampling fluctuation applied during sequencing.
#'
#' @return A data frame of simulated read counts for each sgRNA after sequencing, incorporating stochastic sampling noise.
#' @importFrom stats rnorm rmultinom
#'
sequenceStep <- function(amp_frags, totalDepth, sf_sd=0.03) {

  geneProbs_all <- sweep(amp_frags, MARGIN=2, colSums(amp_frags), `/`)
  # equalizeCells: combine everything in equal proportion.
  totalM <- colSums(amp_frags)
  target <- min(totalM)
  useMean <- target / totalM
  SF <- useMean + stats::rnorm(length(useMean), 0, sd = sf_sd)
  SF[SF>1] <- 1
  SF[SF<0] <- 1e-4
  totalM <- totalM * SF

  inputRange <- totalM
  if (any(inputRange >= .Machine$integer.max)) { # largest value in R, just rescaling.
    SCALEALL <- max(inputRange) / .Machine$integer.max
    inputRange <- inputRange/SCALEALL
  }

  afterEq <- sapply(seq_len(ncol(geneProbs_all)), function(ind) {
    stats::rmultinom(n=1, size=inputRange[ind], prob=geneProbs_all[,ind])
  })

  # mix all the samples, give unique gene_sample names
  Probs_all <- afterEq/sum(afterEq)
  # considering efficiency factors
  efficiencyFactors <- stats::rnorm(nrow(Probs_all), mean=.99, sd=0.01)
  efficiencyFactors <- pmin(efficiencyFactors, 1)
  adjustedProbs_all <- Probs_all * efficiencyFactors

  # begin sequencing
  mycrispr_vec <- stats::rmultinom(n=1, size=totalDepth, prob=as.vector(adjustedProbs_all))
  mycrispr <- matrix(mycrispr_vec, nrow = nrow(adjustedProbs_all), byrow = F)

  colnames(mycrispr) <- colnames(amp_frags)
  rownames(mycrispr) <- rownames(amp_frags)
  return(mycrispr)
}
