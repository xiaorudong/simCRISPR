#' Generate Constant Random Noise Across sgRNAs
#'
#' This function generates a constant random noise value drawn from a normal distribution with mean 0 and standard deviation `SD`,
#' and repeats it across all sgRNAs. It is used to simulate global or shared perturbations across all guides.
#'
#' @param SD Standard deviation of the noise.
#' @param n_total Total number of sgRNAs to apply the noise to.
#'
#' @return A numeric vector of length `n_total` with identical noise values.
#' @importFrom stats rnorm
BLEH <- function(SD, n_total) {
  return( rep(stats::rnorm(1, mean = 0, sd=SD), n_total))
}

#' Generate Independent Random Noise for Each sgRNA
#'
#' This function simulates sgRNA-specific noise by drawing independent values from a normal distribution with mean 0
#' and specified standard deviation.
#'
#' @param sd Standard deviation of the noise.
#' @param n_total Total number of sgRNAs to generate noise for.
#'
#' @return A numeric vector of length `n_total`, each element independently drawn from a normal distribution.
#' @importFrom stats rnorm
#'
sim_noise <- function(sd=0.5, n_total){
  return(stats::rnorm(n_total, mean = 0, sd))
}





#' Simulate Pooled CRISPR Screen Data
#'
#' This function generates synthetic data for pooled CRISPR screening experiments under various experimental conditions.
#' It simulates counts by modeling sgRNA knockout effects, treatment effects, interaction effects, sequencing noise,
#' and disturbances such as DNA damage responses. The function supports both exponential and logistic cell growth models.
#'
#' @param method Growth model used to simulate cell population. Options are `"logit"`, `"exp"`, or `"both"`. Default is `"exp"`. `"logit"` uses the logistic growth formula, modeling saturation effects as cell populations approach a carrying capacity. `"exp"` uses exponential growth, assuming unlimited resources and constant growth rate. `"both"` generates output under both logistic and exponential growth models for comparison.
#' @param samples Experimental design for sample generation. Options are `"independent"` or `"matched"`. Default is `"independent"`. `"independent"` means the control and treatment samples are generated independently. `"matched"` means the control and treatment samples are paired, simulating repeated measurements or shared sample origin.
#' @param L Carrying capacity used in the logistic growth model.
#' @param trt_eff_mu Mean of the treatment effect.
#' @param trt_eff_sd Standard deviation of the treatment effect.
#'
#' @param sg_eff_prop1 Proportion of sgRNAs with a primary knockout effect.
#' @param sg_eff_mu1 Mean of the primary sgRNA knockout effect.
#' @param sg_eff_sd1 Standard deviation of the primary sgRNA knockout effect.
#' @param sg_eff_mu2 Mean of the remaining sgRNA knockout effects (used when not assigned to sg_eff_mu1). If not provided, defaults to -1*`sg_eff_mu1`.
#' @param sg_eff_sd2 Standard deviation of the remaining sgRNA knockout effects. If not provided, defaults to `sg_eff_sd1`.
#'
#' @param distbDNA_mu Mean of the DNA damage response effect.
#' @param distbDNA_sd Standard deviation of the DNA damage response effect.
#'
#' @param sg_trt_eff_prop1 Proportion of sgRNAs with interaction effects between sgRNA and treatment from the primary group.
#' @param sg_trt_eff_mu1 Mean of the sgRNA-treatment interaction effect for the primary group.
#' @param sg_trt_eff_sd1 Standard deviation of the sgRNA-treatment interaction effect for the primary group.
#' @param sg_trt_eff_mu2 Mean of the remaining sgRNA-treatment interaction effects. If not provided, defaults to -1*`sg_trt_eff_mu1`.
#' @param sg_trt_eff_sd2 Standard deviation of the remaining sgRNA-treatment interaction effects. If not provided, defaults to `sg_trt_eff_sd1`.
#' @param prop_sg_trt Proportion of sgRNAs under treatment conditions.
#'
#' @param distb_trt_eff_prop1 Proportion of sgRNAs with DNA damage-treatment interaction effects from the primary group.
#' @param distb_trt_eff_mu1 Mean of the DNA damage-treatment interaction effect (primary group).
#' @param distb_trt_eff_sd1 Standard deviation of the DNA damage-treatment interaction effect (primary group).
#' @param distb_trt_eff_mu2 Mean of the remaining DNA damage-treatment interaction effects.
#' @param distb_trt_eff_sd2 Standard deviation of the remaining DNA damage-treatment interaction effects.
#' @param prop_distb_trt Proportion of sgRNAs assigned DNA damage-treatment interaction effects.
#'
#' @param bleh_SD Standard deviation of random noise across different conditions or replicates.
#' @param noise_SD Standard deviation of additional technical noise applied to knockout samples.
#' @param days Number of time points (e.g., days) across which the experiment is simulated.
#' @param reps Number of biological replicates per condition.
#' @param n_total Total number of sgRNAs, including knockout, non-targeting, and safe harbor controls.
#' @param n_ntgt Number of non-targeting control sgRNAs.
#' @param n_sfhb Number of safe harbor control sgRNAs.
#' @param initial_mu Mean initial abundance of each sgRNA at the beginning of the experiment.
#'
#' @return A list containing simulated read count matrices and the corresponding true effects for each sgRNA.
#' The structure of the returned list depends on the specified `method`:
#' \itemize{
#'   \item If \code{method = "logit"}: returns \code{sim_raw}, \code{true_eff_t}, and \code{sim_full} for logistic growth.
#'   \item If \code{method = "exp"}: returns \code{sim_raw}, \code{true_eff_t}, and \code{sim_full} for exponential growth.
#'   \item If \code{method = "both"}: returns \code{sim_logit_raw}, \code{sim_exp_raw}, \code{true_eff_t}, \code{sim_logit_full}, and \code{sim_exp_full}.
#' }
#'
#' \describe{
#'   \item{\code{sim_raw}, \code{sim_logit_raw}, \code{sim_exp_raw}}{Raw counts from knockout samples only.}
#'   \item{\code{true_eff_t}}{Matrix of true effect values for each sgRNA under different conditions (e.g., knockout, treatment, interaction).}
#'   \item{\code{sim_full}, \code{sim_logit_full}, \code{sim_exp_full}}{Full raw count matrix including initial sgRNA counts, control samples, and knockout samples.}
#' }
#' @import FamilyRank
#' @importFrom stats rnorm

sim_crispr <- function(method = "exp",
                       samples = "independent",
                       L=10000,
                       trt_eff_mu= -0.5,
                       trt_eff_sd= 0.01,

                       sg_eff_prop1 = 0.5, ##Set large KO effect
                       sg_eff_mu1 = -0.3,
                       sg_eff_sd1 = 0.1,
                       sg_eff_mu2 = NULL,
                       sg_eff_sd2 = NULL,

                       distbDNA_sd=0.02,
                       distbDNA_mu=-0.01,

                       sg_trt_eff_prop1 = 0.7,
                       sg_trt_eff_mu1 = -0.2,
                       sg_trt_eff_sd1 = 0.1,
                       sg_trt_eff_mu2 = NULL,
                       sg_trt_eff_sd2 = NULL,
                       prop_sg_trt=0.4,

                       distb_trt_eff_prop1 = 1,
                       distb_trt_eff_mu1 = -0.001,
                       distb_trt_eff_sd1 = 0.1,
                       distb_trt_eff_mu2 = NULL,
                       distb_trt_eff_sd2 = NULL,
                       prop_distb_trt=0.4,

                       bleh_SD=0.1,
                       noise_SD=0.005,
                       days=5,
                       reps=3,
                       n_total=1000,
                       n_ntgt=100,
                       n_sfhb=50,
                       initial_mu=1000) {
  #### Try to distinguish non targetting vs safe harbor (safe harbor make still have a "KO" effect due to cutting hurting the cell's overall fitness)

  valid_methods <- c("logit", "exp", "both")
  if (!method %in% valid_methods) {stop("Invalid method specified. Use 'logistic' or 'exponential'")}

  valid_samples <- c("independent", "matched")
  if (!samples %in% valid_samples) {stop("Invalid samples specified. Use 'independent' or 'matched'")}


  #################################
  ###### Simulation Settings ######
  #################################

  nsgRNA <- n_total-n_ntgt-n_sfhb

  # original growth rate
  orig_gr <- rep(1, n_total) # growth rate

  # initial cells are from a binomial distribution
  # y0_binom_prob <- 0.5
  # y0 <- rbinom(n_total, initial_mu, prob = y0_binom_prob)
  y0 <- round(stats::rnorm(n_total, initial_mu, initial_mu/10))
  if (any(y0 < 10)) {
    warning("Warning: y0 contains values less than 10, setting them to 10.")
    y0[y0 < 10] <- 10
  }

  # knock out efficiency based on beta distribution
  # ko_eff_mode <- 0.8; ko_eff_shape2 <- 5
  # ko_eff <- rbeta(nsgRNA, shape1 = (ko_eff_mode*(ko_eff_shape2-2)+1)/(1-ko_eff_mode), shape2 = ko_eff_shape2)
  # hist(ko_eff)
  ko_eff <- c(rep(0, n_ntgt), rep(1, n_sfhb), rep(1, nsgRNA))

  # sgRNA effects on sgRNA only
  if(is.null(sg_eff_mu2)) sg_eff_mu2 <- (-1)*sg_eff_mu1
  if(is.null(sg_eff_sd2)) sg_eff_sd2 <- sg_eff_sd1
  sg_eff_prop2 <- 1-sg_eff_prop1
  sg_eff <- FamilyRank::rbinorm(n=nsgRNA, mean1=sg_eff_mu1, mean2=sg_eff_mu2, sd1=sg_eff_sd1,
                                sd2=sg_eff_sd2, prop = sg_eff_prop1)
  sg_eff <- c(rep(0, n_ntgt+n_sfhb), sg_eff)

  # DNA disturb effect
  disturb_eff <- c(rep(0, n_ntgt), stats::rnorm(nsgRNA+n_sfhb, mean=distbDNA_mu, sd=distbDNA_sd))

  # treatment or toxin effects # negative value indicates a toxin
  trt_eff <- stats::rnorm(n_total, mean = trt_eff_mu, sd=trt_eff_sd)

  # interactive effect between sgRNA and treatment/toxin from a bimodal distribution
  if(is.null(sg_trt_eff_mu2)) sg_trt_eff_mu2 <- (-1)*sg_trt_eff_mu1
  if(is.null(sg_trt_eff_sd2)) sg_trt_eff_sd2 <- sg_trt_eff_sd1
  sg_trt_eff_prop2 <- 1-sg_trt_eff_prop1
  sg_trt_eff <- FamilyRank::rbinorm(n=nsgRNA, mean1=sg_trt_eff_mu1, mean2=sg_trt_eff_mu2,
                                    sd1=sg_trt_eff_sd1, sd2=sg_trt_eff_sd2, prop = sg_trt_eff_prop1)
  I_sg_trt <- sample(c(0, 1), size = nsgRNA, replace = T, prob = c(1-prop_sg_trt, prop_sg_trt))
  sg_trt_eff <- I_sg_trt*sg_trt_eff
  sg_trt_eff <- c(rep(0, n_ntgt+n_sfhb), sg_trt_eff)


  # interactive effect between distrub DNA and treatment/toxin from a bimodal distribution
  if(is.null(distb_trt_eff_mu2)) distb_trt_eff_mu2 <- (-1)*distb_trt_eff_mu1
  if(is.null(distb_trt_eff_sd2)) distb_trt_eff_sd2 <- distb_trt_eff_sd1
  distb_trt_eff_prop2 <- 1-distb_trt_eff_prop1
  distb_trt_eff <- FamilyRank::rbinorm(n=n_sfhb+nsgRNA, mean1=distb_trt_eff_mu1, mean2=distb_trt_eff_mu2,
                                       sd1=distb_trt_eff_sd1, sd2=distb_trt_eff_sd2, prop = distb_trt_eff_prop1)
  I_distb_trt <- sample(c(0, 1), size = n_sfhb+nsgRNA, replace = T, prob = c(1-prop_distb_trt, prop_distb_trt))
  distb_trt_eff <- I_distb_trt*distb_trt_eff
  distb_trt_eff <- c(rep(0, n_ntgt), distb_trt_eff)


  #########################
  ###### True Effect ######
  #########################

  true_eff <- as.data.frame(rbind(sg_eff*ko_eff, trt_eff, ko_eff*sg_trt_eff, disturb_eff, distb_trt_eff))
  colnames(true_eff)[1:n_ntgt] <- paste("ntgt", 1:n_ntgt, sep = "")
  colnames(true_eff)[(n_ntgt+1):(n_ntgt+n_sfhb)] <- paste("sfhb", 1:n_sfhb, sep = "")
  colnames(true_eff)[(n_ntgt+n_sfhb+1):n_total] <- paste("sg", 1:nsgRNA, sep = "")
  rownames(true_eff) <- c("KO", "TRT", "INT", "Distb", "DistbINT")

  true_eff_t <- as.data.frame(t(true_eff))
  true_eff_t$sgID <- rownames(true_eff_t)


  ##############################
  ###### Start Simulation ######
  ##############################

  contrast_mat <- matrix(data = c(0, 0, 0, 0, 0,
                                  0, 1, 0, 0, 0,
                                  1, 0, 0, 1, 0,
                                  1, 1, 1, 1, 1),
                         nrow = 4, ncol=5, byrow = T)
  rownames(contrast_mat) <- c("ctrl","ctrl_trt", "ko", "ko_trt")
  colnames(contrast_mat) <- c("sgRNA_effect", "trt_effect", "interaction",
                              "distb_effect", "distb_int")


  y0_noise <- function(SD=.1) {
    return(stats::rnorm(length(y0), 0, SD))
  }

  # number of reps
  gr_reps <- list()
  y0_reps <- list()
  ytX_reps <- list()
  ytX_lgst_reps <- list()
  for (i_reps in 1:(reps*2)) {
    extra_gr <- as.matrix(cbind(sg_eff*ko_eff+simCRISPR::BLEH(bleh_SD, n_total),
                                trt_eff+simCRISPR::BLEH(bleh_SD, n_total),
                                ko_eff*sg_trt_eff+simCRISPR::BLEH(bleh_SD, n_total),
                                disturb_eff+simCRISPR::BLEH(bleh_SD, n_total),
                                distb_trt_eff+simCRISPR::BLEH(bleh_SD, n_total))) %*% t(contrast_mat)
    extra_gr[,1] <- sim_noise(sd=noise_SD, n_total)
    extra_gr[,3] <- extra_gr[,3]+sim_noise(sd=noise_SD, n_total)
    gr <- orig_gr + extra_gr

    y0_rep <- floor(y0+10*y0_noise())

    if (method == "exp" | method == "both") {
      ytX <- as.data.frame(floor(y0_rep*exp(gr*days)))
      colnames(ytX) <- paste(colnames(ytX), "_rep", i_reps, sep = "")
      ytX_reps[[i_reps]] <- ytX
    }

    if (method == "logit" | method == "both") {
      ytX_lgst <- as.data.frame(floor(L/(1+((L/y0_rep)-1)*exp((-1)*gr*days))))
      colnames(ytX_lgst) <- paste(colnames(ytX_lgst), "_rep", i_reps, sep = "")
      ytX_lgst_reps[[i_reps]] <- ytX_lgst
    }
    gr_reps[[i_reps]] <- gr

    y0_info <- as.data.frame(y0_rep)
    colnames(y0_info) <- paste("y0_rep", i_reps, sep = "")
    y0_reps[[i_reps]] <- y0_info
  }


  if(samples=="independent") {
    picky0reps <- sample(1:(reps*2), reps)
    pickctrlreps <- sample(1:(reps*2), reps)
    pickctrltrtreps <- sample(1:(reps*2), reps)
    pickkoreps <- sample(1:(reps*2), reps)
    pickkotrtreps <- sample(1:(reps*2), reps)

    picky0cols <- paste("y0_rep", picky0reps, sep = "")
    pickallcols <- c(paste("ctrl_rep", pickctrlreps, sep = ""), paste("ctrl_trt_rep", pickctrltrtreps, sep = ""),
                     paste("ko_rep", pickkoreps, sep = ""), paste("ko_trt_rep", pickkotrtreps, sep = ""))
  }

  if(samples=="matched") {
    pickreps <- sample(1:(reps*2), reps)

    picky0cols <- paste("y0_rep", pickreps, sep = "")
    pickallcols <- c(paste("ctrl_rep", pickreps, sep = ""), paste("ctrl_trt_rep", pickreps, sep = ""),
                     paste("ko_rep", pickreps, sep = ""), paste("ko_trt_rep", pickreps, sep = ""))
  }

  y0s_df <- as.data.frame(do.call(cbind, y0_reps))

  if (method == "exp" | method == "both") {
    raw_exp_count_wide <- do.call(cbind, ytX_reps)
    rownames(raw_exp_count_wide) <- rownames(true_eff_t)

    raw_exp_count_wide <- raw_exp_count_wide[,c(grep("ctrl_rep", colnames(raw_exp_count_wide), fixed = T),
                                                grep("ctrl_trt_rep", colnames(raw_exp_count_wide), fixed = T),
                                                grep("ko_rep", colnames(raw_exp_count_wide), fixed = T),
                                                grep("ko_trt_rep", colnames(raw_exp_count_wide), fixed = T))]
    raw_exp_count_wide <- as.data.frame(raw_exp_count_wide)
    raw_exp_count_wide <- as.data.frame(cbind(y0s_df[,picky0cols], raw_exp_count_wide[,pickallcols]))

    sub_exp_wide_data <- raw_exp_count_wide[,grep("ko", colnames(raw_exp_count_wide))]
  }

  if (method == "logit" | method == "both"){
    raw_lgst_count_wide <- do.call(cbind, ytX_lgst_reps)
    rownames(raw_lgst_count_wide) <- rownames(true_eff_t)

    raw_lgst_count_wide <- raw_lgst_count_wide[,c(grep("ctrl_rep", colnames(raw_lgst_count_wide), fixed = T),
                                                  grep("ctrl_trt_rep", colnames(raw_lgst_count_wide), fixed = T),
                                                  grep("ko_rep", colnames(raw_lgst_count_wide), fixed = T),
                                                  grep("ko_trt_rep", colnames(raw_lgst_count_wide), fixed = T))]
    raw_lgst_count_wide <- as.data.frame()
    raw_lgst_count_wide <- as.data.frame(cbind(y0s_df[,picky0cols], raw_lgst_count_wide[,pickallcols]))

    sub_lgst_wide_data <- raw_lgst_count_wide[,grep("ko", colnames(raw_lgst_count_wide))]
  }



  if (method == "logit") return(list(sim_raw=sub_lgst_wide_data,true_eff_t = true_eff_t,
                                     sim_full=raw_lgst_count_wide))
  if (method == "exp") return(list(sim_raw=sub_exp_wide_data,true_eff_t = true_eff_t,
                                   sim_full=raw_exp_count_wide))
  if (method == "both") return(list(sim_logit_raw=sub_lgst_wide_data, sim_exp_raw=sub_exp_wide_data, true_eff_t = true_eff_t,
                                    sim_logit_full=raw_lgst_count_wide, sim_exp_full=raw_exp_count_wide))
}

