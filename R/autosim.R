#' Auto-simulate the records/traits for plant breeding with AlphaSimR
#'
#' This function would read in an edibble experiment (that must include record factors)
#' and infers the context to determine what the appropriate values or distributions
#' that the record factor takes place using the power of large language models (LLM).
#' This function is similar in spirit with `edibble::autofill_rcrds()`, however, it is
#' purposely built for simulating plant breeding data with `AlphaSimR` and assumes all
#' records/traits are quantitative and have a genetic basis.
#' Based on the context given LLM will infer the parameters for the simulation scheme.
#'
#'
#' ## Limiations
#'
#' Currently all selections were assumed to be phenotypic.
#'
#' This was a prototype implementation for proof of concept, so you should not use this for serious work.
#'
#'
#'
#' @param design An `edibble` design object.
#' @param chat A Chat object from `ellmer`.
#' @param context Further context
#' @param seed Seed number for reproducibility of the simulation (not simulation scheme).
#' @param nsim Number of replicate for simulation. Default is 1.
#'
#' @export
simulate_plant_experiment <- function(design = NULL,
                                      chat = ellmer::chat_openai(),
                                      context = NULL,
                                      seed = NULL,
                                      nsim = 1L) {
  prov <- activate_provenance(design)
  rcrds <- prov$rcrd_names()

  # functions
  ask_yn <- function(question, env = rlang::caller_env()) {
    question <- glue_c(question, .content_envir = env)
    chatyn <- chat$clone()
    chatyn$chat_structured(question, type = ellmer::type_boolean())
  }

  ask_value <- function(question, env = rlang::caller_env()) {
    question <- glue_c(question, .content_envir = env)
    chatval <- chat$clone()
    chatval$chat_structured(question, type = ellmer::type_number())
  }

  ask_factor <- function(question, env = rlang::caller_env()) {
    question <- glue_c(question, .content_envir = env)
    chatfct <- chat$clone()
    out <- chatfct$chat_structured(question, type = ellmer::type_array(items = ellmer::type_string()))
    unlist(out)
  }

  # must have genotype as treatment -----------------------------------------
  geno <- chat$chat_structured(glue_c("What is the name of the variable with genotype in this experiment?"),
                               type = ellmer::type_enum(values = prov$trt_names()))
  n_geno <- length(prov$fct_levels(name = geno)[[1]])
  n <- nrow(design)


  # must have record objects ------------------------------------------------
  if(length(rcrds) == 0) {
    cli::cli_warn("There are no record objects in the design object.")
    return(design)
  }


  # pass the experimental context to LLM ------------------------------------
  chat$chat(paste0("Note that the following experimental context should be used to answer the questions from now on: ", describe_context(design)))
  if(!is.null(context)) chat$chat(context)


  # add expected record type -----------------------------------------------------
  current <- prov$get_validation()$rcrds
  # don't overwrite existing validation schemes
  rcrds_new <- setdiff(rcrds, names(current))

  ask_rcrd_type <- function(rcrd, type) {
    q <- glue_c("Given the experimental context, is {rcrd} commonly recorded as {type} only?")
    ask_yn(q)
  }

  get_rcrd_type <- function(rcrd) {
    ret <- list(type = NA, lower = NA, upper = NA, values = NA)
    ret$type <- "numerical"
    int <- ask_rcrd_type(rcrd, "integer")
    if(int) ret$type <- "integer"
    pnum <- ask_rcrd_type(rcrd, "a positive value")
    if(pnum) ret$lower <- 0
    out <- ask_value("What is the lower bound for {rcrd} if any?")
    if(!is.na(out)) ret$lower <- out
    out <- ask_value("What is the absolute upper bound for {rcrd} if any?")
    if(!is.na(out) && out > ret$lower) ret$upper <- out
    return(ret)
  }

  expect_new <- list()
  for(rcrd in rcrds_new) {
    type <- get_rcrd_type(rcrd)
    if(is.na(type$lower) & is.na(type$upper)) {
      wv <- with_value(between = c(-Inf, Inf))
    } else if(is.na(type$lower)) {
      wv <- with_value("<=", value = type$upper)
    } else if(is.na(type$upper)) {
      wv <- with_value(">=", value = type$lower)
    } else {
      wv <- with_value(between = c(type$lower, type$upper))
    }
    expect_new[[rcrd]] <- switch(type$type,
                                 "numerical" = to_be_numeric(wv),
                                 "integer" = to_be_integer(wv),
                                 "factor" = to_be_factor(levels = type$values))
  }
  rcrds_with_class <- map_chr(expect_new, function(x) x$record)
  design_with_expectations <- do.call("expect_rcrds", c(list(design), expect_new))


  # focus on numeric records ------------------------------------------------
  rcrds_numeric <- rcrds_with_class[rcrds_with_class %in% c("integer", "numeric")]
  rcrds_names <<- names(rcrds_numeric)
  if(length(rcrds_numeric) == 0) {
    cli::cli_warn("There are no numeric record objects in the design object.")
    return(autofill_rcrds(design_with_expectations))
  }

  # start simulation with AlphaSimR -----------------------------------------
  species <- chat$chat_structured("What species does this experiment involve?
                                   If it cannot be identified, return GENERIC.",
                                   type = ellmer::type_enum(values = c("GENERIC", "WHEAT", "MAIZE")))


  ## choose a breeding scheme
  cli::cli_alert_info("The identified species for this experiment is {.field {species}}.")

  ## specify global parameters
  n_chromosomes <- chat$chat_structured("What is the number of chromosomes for this species?",
                                         type = ellmer::type_integer())

  cli::cli_alert_info("Looks like there are {.field {n_chromosomes}} chromosomes.")

  # Using similar numbers to the Bancic et al (2025) paper
  # Not clear how these numbers are best chosen
  # Won't ask AI since I don't think this is common knowledge
  # Assumes these are the same across all traits
  # For other global parameters, it likely depends on traits so choose per trait
  n_sim_reps <- 1 # Number of simulation replicates
  n_burnin <- 20 # Number of years in burn-in phase
  n_future <- 20 # Number of years in future phase
  n_cycles <- n_burnin + n_future
  n_parents <- 50 # keep it tame for now so it's relatively quick
  n_crosses <- 100 # keep it tame for now
  # n_parents <- chat$chat_structured("What number of parents should we start a breeding cycle?",
  #                                   type = ellmer::type_integer())
  # n_crosses <- chat$chat_structured("How many number of crosses per year?",
  #                                   type = ellmer::type_integer())
  n_ploidy <- chat$chat_structured(glue::glue("What is the ploidy level of this species?"),
                                    type = ellmer::type_integer())
  cli::cli_alert_info("The identified species assumes a ploidy level of {.field {n_ploidy}}.")

  n_cohort <- 2 # number of cohorts - do not change this since downstream code is currently dependent on this
  info <- list(n_chromosomes = n_chromosomes,
               n_parents = n_parents,
               n_crosses = n_crosses,
               n_ploidy = n_ploidy,
               n_burnin = n_burnin,
               n_future = n_future,
               is_inbred = TRUE, # are all plants inbred... probably not, should ask AI
               n_cohort = n_cohort,
               n_geno = n_geno,
               n = n,
               geno = geno,
               species = species)

  # add simulation scheme to edibble ----------------------------------------
  design_sim <- design_with_expectations |>
    simulate_process(.joint = function(rcrds_names = rcrds_names,
                                       expect = expect,
                                       info = NULL,
                                       chat = NULL,
                                       initial = TRUE) {
      # this is needed because edibble evaluates the default initially to determine the records
      if(initial) {
        rnames <- get("rcrds_names", envir = .GlobalEnv)
        out <- lapply(rnames,
                      function(x) rnorm(n())) |>
          as.data.frame() |>
          setNames(rnames)
        return(out)
      }
      out <- list()
      for(rname in rcrds_names) {
        cli::cli_h2(rname)
        # specify global parameters
        n_qtls <- chat$chat_structured(glue::glue("What is the average number of QTLs for {rname} in {info$species}? The number should be typically much greater than one for quantitative traits."),
                                       type = ellmer::type_integer())
        n_snps <- chat$chat_structured(glue::glue("What is the average number of SNPs for {rname} in {info$species}? The number should be typically much greater than one for quantitative traits but presumably close to {n_qtls}."),
                                       type = ellmer::type_integer())
        cli::cli_alert_info("We expect an average of {.field {n_qtls}} QTLs and {.field {n_snps}} SNPs per chromosome.")
        mean_g <- chat$chat_structured(glue::glue("The mean phenotypic value for {rname}"),
                                        type = ellmer::type_number())
        var_g <- chat$chat_structured(glue::glue("The genetic variance for {rname}.
                                                 The result needs to be positive."),
                                        type = ellmer::type_number())
        var_ge <- chat$chat_structured(glue::glue("The total genotype-by-environment variance for {rname}.
                                                 The result needs to be positive."),
                                      type = ellmer::type_number())
        var_e <- chat$chat_structured(glue::glue("The environmental variance for {rname}.
                                                 The result needs to be positive."),
                                      type = ellmer::type_number())
        cli::cli_alert_info("Using the following parameters: mean = {.field {mean_g}}, var = {.field {var_g}}, varGxE = {.field {var_ge}}, and varE = {.field {var_e}}.")
        max_rep <- ceiling(info$n / info$n_geno)

        ## simulating genomes and founders
        founder_pop <-  AlphaSimR::runMacs(
          nInd     = info$n_parents,
          nChr     = info$n_chromosomes,
          segSites = n_qtls + n_snps,
          inbred   = info$is_inbred,
          ploidy   = info$n_ploidy,
          species  = info$species
        )
        SP <- AlphaSimR::SimParam$new(founder_pop)
        # Add SNP chip
        SP$restrSegSites(n_qtls, n_snps)
        if (n_snps > 0) {
          SP$addSnpChip(n_snps)
        }
        # Add traits for one or more additive GxE traits
        SP$addTraitAG(nQtlPerChr = n_qtls,
                      mean = rep(mean_g, max_rep),
                      var = rep(var_g, max_rep),
                      varGxE = rep(var_ge, max_rep))
        # Collect pedigree
        SP$setTrackPed(TRUE)
        # Create founder parents
        parents <- AlphaSimR::newPop(founder_pop, simParam = SP)

        ## filling the breeding pipeline
        # TODO: could implement a variation of the scheme here
        for(icohort in 1:info$n_cohort) {
          cli::cli_progress_message("Filling pipeline stage: {icohort} out of {info$n_cohort}")
          if(icohort < (info$n_cohort + 1)) {
            F1 <- AlphaSimR::randCross(parents, info$n_crosses, simParam = SP)
          }
          if(icohort < info$n_cohort){
            F1 <- AlphaSimR::setPheno(F1, varE = rep(var_e, max_rep), simParam = SP)
            selected <- AlphaSimR::selectInd(F1, info$n_parents, simParam = SP)
          }
        }
        cli::cli_progress_done()
        ## running the burn-in phase
        for(year in 1:info$n_burnin) {
          cli::cli_progress_message("Working on burn-in year: {year}")
          parents <- selected
          F1 <- AlphaSimR::setPheno(F1, varE = rep(var_e, max_rep), simParam = SP)
          selected <- AlphaSimR::selectInd(F1, info$n_parents, simParam = SP)
          F1 <- AlphaSimR::randCross(parents, info$n_crosses, simParam = SP)
        }
        cli::cli_progress_done()

        ## running the future phase with competing scenarios
        for(year in (info$n_burnin + 1):(info$n_burnin + info$n_future)) {
          cli::cli_progress_message("Working on future year: {year}")
          parents <- selected
          F1 <- AlphaSimR::setPheno(F1, varE = rep(var_e, max_rep), simParam = SP)
          selected <- AlphaSimR::selectInd(F1, info$n_parents, simParam = SP)
          F1 <- AlphaSimR::randCross(parents, info$n_crosses, simParam = SP)
        }
        cli::cli_progress_done()
        # a bit opinionated way of making sure it fits the exp design
        # maybe need to rethink if this should be the default
        final <- AlphaSimR::randCross(F1, info$n_geno, simParam = SP)
        pheno <- AlphaSimR::setPheno(final, varE = rep(var_e, max_rep), simParam = SP, onlyPheno = TRUE)
        geno <- as.integer(get(info$geno))
        out[[rname]] <- pheno[cbind(geno,
                                    sapply(1:length(geno), function(x) sum(geno[1:x]==geno[x])))]

        # ensure it is in a valid range
        valid_lower <- switch(expect[[rname]]$operator,
                              "greaterThan" = expect[[rname]]$value,
                              "greaterThanOrEqual" = expect[[rname]]$value,
                              "equal" = expect[[rname]]$value,
                              "between" = expect[[rname]]$value[1],
                              "lessThanOrEqual" = NA,
                              "lessThan" = NA)
        valid_upper <- switch(expect[[rname]]$operator,
                              "greaterThan" = NA,
                              "greaterThanOrEqual" = NA,
                              "equal" = expect[[rname]]$value,
                              "between" = expect[[rname]]$value[2],
                              "lessThanOrEqual" = expect[[rname]]$value,
                              "lessThan" = expect[[rname]]$value)
        if(!is.na(valid_upper)) {
          out[[rname]] <- rescale_values(out[[rname]],
                                         lower = valid_lower + 5 * .Machine$double.eps,
                                         upper = valid_upper - 5 * .Machine$double.eps)
        }
      }
      as.data.frame(out)
    })

  simulate_rcrds(design_sim,
                 .joint = with_params(rcrds_names = rcrds_names,
                                      expect = expect_new,
                                      info = info,
                                      chat = chat,
                                      initial = FALSE),
                 .seed = seed,
                 .nsim = nsim)
}
