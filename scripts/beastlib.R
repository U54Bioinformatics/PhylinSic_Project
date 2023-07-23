# Functions:
# run.beast            Do a phylogeny analysis and return the model.
# continue.beast       Resume an analysis, official implementation.
# run.more.beast       Resume an analysis, my implementation.
# read.beast.model     Read in a previously analyzed beast model.
# make.consensus.tips  Return the order of the tips in a consensus tree.
#
# write.summary.stats  Write out a table with the summary statistics.
# write.ess.table      Write out a table with the ESS.
#
# plot.posterior       Make a plot of posterior over iterations.
# plot.tree.height     Make a plot of posterior over iterations.
# plot.bd.birthrate    Make a plot of posterior over iterations.
#
# plot.densitree
#
# name.lineage
# write.tree.metadata
# write.mcc.tree.metadata
#
#
# Smoothing functions:
# make.genotype.and.probability.matrices
# calc.p.dist.matrix
# calc.p.dist
#
# calc.cell.similarity
# assign.genotype.from.p
# sample.genotype.from.p
#
# make.dropout.matrix



run.beast <- function(
  fasta.file, SITE.MODEL="jc69", CLOCK.MODEL="strict", TREE.PRIOR="bd",
  cbs_group_sizes_dimension=5,
  mrca.is_monophyletic=NULL,
  mrca.normal.prior=NULL,
  tipdates.file=NULL,
  num.iter=1000000, sample.interval=1000,
  beast.infile=NULL, beast.outfile=NULL, tree.log=NULL,
  trace.log=NULL, screen.log=NULL, beast2.path=NULL, rng.seed=0) {
  # SITE.MODEL
  # jc69    Simplest.  equal mutation rates.
  # hky     transition, transversion.  base freq
  # tn93    different transitions
  # gtr     most general
  #
  # CLOCK.MODEL
  # strict  Same rate of mutation for every branch
  # rln     Relaxed.
  #
  # TREE.PRIOR
  # bd      birth death
  # ccp     Coalescent Constant Population
  # cep     Coalescent Exponential Population
  # cbd     Coalescent Bayesian Skyline
  # yule    Constant birth rate.  Better for speciation events.
  #
  # cbs_group_sizes_dimension is the group_sizes_dimension paramter
  # for the cbs tree prior.
  # mrca.is_monophyletic should be TRUE or FALSE.
  # mrca.normal.prior should be c(<mean height>, <sd height>).
  # tipdates.file is the name of a file where the first column is the
  # cell, and the second is the date that the cell was profiled
  # (e.g. days).
  library(beautier)
  library(babette)

  unlink.trace <- FALSE
  unlink.tree <- FALSE
  unlink.screen <- FALSE
  unlink.beast.infile <- FALSE
  unlink.beast.outfile <- FALSE
  if(is.null(trace.log)) {
    trace.log <- tempfile(tmpdir=".")
    unlink.trace <- TRUE
  }
  if(is.null(tree.log)) {
    tree.log <- tempfile(tmpdir=".")
    unlink.tree <- TRUE
  }
  if(is.null(screen.log)) {
    screen.log <- tempfile(tmpdir=".")
    unlink.screen <- TRUE
  }
  if(is.null(beast.infile)) {
    beast.infile <- tempfile(tmpdir=".")
    unlink.beast.infile <- TRUE
  }
  if(is.null(beast.outfile)) {
    beast.outfile <- tempfile(tmpdir=".")
    unlink.beast.outfile <- TRUE
  }

  if(SITE.MODEL == "jc69") {
    # Simplest.  equal mutation rates.
    site.model <- create_jc69_site_model()
  } else if(SITE.MODEL == "hky") {
    # transition, transversion.  base freq
    site.model <- create_hky_site_model()
  } else if(SITE.MODEL == "tn93") {
    # different transitions
    site.model <- create_tn93_site_model()
  } else if(SITE.MODEL == "gtr") {
    # most general
    site.model <- create_gtr_site_model()
  } else {
    stop(sprintf("Unknown site model: %s", SITE.MODEL))
  }

  if(CLOCK.MODEL == "strict") {
    clock.model <- create_strict_clock_model()
  } else if(CLOCK.MODEL == "rln") {
    clock.model <- create_rln_clock_model()
  } else {
    stop(sprintf("Unknown clock model: %s", CLOCK.MODEL))
  }

  if(TREE.PRIOR == "bd") {
    # birth death
    tree.prior <- create_bd_tree_prior()
  } else if(TREE.PRIOR == "cbs") {
    # Coalescent Bayesian Skyline
    tree.prior <- create_cbs_tree_prior(
      group_sizes_dimension=cbs_group_sizes_dimension)
  } else if(TREE.PRIOR == "ccp") {
    # Coalescent Constant Population
    tree.prior <- create_ccp_tree_prior()
  } else if(TREE.PRIOR == "cep") {
    # Coalescent Exponential Population
    tree.prior <- create_cep_tree_prior()
  } else if(TREE.PRIOR == "yule") {
    # Yule
    tree.prior <- create_yule_tree_prior()
  } else {
    stop(sprintf("Unknown tree prior: %s", TREE.PRIOR))
  }

  params <- list()
  if(!is.null(mrca.is_monophyletic))
    params[["is_monophyletic"]] <- mrca.is_monophyletic
  if(!is.null(mrca.normal.prior)) {
    if(length(mrca.normal.prior) != 2)
      stop("mrca.normal.prior should be c(mean, sigma).")
    m <- mrca.normal.prior[1]
    s <- mrca.normal.prior[2]
    params[["mrca_distr"]] <- create_normal_distr(mean=m, sigma=s)
  }
  mrca.prior <- do.call(create_mrca_prior, params)
  # mrca_prior= create_mrca_prior(
  #   taxa_names=subset_nms[cancer.immune.normal == "cancer"]$Cell ,
  #   is_monophyletic=T,
  #   mrca_distr=create_normal_distr(mean=-treeheightPriormu, sigma=treeheightPriorSigma)  )

  x1 <- create_tracelog(trace.log, log_every=sample.interval)
  x2 <- create_treelog(tree.log, log_every=sample.interval)
  x3 <- create_screenlog(screen.log, log_every=sample.interval)

  mcmc <- create_mcmc(
    chain_length=num.iter,
    store_every=NA,
    tracelog=x1,
    treelog=x2,
    screenlog=x3
    )

  params <- list()
  params[["site_model"]] <- site.model
  params[["clock_model"]] <- clock.model
  params[["tree_prior"]] <- tree.prior
  params[["mrca_prior"]] <- mrca.prior
  params[["mcmc"]] <- mcmc
  if(!is.null(tipdates.file))
    params[["tipdates_filename"]] <- tipdates.file
  inference_model <- do.call(create_inference_model, params)

  # specify input and output file names
  if(is.null(beast2.path))
    beast2.path <- get_default_beast2_path()
  beast2_options <- create_beast2_options(
    input_filename=beast.infile,
    output_state_filename=beast.outfile,
    beast2_path=beast2.path,
    rng_seed=rng.seed)
  #beast2_options$verbose <- TRUE

  ### Fit tree
  # Will create beast2_options$input_filename
  beast.model <- bbt_run_from_model(
    fasta_filename=fasta.file,
    inference_model=inference_model,
    beast2_options=beast2_options
  )

  # Clean up the model.
  x <- names(beast.model)
  if(length(x) != 4) stop("unexpected model format")
  x <- x[!(x %in% c("estimates", "operators", "output"))]
  if(length(x) != 1) stop("bad")
  beast.model[["trees"]] <- beast.model[[x]]
  beast.model[[x]] <- NULL

  if(unlink.trace)
    unlink(trace.log)
  if(unlink.screen)
    unlink(screen.log)
  if(unlink.tree)
    unlink(tree.log)
  if(unlink.beast.infile)
    unlink(beast.infile)
  if(unlink.beast.outfile)
    unlink(beast.outfile)

  return(beast.model)
}


continue.beast <- function(
  fasta.file, SITE.MODEL="jc69", CLOCK.MODEL="strict", TREE.PRIOR="bd",
  cbs_group_sizes_dimension=5,
  mrca.is_monophyletic=NULL,
  mrca.normal.prior=NULL,
  tipdates.file=NULL,
  num.iter=1000000, sample.interval=1000,
  beast.infile=NULL, beast.outfile=NULL,
  tree.log=NULL, trace.log=NULL, screen.log=NULL, rng.seed=0) {
  # beast.infile and beast.outfile should be the files from the first
  # analysis.
  unlink.trace <- FALSE
  unlink.tree <- FALSE
  unlink.screen <- FALSE
  unlink.beast.infile <- FALSE
  unlink.beast.outfile <- FALSE
  if(is.null(trace.log)) {
    trace.log <- tempfile(tmpdir=".")
    unlink.trace <- TRUE
  }
  if(is.null(tree.log)) {
    tree.log <- tempfile(tmpdir=".")
    unlink.tree <- TRUE
  }
  if(is.null(screen.log)) {
    screen.log <- tempfile(tmpdir=".")
    unlink.screen <- TRUE
  }
  if(is.null(beast.infile)) {
    beast.infile <- tempfile(tmpdir=".")
    unlink.beast.infile <- TRUE
  }
  if(is.null(beast.outfile)) {
    beast.outfile <- tempfile(tmpdir=".")
    unlink.beast.outfile <- TRUE
  }

  if(SITE.MODEL == "jc69") {
    # Simplest.  equal mutation rates.
    site.model <- create_jc69_site_model()
  } else if(SITE.MODEL == "hky") {
    # transition, transversion.  base freq
    site.model <- create_hky_site_model()
  } else if(SITE.MODEL == "tn93") {
    # different transitions
    site.model <- create_tn93_site_model()
  } else if(SITE.MODEL == "gtr") {
    # most general
    site.model <- create_gtr_site_model()
  } else {
    stop(sprintf("Unknown site model: %s", SITE.MODEL))
  }

  if(CLOCK.MODEL == "strict") {
    clock.model <- create_strict_clock_model()
  } else if(CLOCK.MODEL == "rln") {
    clock.model <- create_rln_clock_model()
  } else {
    stop(sprintf("Unknown clock model: %s", CLOCK.MODEL))
  }

  if(TREE.PRIOR == "bd") {
    # birth death
    tree.prior <- create_bd_tree_prior()
  } else if(TREE.PRIOR == "cbs") {
    # Coalescent Bayesian Skyline
    tree.prior <- create_cbs_tree_prior(
      group_sizes_dimension=cbs_group_sizes_dimension)
  } else if(TREE.PRIOR == "ccp") {
    # Coalescent Constant Population
    tree.prior <- create_ccp_tree_prior()
  } else if(TREE.PRIOR == "cep") {
    # Coalescent Exponential Population
    tree.prior <- create_cep_tree_prior()
  } else if(TREE.PRIOR == "yule") {
    # Yule
    tree.prior <- create_yule_tree_prior()
  } else {
    stop(sprintf("Unknown tree prior: %s", TREE.PRIOR))
  }

  # TODO: Should merge the creation of the priors with the code in
  # run.beast.
  params <- list()
  if(!is.null(mrca.is_monophyletic))
    params[["is_monophyletic"]] <- mrca.is_monophyletic
  if(!is.null(mrca.normal.prior)) {
    if(length(mrca.normal.prior) != 2)
      stop("mrca.normal.prior should be c(mean, sigma).")
    m <- mrca.normal.prior[1]
    s <- mrca.normal.prior[2]
    params[["mrca_distr"]] <- create_normal_distr(mean=m, sigma=s)
  }
  mrca.prior <- do.call(create_mrca_prior, params)

  x1 <- create_tracelog(trace.log, log_every=sample.interval)
  x2 <- create_treelog(tree.log, log_every=sample.interval)
  x3 <- create_screenlog(screen.log, log_every=sample.interval)

  mcmc <- create_mcmc(
    chain_length=num.iter,
    store_every=NA,
    tracelog=x1,
    treelog=x2,
    screenlog=x3
    )
  params <- list()
  params[["site_model"]] <- site.model
  params[["clock_model"]] <- clock.model
  params[["tree_prior"]] <- tree.prior
  params[["mrca_prior"]] <- mrca.prior
  params[["mcmc"]] <- mcmc
  if(!is.null(tipdates.file))
    params[["tipdates_filename"]] <- tipdates.file
  inference_model <- do.call(create_inference_model, params)

  # specify input and output file names
  beast2_options <- create_beast2_options(
    input_filename=beast.infile,
    output_state_filename=beast.outfile,
    rng_seed=rng.seed)
  create_beast2_input_file_from_model(
    input_filename=fasta.file,
    output_filename=beast2_options$input_filename,
    inference_model=inference_model)

  # Fit tree
  beast.model <- bbt_continue(
    fasta_filename=fasta.file,
    inference_model=inference_model,
    beast2_options=beast2_options
  )

  # Clean up the model.
  x <- names(beast.model)
  if(length(x) != 4) stop("unexpected model format")
  x <- x[!(x %in% c("estimates", "operators", "output"))]
  if(length(x) != 1) stop("bad")
  beast.model[["trees"]] <- beast.model[[x]]
  beast.model[[x]] <- NULL

  if(unlink.trace)
    unlink(trace.log)
  if(unlink.screen)
    unlink(screen.log)
  if(unlink.tree)
    unlink(tree.log)
  if(unlink.beast.infile)
    unlink(beast.infile)
  if(unlink.beast.outfile)
    unlink(beast.outfile)

  return(beast.model)
}


run.more.beast <- function(beast.infile, beast.outfile, num.iter=NULL,
  rng.seed=0, beast2.folder=NULL) {
  # If num.iter is NULL, will run with the same iterations as before.
  # Otherwise, will run for this number of iterations.
  # beast2.folder is where beast2 is installed.  Default is
  # "/usr/local".

  if(!file.exists(beast.infile))
    stop(sprintf("Cannot find infile: %s", beast.infile))
  if(!file.exists(beast.outfile))
    stop(sprintf("Cannot find infile: %s", beast.outfile))

  # Read the trace.log, screen.log, and tree.log files.
  # <logger id="tracelog" fileName="./file3db6a7074c3" logEvery="1000" model="@posterior" sanitiseHeaders="true" sort="smart">
  # <logger id="screenlog" fileName="./file3db440203b0" logEvery="1000">
  # <logger id="treelog.t:temp.987.1.125" fileName="test.tree.log" logEvery="1000" mode="tree">

  tree.log <- .get.log.filename(beast.infile, "treelog")
  trace.log <- .get.log.filename(beast.infile, "tracelog")
  screen.log <- .get.log.filename(beast.infile, "screenlog")

  # If tree.log doesn't exist, BEAST2 will give a warning and create a
  # new one.
  #
  # WARNING: Resuming, but file test.tree.log does not exist yet
  # (perhaps the seed number is not the same as before?).
  # Writing new file test.tree.log
  #
  # But this function is useless if the tree log doesn't exist.  Make
  # sure it exists.
  if(!file.exists(tree.log))
    stop("tree.log doesn't exist.  Specify tree.log when running run.beast.")

  # If screen.log doesn't exist, BEAST2 will create a new one.
  # If trace.log doesn't exist, BEAST2 will create a new one.
  unlink.trace <- !file.exists(trace.log)
  unlink.screen <- !file.exists(screen.log)

  # Update num.iter if needed.
  if(!is.null(num.iter)) {
    if((num.iter < 1) || (num.iter > 1E10)) stop("bad iterations")
    lines <- readLines(beast.infile)

    # Change the chainLength.
    # <run id="mcmc" spec="MCMC" chainLength="50000">
    I <- grep("chainLength=", lines)
    if(is.null(I)) stop(sprintf("Could not find chainLength %s", beast.infile))
    if(length(I) > 1) stop("multiple chainLengths")
    x <- gsub('.+chainLength="(\\d+)".+', "\\1", lines[I], perl=TRUE)
    if(x != num.iter) {
      x <- sprintf('chainLength="%s"', format(num.iter, scientific=FALSE))
      x <- gsub('chainLength="\\d+"', x, lines[I], perl=TRUE)
      if(x == lines[I]) stop("could not parse")
      lines[I] <- x
      writeLines(lines, beast.infile)
    }
  }

  if(is.null(beast2.folder))
    beast2.folder <- get_default_beast2_folder()
  launcher.file <- sprintf("%s/beast/lib/launcher.jar", beast2.folder)
  if(!file.exists(launcher.file))
    stop(sprintf("Cannot find launcher file: %s", launcher.file))

  # Run BEAST2.
  # java -cp /usr/local/beast/lib/launcher.jar
  #   beast.app.beastapp.BeastLauncher \
  #  -seed 1 -overwrite -resume -statefile /data/jchang4/biocore5/outfile.txt
  #  /data/jchang4/biocore5/infile.txt
  x1 <- gsub(" ", "\\\\ ", beast.infile)
  x2 <- gsub(" ", "\\\\ ", beast.outfile)
  cmd <- c(
    "-cp", launcher.file,
    "beast.app.beastapp.BeastLauncher",
    "-seed", rng.seed, "-overwrite", "-resume", "-statefile",
    x2, x1)
  output <- system2("java", cmd, stdout=TRUE, stderr=TRUE)

  out <- read.beast.model(tree.log, trace.log, beast.outfile)
  out$output <- output

  if(unlink.trace)
    unlink(trace.log)
  if(unlink.screen)
    unlink(screen.log)

  out
}


.get.log.filename <- function(beast.infile, logid) {
  # Pull out the line with the filename.
  lines <- readLines(beast.infile)
  x <- sprintf('<logger id="%s', logid)
  I <- grep(x, lines)
  if(is.null(I)) stop(sprintf("Cannot find log %s", logid))
  if(length(I) > 1) stop("multiple logs")
  l <- lines[I]

  # <logger id="screenlog" fileName="./file3db440203b0" logEvery="1000">
  x <- gsub('.+fileName="([^\\"]+)".+', "\\1", l, perl=TRUE)
  if(l == x) stop("could not find fileName")
  x
}


read.beast.model <- function(tree.file, trace.file, beast.outfile) {
  out <- tracerer::parse_beast_output_files(
    trees_filenames=tree.file,
    log_filename=trace.file,
    state_filename=beast.outfile)

  ## Get the alignment ID.
  #lines <- readLines(beast.infile)
  ## id="temp.30707.1.125"
  #I <- grep('^id="([\\w.]+?)"$', lines, perl=TRUE)
  #if(is.null(I)) stop(sprintf("Cannot find alignment id %s", beast.infile))
  #if(length(I) > 1) stop("multiple alignment IDs")
  #x <- gsub('id="([\\w.]+?)"', "\\1", lines[I], perl=TRUE)
  #if(x == lines[I]) stop("could not parse")
  #alignment.id <- x
  #names(out)[1] <- paste0(alignment.id, "_trees")
  names(out)[1] <- "trees"
  #out$output <- output

  out
}


make.consensus.tips <- function(beast.model, burnin.frac=0.1, max.trees=200) {
  # Return a vector of the tip labels.
  x <- beast.model$trees
  trees <- x[floor(burnin.frac*length(x)):length(x)]
  if(length(trees) > max.trees) {
    I <- seq(1, length(trees), length(trees)/max.trees)
    I <- round(I)
    I <- I[!duplicated(I)]
    trees <- trees[I]
  }


  if(0) {
    st <- superTree(trees, method="MRP", rooted=TRUE)
  } else {
    st <- consensus(trees, p=1)
  }
  x <- st$edge[,2]
  x <- x[x <= length(st$tip.label)]
  tips <- st$tip.label[x]
  return(tips)
}


write.ess.table <- function(filename, beast.model, sample.interval=1000,
  burnin.frac=0.1) {
  traces <- remove_burn_ins(traces=beast.model$estimates,
    burn_in_fraction=burnin.frac)

  esses <- t(calc_esses(traces, sample_interval=sample.interval))
  colnames(esses) <- "ESS"
  data.out <- cbind(rownames(esses), esses)
  colnames(data.out)[1] <- "Name"
  my.write(data.out, filename, col.names=TRUE)
  return(esses)
}


write.summary.stats <- function(filename, beast.model, sample.interval=1000,
  burnin.frac=0.1) {
  traces <- remove_burn_ins(traces=beast.model$estimates,
    burn_in_fraction=burnin.frac)
  stats <- t(
    calc_summary_stats(traces$posterior, sample_interval=sample.interval))
  colnames(stats) <- "Statistic"
  data.out <- cbind(rownames(stats), stats)
  colnames(data.out)[1] <- "Name"
  my.write(data.out, filename, col.names=TRUE)
  return(stats)
}


plot.posterior <- function(filename, beast.model, burnin=NULL) {
  .plot.estimate(filename, beast.model, "posterior", "red", burnin)
}

plot.tree.height <- function(filename, beast.model, burnin=NULL) {
  .plot.estimate(filename, beast.model, "TreeHeight", "green", burnin)
}

plot.bd.birthrate <- function(filename, beast.model, burnin=NULL) {
  .plot.estimate(filename, beast.model, "BDBirthRate", "blue", burnin)
}


.plot.estimate <- function(filename, beast.model, name, color, burnin) {
  p <- ggplot(data=beast.model$estimates, aes(x=Sample)) +
    theme_classic() +
    xlab("MCMC iteration")
  if(!is.null(burnin))
    p <- p + geom_vline(xintercept=burnin, linetype="dotted")
  pdf(filename, height=8, width=11)
  p <- p + geom_line(aes_string(y=name), color=color)
  print(p)
  dev.off()
}


plot.densitree <- function(filename, beast.model, height=8, width=11,
  alpha=NULL, burnin.frac=0.1, consensus=NULL, color.bar=NULL,
  sequence.mat=NULL) {
  # consensus is a vector of the names of cells.  Should match
  # tip.label in the trees in beast.model$trees.
  # color.bar should be a list of colors, same length as consensus.
  # sequence.mat should be a sequence by base matrix
  # containing "a", "c", "g", "t", or "-".  The row names should
  # be the names of the cells.
  show.color.bar <- FALSE
  if(!is.null(color.bar)) {
    if(is.null(consensus)) stop("color.bar must be in same order as consensus")
    if(length(color.bar) != length(consensus)) stop("length mismatch")
    show.color.bar <- TRUE
  }
  show.matrix <- FALSE
  if(!is.null(sequence.mat)) {
    if(is.null(consensus)) stop("consensus should be provided")
    I <- match(consensus, rownames(sequence.mat))
    if(any(is.na(I))) stop("rownames(sequence.mat) should match consensus")
    nuc.mat <- sequence.mat[I,]
    nuc.mat[is.na(nuc.mat)] <- 0
    nuc.mat[nuc.mat %in% c("-", "?")] <- 0
    nuc.mat[nuc.mat=="c"] <- 1
    nuc.mat[nuc.mat=="g"] <- 2
    nuc.mat[nuc.mat=="a"] <- 3
    nuc.mat[nuc.mat=="t"] <- 4
    class(nuc.mat) <- "numeric"
    nuc.col <- c("white", BREWER_S1_ORANGE, BREWER_S1_BLUE, BREWER_S1_GREEN,
      BREWER_S1_YELLOW)
    show.matrix <- TRUE
  }

  x <- beast.model$trees
  trees <- x[floor(burnin.frac*length(x)):length(x)]

  MAX.TREES <- 100
  if(length(trees) > MAX.TREES) {
    I <- seq(1, length(trees), length(trees)/MAX.TREES)
    I <- round(I)
    I <- I[!duplicated(I)]
    trees <- trees[I]
  }

  if(is.null(alpha)) {
    # Guess a good alpha.
    # TREES  alpha
    #  200    0.01
    #   91    0.03    0.057 is too high.
    #   40    0.06
    #   20     0.1
    #alpha <- 0.1-length(trees)/2300
    alpha <- 1/length(trees)   # default
    alpha <- min(max(alpha, 0.001), 1)
  }
  pdf(filename, height=height, width=width)

  x <- c(1)
  widths <- c(20)
  if(show.color.bar) {
    x <- c(x, max(x)+1)
    widths <- c(widths, 1)
  }
  if(show.matrix) {
    x <- c(x, max(x)+1)
    widths <- c(widths, 8)
  }
  if(length(x) > 1) {
    layout.matrix <- matrix(x, nrow=1)
    layout(mat=layout.matrix, heights=c(1), widths=widths)
  }

  # Make the densitree plot.
  if(show.color.bar | show.matrix) {
    #par(mar=c(5, 4, 4, 0))
    #par(oma=c(0, 0, 0, 0))
    par(mar=c(0.1, 0.1, 0.1, 0))
  }
  plot_densitree(trees, type="cladogram", width=1, scaleX=T, cex=0.01,
    alpha=alpha, consensus=consensus)
  mtext(expression(bold("Bayesian inference")), side=3, line=0.5, cex=1.0)

  # Make the color bar, if needed.
  if(show.color.bar) {
    M <- 2.0
    #par(mar=c(5+M, 0, 4+M, 2))
    #M <- 0.1
    #par(oma=c(0, 0, 0, 0))
    par(mar=c(0.1+M, 0, 0.1+M, 0.1))

    col <- sort(unique(color.bar))
    x <- match(color.bar, col)
    m <- matrix(x, ncol=1)
    image(t(m), axes=FALSE, col=col)
  }

  # Make the matrix, if needed.
  if(show.matrix) {
    M <- 2.0
    par(mar=c(0.1+M, 0, 0.1+M, 0.1))
    # Need to reverse order and transpose for image.
    x <- nuc.mat[seq(nrow(nuc.mat), 1, -1),]
    image(t(x), col=nuc.col, xlab="", ylab="", axes=FALSE)
  }

  # Done
  dev.off()
}


# n               length(tree$tip.label)
# tree$Nnode      Number of internal nodes.  n-1.  (n-2 for unrooted?)
# root node       n+1  (I think)
# tree$root.edge  Missing for rerooted.tree.
#
# - Edges 1-n correspond to trees$tip.label.
# - 1 is on the "left" (root on top, tips on bottom) of the tree.  n
#   is on the top.
# - Lower numbers are on the "left" of the tree.  (I think).
# - Builds internal nodes with depth-first search, left first
#   traversal.

.name.lineage.h <- function(tree, node, parent.name, my.name) {
  lineage.names <- c()
  if(parent.name == "") {
    name <- my.name
  } else {
    name <- sprintf("%s.%s", parent.name, my.name)
  }
  lineage.names <- rbind(lineage.names, c(node, name))

  # Name each of the children.
  I <- which(tree$edge[,1] == node)
  if(length(I) == 0) return(lineage.names)
  J <- sort(tree$edge[I,][,2])
  for(i in 1:length(J)) {
    x <- .name.lineage.h(tree, J[i], name, sprintf("%s", i))
    lineage.names <- rbind(lineage.names, x)
  }
  return(lineage.names)
}

name.lineage <- function(tree) {
  if(!is.rooted(tree)) stop("Cannot name unrooted tree")
  # Find the root of a tree.  Is a node with no parents.
  I <- match(1:max(tree$edge), tree$edge[,2])
  x <- which(is.na(I))
  if(!length(x)) stop("can't find root")
  if(length(x) > 1) stop("found multiple roots")
  root.node <- x[1]
  #library(phytools)
  #heights <- nodeHeights(tree)  # correspond to tree$edge
  #I <- which(heights[,1] == min(heights[,1]))
  #x <- tree$edge[I,1]
  #root.node <- unique(x)
  #if(length(root.node) != 1) stop("can't find root")

  name.tab <- .name.lineage.h(tree, root.node, "", "X")
  lineage.names <- rep(NA, max(tree$edge))
  I <- as.numeric(name.tab[,1])
  lineage.names[I] <- name.tab[,2]
  return(lineage.names)
}

write.tree.metadata <- function(tree, outfile) {
  # Write out the metadata for each of the trees.
  library(phytools)
  num.nodes <- max(tree$edge)  # leaf and internal nodes
  index <- 1:num.nodes
  cell <- rep("", num.nodes)
  cell[1:length(tree$tip.label)] <- tree$tip.label
  lineage.names <- rep("", num.nodes)
  heights <- rep("", num.nodes)
  if(is.rooted(tree)) {
    lineage.names <- name.lineage(tree)
    heights <- rep(NA, num.nodes)
    x <- nodeHeights(tree)  # correspond to tree$edge
    heights[tree$edge[,1]] <- x[,1]
    heights[tree$edge[,2]] <- x[,2]
  }

  metadata <- cbind(index, cell, lineage.names, heights)
  colnames(metadata) <- c("Index", "Cell", "Lineage", "Height")
  my.write(metadata, outfile, col.names=TRUE)
}

write.mcc.tree.metadata <- function(mcc.tree, outfile) {
  # Write out the metadata for each of the trees.
  library(phytools)
  library(tidytree)

  tree <- mcc.tree@phylo  # S3 object
  # tree$edge  matrix of c(<parent>, <child>)
  num.nodes <- max(tree$edge)  # leaf and internal nodes
  index <- 1:num.nodes
  cell <- rep("", num.nodes)
  cell[1:length(tree$tip.label)] <- tree$tip.label
  lineage.names <- rep("", num.nodes)
  heights <- rep("", num.nodes)
  if(is.rooted(tree)) {
    lineage.names <- name.lineage(tree)
    heights <- rep(NA, num.nodes)
    x <- nodeHeights(tree)  # correspond to tree$edge
    heights[tree$edge[,1]] <- x[,1]
    heights[tree$edge[,2]] <- x[,2]
  }

  # Add information from the mcc tree.
  treedata <- as_tibble(mcc.tree)
  # Make sure treedata is aligned with the tree.
  # Make sure the cell names are the same.
  x <- treedata$label == cell
  if(!all(x[!is.na(x)])) stop("unaligned 1")
  # Make sure the parents are the same.
  I <- match(treedata$node, tree$edge[,2])
  parent <- tree$edge[I,1]
  i <- which(is.na(parent))
  parent[i] <- treedata$node[i]
  if(!all(parent == treedata$parent)) stop("unaligned 2")
  # Make sure the heights are the same.
  # cor should be -1.  heights has root at height 0, but treedata has
  # leaves at 0.
  x <- cor(heights, treedata$height)
  if(abs(x+1) > 1E-5) stop("unaligned 3")

  # Assemble metadata file.
  # ignore nodeHeights.  Just write the one from the mcc tree.
  metadata <- cbind(index, cell, lineage.names)
  n <- c("Index", "Cell", "Lineage")

  if(any(names(treedata) == "branch.length")) {
    metadata <- cbind(metadata, treedata$branch.length)
    n <- c(n, "Branch Length")
  }
  if(any(names(treedata) == "height")) {
    x1 <- as.numeric(treedata$height)
    x2 <- unlist(lapply(treedata$height_0.95_HPD, function(x) x[1]))
    x3 <- unlist(lapply(treedata$height_0.95_HPD, function(x) x[2]))
    x4 <- as.numeric(treedata$height_median)
    x5 <- unlist(lapply(treedata$height_range, function(x) x[1]))
    x6 <- unlist(lapply(treedata$height_range, function(x) x[2]))
    metadata <- cbind(metadata, x1, x2, x3, x4, x5, x6)
    n <- c(n, "Height", "Height 0.95 HPD Lo", "Height 0.95 HPD Hi",
      "Height Median", "Height Range Lo", "Height Range Hi")
  }
  if(any(names(treedata) == "length")) {
    x1 <- as.numeric(treedata$length)
    x2 <- unlist(lapply(treedata$length_0.95_HPD, function(x) x[1]))
    x3 <- unlist(lapply(treedata$length_0.95_HPD, function(x) x[2]))
    x4 <- as.numeric(treedata$length_median)
    x5 <- unlist(lapply(treedata$length_range, function(x) x[1]))
    x6 <- unlist(lapply(treedata$length_range, function(x) x[2]))
    metadata <- cbind(metadata, x1, x2, x3, x4, x5, x6)
    n <- c(n, "Length", "Length 0.95 HPD Lo", "Length 0.95 HPD Hi",
      "Length Median", "Length Range Lo", "Length Range Hi")
  }
  if(any(names(treedata) == "posterior")) {
    x <- treedata$posterior
    x[is.na(x)] <- ""
    metadata <- cbind(metadata, x)
    n <- c(n, "Posterior")
  }
  if(any(names(treedata) == "rate")) {
    x1 <- as.numeric(treedata$rate)
    x2 <- unlist(lapply(treedata$rate_0.95_HPD, function(x) x[1]))
    x3 <- unlist(lapply(treedata$rate_0.95_HPD, function(x) x[2]))
    x4 <- as.numeric(treedata$rate_median)
    x5 <- unlist(lapply(treedata$rate_range, function(x) x[1]))
    x6 <- unlist(lapply(treedata$rate_range, function(x) x[2]))
    metadata <- cbind(metadata, x1, x2, x3, x4, x5, x6)
    n <- c(n, "Rate", "Rate 0.95 HPD Lo", "Rate 0.95 HPD Hi",
      "Rate Median", "Rate Range Lo", "Rate Range Hi")
  }
  # The root node has a bunch of NA.
  metadata[is.na(metadata)] <- ""
  colnames(metadata) <- n
  my.write(metadata, outfile, col.names=TRUE)
}




REF.CUT <- 0.30
ALT.CUT <- 0.70


make.genotype.and.probability.matrices <- function(M.ref, M.alt,
  K=NULL, DELTA=NULL, num.features=NULL, num.samples=NULL, sample=FALSE,
  rng.seed=0, num.cores=1) {
  # M.ref and M.alt are matrices of the ref and alt read counts.
  p.dist <- calc.p.dist.matrix(M.ref, M.alt, K=K, DELTA=DELTA,
    num.features=num.features, num.samples=num.samples, rng.seed=rng.seed,
    num.cores=num.cores)
  geno.mat <- matrix(NA, nrow(M.ref), ncol(M.ref))
  p.mat <- matrix(NA, nrow(M.ref), ncol(M.ref))
  for(i in 1:nrow(M.ref)) {
    p <- t(matrix(unlist(p.dist[[i]]), ncol=ncol(M.ref)))
    if(sample) {
      geno <- sample.genotype.from.p(p)
    } else {
      geno <- assign.genotype.from.p(p)
    }
    geno.mat[i,] <- geno
    p.mat[i,] <- apply(p, 1, max)
  }
  list(geno.mat, p.mat)
}

calc.p.dist.matrix <- function(M.ref, M.alt, K=NULL, DELTA=NULL,
  num.features=NULL, num.samples=NULL, rng.seed=0, num.cores=1) {
  # M.ref and M.alt are n x m matrices of the reference and alt
  # counts, where n is the number of sites, and m is the number of
  # cells.
  # Return an n x m x 3 list of lists of probability distributions.
  # K is the number of nearest neighbors to consider, and DELTA is the
  # amount of probability assigned to the neighboring cells.  0 does
  # no smoothing, and 1 does the most smoothing.
  # num.features is the number of features to use to calculate the
  # neighbors.
  # num.samples is the number of genotypes to sample when comparing
  # two cells.
  if(is.null(K))
    K <- 5
  if(is.null(DELTA))
    DELTA <- 0                      # by default, do no smoothing
  if(is.null(num.features))
    num.features <- nrow(M.ref)   # use all features by default
  if(is.null(num.samples))
    num.samples <- 100

   ## Choose the features with the highest VAF.
  #M.vaf <- M.alt / (M.ref + M.alt)
  #M.vaf[is.nan(M.vaf)] <- 0
  #x <- apply(M.vaf, 1, var)
  # Choose the features with the highest coverage.
  x <- apply(M.ref, 1, sum) + apply(M.alt, 1, sum)
  O <- order(x, decreasing=T)
  I <- O[1:num.features]

  # Calculate the probability distributions from cell j.
  # m-length list of n x 3 matrices of probabilities.
  p.j <- mclapply(1:ncol(M.ref), function(j) {
    calc.p.dist(M.ref[,j], M.alt[,j])
    }, mc.cores=num.cores)

  # If DELTA > 0, then smooth the probabilities using k nearest
  # neighbors.
  if(DELTA > 1E-6) {
    # Make a ncol x ncol matrix of cell similarities.  Optimization:
    # calculate the upper diagonal only.
    x <- mclapply(1:ncol(M.ref), function(j1) {
      set.seed(rng.seed)  # make sure this is reproducible.
      S.j <- rep(0, ncol(M.ref))
      for(j2 in j1:ncol(M.ref)) {
        S.j[j2] <- calc.cell.similarity(
          p.j[[j1]][I,], p.j[[j2]][I,], num.samples)
      }
      return(S.j)
    }, mc.cores=num.cores)
    S <- matrix(unlist(x), ncol(M.ref), ncol(M.ref))
    # Make symmetric.
    S <- pmax(S, t(S))

    # Find the nearest neighbors for each cell.
    NN <- matrix(NA, K, ncol(S))
    for(j in 1:ncol(S)) {
      O <- order(S[,j], decreasing=TRUE)
      NN[,j] <- sort(O[2:(2+K-1)])  # ignore highest score (matches itself)
    }

    # For each cell, average the probability distribution of the
    # neighbors.
    p.avgn <- list()   # m-length list of num.features x 3 matrices.
    for(j in 1:ncol(NN)) {
      p.avgn.j <- matrix(0, nrow(M.ref), 3)
      for(j2 in NN[,j])
        p.avgn.j <- p.avgn.j + p.j[[j2]]
      p.avgn[[j]] <- p.avgn.j / nrow(NN)
    }

    # Smooth the probability distribution.
    for(j in 1:length(p.j))
      p.j[[j]] <- (1-DELTA)*p.j[[j]] + DELTA*p.avgn[[j]]
  }

  dist.matrix <- list()
  for(i in 1:nrow(M.ref))
    dist.matrix[[i]] <- list()
  for(j in 1:length(p.j)) {
    p.jj <- p.j[[j]]
    colnames(p.jj) <- NULL
    for(i in 1:nrow(p.jj))
      dist.matrix[[i]][[j]] <- as.numeric(p.jj[i,])
  }
  dist.matrix
}

calc.p.dist <- function(x.ref, x.alt) {
  # Calculate a probability distribution over HOMO-R, HET, HOMO-A.
  # for each cell.  x.ref and x.alt are n-length vectors of read
  # counts for each site.  Return an n x 3 matrix of probabilities.
  # Each row sums to 1.

  x1 <- pbeta(REF.CUT, x.alt+1, x.ref+1)
  x2 <- pbeta(ALT.CUT, x.alt+1, x.ref+1)
  p1 <- x1
  p2 <- x2 - x1
  p3 <- 1 - x2
  probs <- cbind(p1, p2, p3)

  probs
}

calc.cell.similarity <- function(p.geno.1, p.geno.2, num.samples) {
  # Return a score of the similarity between two cells.
  scores <- rep(NA, num.samples)
  for(i in 1:num.samples) {
    geno.1 <- sample.genotype.from.p(p.geno.1)
    geno.2 <- sample.genotype.from.p(p.geno.2)
    scores[i] <- sum(geno.1 == geno.2)
  }
  mean(scores) / nrow(p.geno.1) * 100
}

assign.genotype.from.p <- function(probs) {
  # Given a probability matrix from calc.p.dist, sample the genotypes.
  # Return a vector containing "R", "H", or "A" for homozygous
  # reference, heterozygous, or homozygous alt.
  x <- apply(probs, 1, which.max)
  geno <- rep(NA, nrow(probs))
  geno[x == 1] <- "R"
  geno[x == 2] <- "H"
  geno[x == 3] <- "A"
  geno
}

sample.genotype.from.p <- function(probs) {
  # Given a probability matrix from calc.p.dist, sample the genotypes.
  # Return a vector containing "R", "H", or "A" for homozygous
  # reference, heterozygous, or homozygous alt.
  p <- runif(nrow(probs))
  geno <- rep("A", nrow(probs))
  geno[p<(probs[,1]+probs[,2])] <- "H"
  geno[p<probs[,1]] <- "R"
  geno
}

make.dropout.matrix <- function(M.ref, M.alt) {
  # Return a matrix of 1/0, indicating whether that position is
  # dropped out.
  m <- matrix(0, nrow(M.ref), ncol(M.alt))
  m[(M.ref == 0) & (M.alt == 0)] <- 1
  m
}
