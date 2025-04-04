get_nsga_fitness = function(
    seed = NULL,
    fitness.fun, 
    n.objectives = NULL, 
    n.dim = NULL, 
    minimize = NULL, 
    lower = NULL, 
    upper = NULL, 
    mu = 100L, 
    lambda = 1,
    mutator = setup(mutPolynomial, eta = 25, p = 0.2, lower = lower, upper = upper), 
    recombinator = setup(recSBX, eta = 15, p = 0.7, lower = lower, upper = upper), 
    terminators = list(stopOnIters(100L)), 
    log.pop = TRUE,
    ...
){
  
  
  optimizer = nsga2(
    seed = seed,
    fitness.fun = fitness.fun, 
    n.objectives = n.objectives, 
    n.dim = n.dim, 
    minimize = minimize, 
    lower = lower, 
    upper = upper, 
    mu = mu, 
    lambda = lambda,
    mutator = mutator, 
    recombinator = recombinator, 
    terminators = terminators, 
    log.pop = log.pop
  )
  
  optimizer$log$env$log.pop
  pop <- optimizer$log$env$pop
  
  clean_pop <- Filter(Negate(is.null), pop)
  
  all_sol <- clean_pop[[1]]$fitness
  
  
  for (i in 2:length(clean_pop)) {
    all_sol <- cbind(all_sol, clean_pop[[i]]$fitness)
  }
  
  
  return(all_sol)
  
}


nsga2 = function(
    seed = NULL,
    fitness.fun,
    n.objectives = NULL,
    n.dim = NULL,
    minimize = NULL,
    lower = NULL,
    upper = NULL,
    mu = 100L,
    lambda = mu,
    mutator = setup(mutPolynomial, eta = 25, p = 0.2, lower = lower, upper = upper),
    recombinator = setup(recSBX, eta = 15, p = 0.7, lower = lower, upper = upper),
    terminators = list(stopOnIters(100L)),
    ...) {
  
  res = ecr_nsga2(seed = seed, fitness.fun = fitness.fun, n.objectives = n.objectives,
                  n.dim = n.dim, minimize = minimize, lower = lower, upper = upper,
                  mu = mu, lambda = lambda, representation = "float", survival.strategy = "plus",
                  parent.selector = selSimple,
                  mutator = mutator,
                  recombinator = recombinator,
                  survival.selector = selNondom,
                  terminators = terminators, ...)
  return(res)
}

ecr_nsga2 = function(
    seed = NULL, fitness.fun, minimize = NULL, n.objectives = NULL,
    n.dim = NULL, lower = NULL, upper = NULL, n.bits,
    representation, mu, lambda, perm = NULL,
    p.recomb = 0.7, p.mut = 0.3,
    survival.strategy = "plus", n.elite = 0L,
    log.stats = list(fitness = list("min", "mean", "max")),
    log.pop = FALSE,
    monitor = NULL,
    initial.solutions = NULL,
    parent.selector = NULL,
    survival.selector = NULL,
    mutator = NULL,
    recombinator = NULL,
    terminators = list(stopOnIters(100L)),
    ...) {
  
  
  
  if (!is.null(seed)){
    
    set.seed(seed)
  }
  
  
  if (!isSmoofFunction(fitness.fun)) {
    n.objectives = asInt(n.objectives, lower = 1L)
    if (is.null(minimize))
      minimize = rep(TRUE, n.objectives)
  }
  
  if (isSmoofFunction(fitness.fun)) {
    n.objectives = getNumberOfObjectives(fitness.fun)
    n.dim = getNumberOfParameters(fitness.fun)
    par.set = getParamSet(fitness.fun)
    upper = getUpper(par.set)
    lower = getLower(par.set)
  }
  
  assertChoice(representation, c("binary", "float", "permutation", "custom"))
  assertChoice(survival.strategy, c("comma", "plus"))
  assertNumber(p.recomb, lower = 0, upper = 1)
  assertNumber(p.mut, lower = 0, upper = 1)
  assertFlag(log.pop)
  assertList(terminators, any.missing = FALSE, all.missing = FALSE, types = "ecr_terminator")
  mu = asInt(mu, lower = 1L)
  lambda.lower = if (survival.strategy == "plus") 1L else mu
  lambda = asInt(lambda, lower = lambda.lower)
  
  control = initECRControl(fitness.fun, n.objectives = n.objectives, minimize = minimize)
  #FIXME: ugly! get rid of the following line
  control$type = representation
  
  n.objectives = control$task$n.objectives
  
  if (representation != "custom" | !is.null(mutator))
    
    control = registerECROperator(control, "mutate", coalesce(mutator, getDefaultEvolutionaryOperators(representation, "mutator", n.objectives, control)))
  if (representation != "custom" | !is.null(recombinator))
    control = registerECROperator(control, "recombine", coalesce(recombinator, getDefaultEvolutionaryOperators(representation, "recombinator", n.objectives, control)))
  control = registerECROperator(control, "selectForSurvival", coalesce(survival.selector, getDefaultEvolutionaryOperators(representation, "survival.selector", n.objectives, control)))
  control = registerECROperator(control, "selectForMating", coalesce(parent.selector, getDefaultEvolutionaryOperators(representation, "parent.selector", n.objectives, control)))
  
  # init logger
  log = initLogger(control,
                   log.stats = log.stats,
                   log.pop = log.pop, init.size = 20000L)
  
  # generate population (depends on representation)
  gen.fun = NULL
  gen.pars = list()
  if (representation == "binary") {
    gen.fun = genBin; gen.pars = list(n.dim = n.bits)
  } else if (representation == "float") {
    gen.fun = genReal; gen.pars = list(n.dim = n.dim, lower = lower, upper = upper)
  } else if (representation == "permutation") {
    gen.fun = genPerm; gen.pars = list(n.dim = perm)
  } else {
    if (!is.null(initial.solutions)) {
      if (length(initial.solutions) != mu) {
        stopf("For custom representations the number of initial solutions need to be equal to mu.")
      }
    } else {
      stopf("For custom representations intial solutions need to be passed.")
    }
  }
  
  population = initial.solutions
  if (representation != "custom")
    population = do.call(initPopulation, c(list(mu = mu, gen.fun = gen.fun, initial.solutions = initial.solutions), gen.pars))
  fitness = evaluateFitness(control, population, ...)
  
  for (i in seq_along(population)) {
    attr(population[[i]], "fitness") = fitness[, i]
  }
  
  updateLogger(log, population, fitness = fitness, n.evals = mu)
  
  
  repeat {
    # generate offspring
    offspring = generateOffspring(control, population, fitness, lambda = lambda, p.recomb = p.recomb, p.mut = p.mut)
    fitness.offspring = evaluateFitness(control, offspring, ...)
    for (i in seq_along(offspring)) {
      attr(offspring[[i]], "fitness") = fitness.offspring[, i]
    }
    
    sel = if (survival.strategy == "plus") {
      replaceMuPlusLambda(control, population, offspring, fitness, fitness.offspring)
    } else {
      replaceMuCommaLambda(control, population, offspring, fitness, fitness.offspring, n.elite = n.elite)
    }
    
    population = sel$population
    fitness = sel$fitness
    
    # do some logging
    updateLogger(log, offspring, fitness.offspring, n.evals = lambda)
    
    
    stop.object = doTerminate(log, terminators)
    if (length(stop.object) > 0L)
      break
  }
  return(makeECRResult(control, log, population, fitness, stop.object))
}


doTerminate = function(log, stop.conds) {
  stop.object = list()
  # if we have not specified any stopping conditions always return the empty object
  if (!length(stop.conds)) {
    return(stop.object)
  }
  
  # otherwise iterate over stopping conditions and check
  for (stop.conds in stop.conds) {
    shouldStop = stop.conds(log = log)
    if (shouldStop) {
      stop.object$name = attr(stop.conds, "name")
      stop.object$message = attr(stop.conds, "message")
      break
    }
  }
  return(stop.object)
}


makeECRResult = function(control, log, population, fitness, stop.object, ...) {
  n.objectives = control$task$n.objectives
  if (n.objectives == 1L)
    return(setupResult.ecr_single_objective(population, fitness, control, log, stop.object, ...))
  moo.res = setupResult.ecr_multi_objective(population, fitness, control, log, stop.object, ...)
  moo.res = filterDuplicated(moo.res)
  return(moo.res)
}

transformFitness = function(fitness, task, selector) {
  # logical vector of opt directions
  task.dir = task$minimize
  # "vectorize" character indicating supported opt direction by selector
  sup.dir = rep(attr(selector, "supported.opt.direction"), task$n.objectives)
  # "logicalize" selector opt direction
  sup.dir = (sup.dir == "minimize")
  
  fn.scale = ifelse(xor(task.dir, sup.dir), -1, 1)
  
  # build transformation matrix
  fn.scale = if (task$n.objectives == 1L) {
    #FIXME: R BUG?!?!
    # diag(ifelse(xor(task.dir, sup.dir), -1, 1)) breaks with message
    # Fehler in diag(ifelse(xor(task.dir, sup.dir), -1, 1)) : ung"ultiger 'nrow' Wert (< 0)
    # if n.objectives is 1! -.-
    # Weird R bug??? diag(1) works!
    as.matrix(fn.scale)
  } else {
    diag(fn.scale)
  }
  
  # transform fitness
  return(fn.scale %*% fitness)
}


setupResult.ecr_multi_objective = function(population, fitness, control, log, stop.object) {
  fitness = transformFitness(fitness, control$task, control$selectForMating)
  
  makeS3Obj(
    task = control$task,
    log = log,
    last.population = population,
    message = stop.object$message,
    classes = c("ecr_multi_objective_result", "ecr_result")
  )
}

selNondom = makeSelector(
  selector = function(fitness, n.select) {
    nondom.layers = doNondominatedSorting(fitness)
    
    # storage for indizes of selected individuals
    new.pop.idxs = integer()
    
    # get maximum rank, i.e., the number of domination layers
    max.rank = max(nondom.layers$ranks)
    
    # get the indizes of points for each domination layer
    idxs.by.rank = lapply(seq(max.rank), function(r) which(nondom.layers$ranks == r))
    
    # get the number of points in each domination layer ...
    front.len = sapply(idxs.by.rank, length)
    
    # ... cumulate the number of points of the domination layers ...
    cum.front.len = cumsum(front.len)
    
    # ... and determine the first domination layer, which does not fit as a whole
    front.first.nonfit = which.first(cum.front.len > n.select)
    
    if (front.first.nonfit > 1L) {
      # in this case at least one nondominated front can be added
      new.pop.idxs = unlist(idxs.by.rank[1:(front.first.nonfit - 1L)])
    }
    
    # how many points to select by second criterion, i.e., crowding distance?
    n.diff = n.select - length(new.pop.idxs)
    
    if (n.diff > 0L) {
      idxs.first.nonfit = idxs.by.rank[[front.first.nonfit]]
      cds = computeCrowdingDistance(fitness[, idxs.first.nonfit, drop = FALSE])
      idxs2 = order(cds, decreasing = TRUE)[1:n.diff]
      new.pop.idxs = c(new.pop.idxs, idxs.first.nonfit[idxs2])
    }
    
    # merge the stuff and return
    return(new.pop.idxs)
  },
  supported.objectives = "multi-objective")


setup = function(operator, ...) {
  assertClass(operator, "ecr_operator")
  args = list(...)
  attrs = attributes(operator)
  fn = function(x, ...) {
    do.call(operator, c(list(x), BBmisc::insert(args, list(...))))
  }
  attributes(fn) = attrs
  return(fn)
}

getDefaultEvolutionaryOperators = function(representation, type, n.objectives, control) {
  if (n.objectives == 1L) {
    return(getSingleObjectiveDefaults(representation, type, control))
  }
  return(getMultiObjectiveDefaults(representation, type, control))
}

getMultiObjectiveDefaults = function(representation, type, control) {
  defaults = list(
    "float" = list(
      "parent.selector" = setup(selSimple),
      "mutator" = try(setup(mutGauss), silent = TRUE),
      "recombinator" = setup(recIntermediate),
      "survival.selector" = setup(selNondom)
    ),
    "binary" = list(
      "parent.selector" = setup(selSimple),
      "mutator" = setup(mutBitflip),
      "recombinator" = setup(recCrossover),
      "survival.selector" = setup(selNondom)
    ),
    "permutation" = list(
      "parent.selector" = setup(selSimple),
      "mutator" = setup(mutSwap),
      "recombinator" = setup(recPMX),
      "survival.selector" = setup(selNondom)
    ),
    "custom" = list(
      "parent.selector" = setup(selSimple),
      "mutator" = NULL,
      "recombinator" = NULL,
      "survival.selector" = setup(selNondom)
    )
  )
  
  if (representation %in% names(defaults)) {
    return(defaults[[representation]][[type]])
  }
  stopf("No defaults availiable for custom representation. You need to specify all
    operators by hand.")
}