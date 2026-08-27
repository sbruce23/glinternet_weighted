# Demonstrate recovery of a target-distribution interaction under biased sampling.

set.seed(20260825)

interaction_slope <- function(fit) {
  contrast_grid <- rbind(
    c(0, -1), c(1, -1),
    c(0,  1), c(1,  1)
  )
  prediction <- drop(predict(fit, contrast_grid, type = "response"))
  ((prediction[4] - prediction[3]) -
    (prediction[2] - prediction[1])) / 2
}

one_replicate <- function(replicate_id, population_n = 4000L) {
  stratum <- rbinom(population_n, 1, 0.25)
  modifier <- rnorm(population_n)
  treatment <- rbinom(population_n, 1, 0.5)
  outcome <- 0.4 + 0.5 * modifier + 0.3 * treatment +
    (1 + 2 * stratum) * treatment * modifier +
    rnorm(population_n, sd = 0.8)

  inclusion_probability <- ifelse(stratum == 1, 0.75, 0.15)
  included <- rbinom(population_n, 1, inclusion_probability) == 1
  X <- cbind(treatment[included], modifier[included])
  y <- outcome[included]
  inverse_sampling_weight <- 1 / inclusion_probability[included]

  common <- list(
    X = X,
    Y = y,
    numLevels = c(2, 1),
    lambda = 1e-5,
    interactionPairs = matrix(c(1, 2), nrow = 1),
    tol = 1e-8,
    maxIter = 12000
  )
  unweighted <- do.call(glinternet, common)
  weighted <- do.call(
    glinternet,
    c(common, list(weights = inverse_sampling_weight))
  )

  data.frame(
    replicate = replicate_id,
    unweighted = interaction_slope(unweighted),
    weighted = interaction_slope(weighted)
  )
}

target_interaction <- 1 + 2 * 0.25
results <- do.call(rbind, lapply(seq_len(30), one_replicate))
summary_table <- data.frame(
  fit = c("unweighted", "weighted"),
  mean = c(mean(results$unweighted), mean(results$weighted)),
  bias = c(
    mean(results$unweighted - target_interaction),
    mean(results$weighted - target_interaction)
  ),
  rmse = c(
    sqrt(mean((results$unweighted - target_interaction)^2)),
    sqrt(mean((results$weighted - target_interaction)^2))
  )
)

print(summary_table, row.names = FALSE)
cat(
  "Weighted estimate was closer to the target in",
  mean(abs(results$weighted - target_interaction) <
    abs(results$unweighted - target_interaction)),
  "of replicates.\n"
)

stopifnot(
  abs(summary_table$bias[summary_table$fit == "weighted"]) <
    abs(summary_table$bias[summary_table$fit == "unweighted"]),
  summary_table$rmse[summary_table$fit == "weighted"] <
    summary_table$rmse[summary_table$fit == "unweighted"]
)
