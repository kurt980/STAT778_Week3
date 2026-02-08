# run a simulation for 1 df and save result to a single result file

# parse args: options + last positional task id
args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

# options
R <- as.integer(get_opt("--repl", "1000"))
n <- as.integer(get_opt("--n", "100"))

# last arg = task id
task_id <- as.integer(tail(args, 1))
if (is.na(task_id)) stop("missing task id as last command-line argument")

# fixed parameters
p <- 3
mu0 <- c(-0.7, 0, 0.7)
mu1 <- -mu0
alpha <- 0.05
R <- 1000
n <- 100

# map array task id -> df
dfs <- c(2, 3, 4, 5, 10, 100)
if (task_id < 1 || task_id > length(dfs)) stop("task id out of range")


df <- dfs[task_id]

set.seed(1000 + df)

################################################################################
# generate test sets once
X0_test <- matrix(rt(10000 * p, df = df), ncol = p) + rep(mu0, each = 10000)
X1_test <- matrix(rt(10000 * p, df = df), ncol = p) + rep(mu1, each = 10000)

typeIerrors <- numeric(R)
typeIIerrors <- numeric(R)

phi <- qnorm(1 - alpha)
#-------------------------------------------------------------------------------
for (r in seq_len(R)) {
  set.seed(1e6 * df + r)
  
  # generate train sets
  X0_train <- matrix(rt(n * p, df = df), ncol = p) + rep(mu0, each = n)
  X1_train <- matrix(rt(n * p, df = df), ncol = p) + rep(mu1, each = n)
  
  # estimate parameters
  mu0_hat <- colMeans(X0_train)
  mu1_hat <- colMeans(X1_train)
  Sigma0_hat <- stats::cov(X0_train)
  Sigma1_hat <- stats::cov(X1_train)
  
  # pool covariances
  Sigma_hat <- ((n - 1) * Sigma0_hat + (n - 1) * Sigma1_hat) / (2 * n - 2)
  
  # group some terms
  delta_hat <- mu1_hat - mu0_hat
  beta_hat <- solve(Sigma_hat, delta_hat)
  quad_term <- drop(t(delta_hat) %*% beta_hat)
  
  # calculate simple version of bound
  simple_bound <- phi * sqrt(quad_term) + drop(t(beta_hat) %*% mu0_hat)
  
  typeIerrors[r] <- mean(drop(X0_test %*% beta_hat) > simple_bound)
  typeIIerrors[r] <- mean(drop(X1_test %*% beta_hat) <= simple_bound)
}

out <- data.frame(
  df = df,
  mean_typeI = mean(typeIerrors),
  sd_typeI = sd(typeIerrors),
  mean_typeII = mean(typeIIerrors),
  sd_typeII = sd(typeIIerrors),
  R = R,
  n = n,
  p = p
)

dir.create("sim_out", showWarnings = FALSE)
saveRDS(out, file.path("sim_out", sprintf("result_df_%s_task_%02d.rds", df, task_id)))

cat("saved:", outfile, "\n")
print(out)
