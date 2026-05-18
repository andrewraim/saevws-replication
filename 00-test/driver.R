library(vws)

source("../01-sim-cond/functions.R", chdir = TRUE)

set.seed(1234)

n = 5000
max_rejects = 1e6

mu = 0
lambda = 1
tol_suff = 0.70
tol_merge = 0.01
kappa = 10
tau = 0.5

# Original code
out1 = r_target(n, mu, tau, kappa, lambda, tol_suff, tol_merge, max_rejects)

hist(out1$draws)
out1$regions
out1$log_bounds
out1$rejects
out1$elapsed

# New code that uses vws package
out2 = r_target_new(n, mu, tau, kappa, lambda, tol_suff, tol_merge, max_rejects, report = 1000)

hist(out2$draws)
out2$regions
out2$log_bounds
out2$rejects
out2$elapsed

# TBD: Are we counting knots and regions differently? Make sure we can match
# results from manuscript!
plot(out1$regions)
plot(out2$regions)

plot(cumsum(out1$rejects) / (cumsum(out1$rejects) + cumsum(seq_len(n))), ylim = c(0,0.01), type = "l")
points(cumsum(out2$rejects) / (cumsum(out2$rejects) + cumsum(seq_len(n))), ylim = c(0,0.01), type = "l")
