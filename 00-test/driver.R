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
out1 = r_target_old(n, mu, tau, kappa, lambda, tol_suff, tol_merge, max_rejects)

hist(out1$draws)
out1$knots
out1$log_bounds
out1$rejections
out1$elapsed

# New code that uses vws package
out2 = r_target(n, mu, tau, kappa, lambda, tol_suff, tol_merge, max_rejects, report = 1000)

hist(out2$draws)
out2$knots
out2$log_bounds
out2$rejections
out2$elapsed

# TBD: Are we counting knots and regions differently? Make sure we can match
# results from manuscript!
plot(out1$knots)
plot(out2$regions)
