library(tidyverse)

set.seed(1234)

m = 2000
X = cbind(1, rnorm(m))
sigma2 = rgamma(m, 1.25, 1/20)
sigma = sqrt(sigma2)

beta_true = c(1, -1)
Xbeta_true = X %*% beta_true
tau_true = 1.25
mu_true = rlnorm(m, Xbeta_true, tau_true)
y = rnorm(m, mu_true, sigma)

data.frame(mu = mu_true, y = y) %>%
	ggplot() +
	geom_point(aes(mu, y)) +
	theme_minimal()

unmatch = tibble(y, sigma, x = X[,2])
write_csv(unmatch, file = "unmatch.csv")
