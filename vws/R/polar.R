#' Polar transformation
#'
#' The polar transformation and its inverse, as expressed in Section 2.7 of
#' Anderson (2003), restricted to the unit sphere.
#'
#' @param u A point in the unit cube \eqn{[0,1]^{d-1}}.
#' @param z A point in the unit sphere
#' \eqn{\mathbb{S}^d = \{ \bm{z} \in \mathbb{R}^d: \bm{z}^\top \bm{z} = 1 \}}.
#'
#' @details
#' The `inv_polar` transformation `z = inv_polar(u)`, from
#' \eqn{[0,1]^{d-1}} to \eqn{\mathbb{S}^d}, can be described as a
#' composition of two transformations. The result
#' \deqn{
#' z_j =
#' \begin{cases}
#' \sin \phi_j \prod_{\ell=1}^{j-1} \cos \phi_\ell, & j = 1, \ldots, d-1, \\
#' \prod_{\ell=1}^{d-1} \cos \phi_\ell & j = d,
#' \end{cases}
#' }
#' is computed from angles
#' \deqn{
#' \phi_j =
#' \begin{cases}
#' \pi u_j - \pi / 2 , & j = 1, \ldots, d-2, \\
#' 2 \pi u_j - \pi, & j = d-1.
#' \end{cases}
#' }
#'
#' The `polar` transformation `u = polar(z)`, from \eqn{\mathbb{S}^d}
#' to \eqn{[0,1]^{d-1}}, is the inverse of `inv_polar`:
#' \deqn{
#' u_j =
#' \begin{cases}
#' (\phi_j + \pi / 2) / \pi, & j = 1, \ldots, d-2, \\
#' (\phi_j + \pi) / (2\pi), & j = d-1,
#' \end{cases}
#' }
#' and
#' \deqn{
#' \phi_j =
#' \begin{cases}
#' \arcsin \left\{ z_j \Big/ \sqrt{1 - z_1^2 - \cdots - z_{j-1}^2} \right\}, & j = 1, \ldots, d-1, \\
#' \arctan\{ z_{d-1} / z_d \}, & j = d.
#' \end{cases}
#' }
#'
#' @name polar
#' @examples
#' x = 1:5 / sum(1:5)
#' z = inv_polar(x)
#' polar(z)
#'
#' @references T.W. Anderson (2003). An Introduction to Multivariate
#' Statistical Analysis, Wiley-Interscience, 3rd Edition.
#'
NULL

#' @name polar
#' @export
polar = function(z)
{
	d = length(z)
	phi = numeric(d-1)

	for (j in seq_len(d-2)) {
		idx = seq_len(j-1)
		phi[j] = asin( z[j] / sqrt(1 - sum(z[idx]^2)) )
	}
	phi[d-1] = atan2(z[d-1], z[d])

	# Unshift and unscale
	idx = seq_len(d-2)
	u = numeric(d-1)
	u[idx] = (phi[idx] + pi/2) / pi
	u[d-1] = (phi[d-1] + pi) / (2*pi)

	return(u)
}

#' @name polar
#' @export
inv_polar = function(u)
{
	d = length(u) + 1
	stopifnot(all(0 <= u & u <= 1))

	# Compute angles phi from x.
	# Shift and scale x_1, ... x_{d-2} to be between -pi/2 and pi/2.
	# Shift and scale x_{d-1} to be between -pi and pi.
	idx = seq_len(d-2)
	phi = numeric(d-1)
	phi[idx] = pi * u[idx] - pi/2
	phi[d-1] = 2*pi * u[d-1] - pi

	z = numeric(d)
	for (j in seq_len(d-1)) {
		idx = seq_len(j-1)
		z[j] = sin(phi[j]) * prod(cos(phi[idx]))
	}
	z[d] = prod(cos(phi))

	return(z)
}

