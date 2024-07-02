# Effort
e1 <- 1
e2 <- 10

# Catchability
q1 <- 1
q2 <- 1

# Selectivity at length
s1 <- c(0.05, 0.1, 0.3, 0.5, 0.9)
s2 <- c(0.1, 0.25, 0.4, 0.7, 0.5)

plot(s1, type = 'l')
lines(s2, col = "red")

# Numerical density available to a unit of sampling effort
n1 <- rep(10, 5)
n2 <- rep(10, 5)

# Selectivity ratio
S12 <- s1/s2

# Catch - If this is in terms of CPUE, effort is used as a normalizing factor
c1 <- (q1*s1*n1)*e1
c2 <- (q2*s2*n2)*e2

# Catch comparison rate
p12 <-  (c1/e1)/(c1/e1+c2/e2)

round((c1/e1)/(c2/e2), 15) == round(p12/(1-p12), 15)

# Proportion at length
r1 <- c1/sum(c1)
r2 <- c2/sum(c2)

p12_prop <- r1/(r1+r2)

# Selectivity ratio - 
r12 <- (r1/r2)
S12_a <- p12/(1-p12) # Absolute selectivity ratio
S12_r <- p12_prop/(1-p12_prop) # Relative selectivity ratio

round(p12_prop/(1-p12_prop), 15)/S12

# Check that selectivity ratios are equal
S12 - S12_a

# Conversion; note that this still requires effort comparison
c1_fit <- S12_a * c2 * (e1/e2)

c1 - c1_fit

c2_fit <- 1/S12_a * c1 * (e2/e1)

# Undefined for the smallest length since gear 1 has selectivity 0 for the smallest size
c2 - c2_fit

S12_r/S12

(q1*e1)/(q2*e2)*mean(c2)/mean(c1)

mean(s2)/mean(s1)


# Webster et al. (2020)
rc1 <- cumsum(c1)/sum(c1)
rc2 <- cumsum(c2)/sum(c2)

rc_test1 <- cumsum(s1*n1)/sum((s1*n1))
rc_test2 <- cumsum(s2*n2)/sum((s2*n2))

rc1/rc2

