#packages
library(rstan)

#parameters
N <- 1e4 #number of points to plot
L <- 10 #r in (0,L)
r <- runif(N, 0, L)
sigma <- 3 #sd of X and Y
D <- 0.9999323


### Compute expectation as integral over [0,1]
integrand <- function(t, r, sigma){
  result <- exp(-(r^2 * t / sigma^2 + sigma^2 / t)/2)/t
  return(result)
}

A <- function(r, sigma){
  result <- integrate(integrand, r=r, sigma=sigma, lower = 0, upper = 1)$value
  result <- exp(sigma^2/2) * result / 2
  return(result)
}

r <- r[order(r)]
y<-0
for(i in 1:length(r)){
  y[i] <- (1-D)*A(r[i], sigma) + D
}
df <- data.frame(r,y)


### Add noise
n <- 100 #number of points to sample
sample_indices <- sample.int(N, size = n)
sample_indices <- sample_indices[order(sample_indices)]
sample_df <- df[sample_indices,]
plot(sample_df)


#Estimate variance
z <- rbeta(n, (10-sample_df$r)*.05,5*sample_df$r)*(1-D)
epsilon <- rnorm(n, sample_df$y, sd = (10 - sample_df$r)*5e-7 + sample_df$r*1e-7)
sample_df$y_noise <- z + epsilon

plot(sample_df$y_noise ~ sample_df$r, ylim = c(0.99992,1))
lines(x = sample_df$r, y = sample_df$y, col = 'red')


### RSTAN model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)
model_code <- '
functions {
  real integrand(real t,real r, real sigma){
    real result;
    result = exp(-(sigma^2/t + r^2 * t/sigma^2)/2)/t;
    return result;
  }

  real B(real r, real sigma){
    real h;
    real x;
    real result;
    h = 1e-2;
    x = h;
    result = h * integrand(h,r,sigma);
    while(x < 1){
      result += h * (integrand(x, r, sigma) + integrand(x+h, r, sigma)) / 2;
      x += h;
    }
    result *= exp(sigma^2/2)/2;
    return result;
  }
}
data {
  int N; //number of observations
  real r[N]; //geographic distances
  real L; //max r value
  real Y[N]; //proportion identical
}
parameters {
  real<lower=0> sigma;
  real<lower=0,upper=1> D;
  real<lower=0> alpha[3];
  real<lower=0, upper=1> Z[N];
}
transformed parameters {
  real epsilon[N];
  for(n in 1:N){
    epsilon[n] = Y[n] - (1-D)*Z[n] - D;
  }
}
model {
  sigma ~ uniform(0,10);
  D ~ beta(5,1);
  alpha ~ exponential(1000);
  for(n in 1:N){
    Z[n] ~ beta( (L - r[n])*alpha[1], L*r[n]/2 );
    epsilon[n] ~ normal( B(r[n], sigma), (L - r[n])*alpha[2] + r[n]*alpha[3] );
  }
}
'

model_fit <- stan(model_code=model_code,
                  data = list(N = n,
                              r = sample_df$r,
                              L = L,
                              Y = sample_df$y_noise))
saveRDS(model_fit, file = "model.rds")