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
  y[i] <- A(r[i], sigma)
}
df <- data.frame(r,y)


### Add noise
n <- 1e3 #number of points to sample
sample_indices <- sample.int(N, size = n)
sample_indices <- sample_indices[order(sample_indices)]
sample_df <- df[sample_indices,]
plot(sample_df)

#Estimate variance
gamma <- 20*sample_df$r
z <- rbeta(n, gamma*sample_df$y, gamma*(1-sample_df$y))
epsilon <- rnorm(n, 0, sd = 1e-6)
sample_df$y_noise <- D + (1-D)*z + epsilon

plot(sample_df$y_noise ~ sample_df$r, cex = 0.5, pch = 20)
lines(x = sample_df$r, y = D + (1-D)*sample_df$y, col = 'red', lwd = 2)


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
  real<lower=0> sigma; //fuzzy nbhd sd
  real<lower=0,upper=1> D; //prop of genome at which each pair of chroms differ
  real<lower=0> eta; //normal noise
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
  epsilon ~ normal( 0, eta );
  eta ~ exponential(1e3);
  for(n in 1:N){
    Z[n] ~ beta( 20*r[n]*B( r[n], sigma ), 20*r[n]*(1 - B( r[n], sigma ));
  }
}
'

model_fit <- stan(model_code=model_code,
                  data = list(N = n,
                              r = sample_df$r,
                              L = L,
                              Y = sample_df$y_noise))
summary(model_fit, pars = c('sigma', 'D', 'eta'))
stan_hist(model_fit, pars = c('sigma', 'D', 'eta'))
stan_trace(model_fit)
saveRDS(model_fit, file = "model.rds")