rm(list=ls()) # clear the global environment

period = 30           # 1 year period
nbsteps = period*250  # number of steps in terms of work days
step = period/nbsteps # numerical representation of a step
timesequence = seq(0,period,length.out = nbsteps) #time range (x axis)

# Simulation of a Vasicek/Ornstein-Uhlenbeck process (i.e. the equivalent of an AR(1) process with mean reversing behavior)
a = 0.1          # convergence speed factor
b = 0.0276       # convergence limit (fundamental value of the spot rate)
sigma = 0.005   # daily instantaneous volatility or randomness impacting the evolution of the interest rate as proportion of the long-term value (this is a modification of Vasicek's model)
r_0 =  0.01      # initial value of the interest rate
r_t = timesequence #initialization of r_t
r_t[1] = r_0     # setting of the initial value of the vector
epsilon = rnorm(nbsteps)   #vector of 300 normally distributed values
#w = cumsum(epsilon)   # vector containing the state of the Wiener process at each time value


#Vasicek
for (i in 2:nbsteps){ #realization of the exact solution of the Vasicek model
  r_t[i]=r_t[i-1]*exp(-a*step)+b*(1-exp(-a*step))+sigma*sqrt((1-exp(-2*a*step))/2*a)*epsilon[i]
}         
plot(timesequence,r_t,type='l',col='red', ylab = "Interest Rate", xlab = "Years to maturity") # plot of the convergence process


