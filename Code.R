rm(list=ls()) # clears the global environment

period = 30           # 1 year period
T = period*250        # number of steps in terms of work days
step = period/T # numerical representation of a step
timesequence = seq(0,period,length.out = T) #time range (x axis)

# Simulation of a Vasicek/Ornstein-Uhlenbeck process (i.e. the equivalent of an AR(1) process with mean reversing behavior)
a = 0.1          # convergence speed factor
b = 0.0276       # convergence limit (fundamental value of the spot rate)
sigma = 0.005    # annualized instantaneous volatility 
r_0 =  0.01      # initial value of the interest rate
r_t = timesequence #initialization of r_t
r_t[1] = r_0     # setting of the initial value of the vector
epsilon = rnorm(T)   #vector of 300 normally distributed values
#w = cumsum(epsilon)   # vector containing the state of the Wiener process at each time value


# Vasicek term structure simulation
for (i in 2:T){
  r_t[i]=r_t[i-1]*exp(-a*step)+b*(1-exp(-a*step))+sigma*sqrt((1-exp(-2*a*step))/2*a)*epsilon[i]
}         
plot(timesequence,r_t,type='l',col='red', ylab = "Interest Rate", xlab = "Years to maturity",main = "Vasicek model realization") 

# CIR term structure simulation
for (i in 2:T){ #realization of the exact solution of the CIR model
  r_t[i]=r_t[i-1]*exp(-a*step)+b*(1-exp(-a*step))+sigma*sqrt((1-exp(-2*a*step))/2*a)*sqrt(r_t[i-1])*epsilon[i]
}         
plot(timesequence,r_t,type='l',col='red', ylab = "Interest Rate", xlab = "Years to maturity",main = "CIR model realization") 


B_t_T=timesequence
A_t_T=timesequence
P_t_T=timesequence

# Vasicek price of a discount bound
for (i in 2:T){
  B_t_T[i] = (1-exp(-a*(T-timesequence[i])))/a  #due to the high number of steps, the exponential term is practically 0
  A_t_T[i] = exp((b-sigma^2/(2*a^2))*(B_t_T[i]-30+timesequence[i])-sigma^2*B_t_T[i]^2/4*a)
  P_t_T[i] = A_t_T[i]*exp(-r_t[i]*B_t_T[i])
}
P_t_T[1]=P_t_T[2]
plot(timesequence,P_t_T,type='l',col='red', ylab = "Price of a discount bond", xlab = "Years to maturity",main = "Discount bond - Vasicek model") 

# CIR price of a discount bound
h=sqrt(a^2+2*sigma^2)
for (i in 2:T){
  B_t_T[i] = 2*exp(h*(T-timesequence[i])-1)/(2*h+(a+h)*exp(h*(T-timesequence[i])-1))
  A_t_T[i] = (2*h*exp((a+h)*(T-timesequence[i])/2)/(2*h+(a+h)*exp(h*(T-timesequence[i])-1)))^(2*a*b/sigma^2)
  P_t_T[i] = A_t_T[i]*exp(-r_t[i]*B_t_T[i])
}
P_t_T[1]=P_t_T[2]
plot(timesequence,P_t_T,type='l',col='red', ylab = "Price", xlab = "Years to maturity",main = "Discount bond - CIR model") 

