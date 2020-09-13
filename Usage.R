############################ tviGARCH #################################
set.seed(1)
n=500;
c1=3;
c2=1*0.1;
c3=-4*0.05;

c01=5#0.8;
c02=0.2;
c03=0.12;

resolution=n #means we compute at every 1/n

g=(1:resolution)/resolution;

cilmu=numeric(0);cila1=numeric(0);cilb1=numeric(0)
accumu=numeric(0);accua1=numeric(0);accub1=numeric(0)

a0 <- 25*exp(-(g-0.5)^2/0.1)
a1 <- 0.3*(g-1)^2+0.1
b1 <- 1-a1

e=rnorm(n);x=e;sigma2=array(0,n)
sigma20=a0[1]/(1-b1[1]);x0=rpois(1, sigma20);
sigma2[1]=a0[1]+a1[1]*x0+b1[1]*sigma20;x[1]=rpois(1, sigma2[1])

for (i in 2:n)
{
  sigma2[i]=a0[i]+a1[i]*x[i-1]+b1[i]*sigma2[i-1]
  x[i]= rpois(1, sigma2[i])
}  

fitS <- fit.tvINGARCHMCMCcombo(as.numeric(c(x0,x)), order1 = 1, order2 = 1, norder = 4, knot = 6)#, sigma2in = sigma2[1], sigup = T)

############################ tvGARCH #################################
n=500;
print(n);
print(seedc)

resolution=n #means we compute at every 1/n

g=(0:resolution)/resolution;
#a1=a1fun(g);a0=a0fun(g);b1=b1fun(g);


a0= -4*sin(0.5*pi*g)+5
a1 = -1*(g-0.3)^2 + 0.5
b1 = 0.2*sin(1*pi*g)+0.2

e=rnorm(n);x=e;sigma2=array(0,n)
sigma20=a0[1]/(1-b1[1]);x0=rnorm(1, 0, sqrt(sigma20));
sigma2[1]=a0[1]+a1[1]*x0^2+b1[1]*sigma20;x[1]=rnorm(1, 0, sqrt(sigma2[1]))

for (i in 2:n)
{
  sigma2[i]=a0[i]+a1[i]*x[i-1]^2+b1[i]*sigma2[i-1]
  x[i]= rnorm(1, 0, sqrt(sigma2[i]))
}  
#data <- as.numeric(c(x0,x))
data=x

fitS <-  fit.tvGARCHMCMCcomboN(data, order1 = 1, order2=1, knot = 5, norder = 4, P=10, Total_itr = itr, burn = brn)

################################ tvARCH #######################################
Ti <- n
t <- 0:Ti/Ti
mut0 <- 10*exp(-(t-0.5)^2/0.1)
at10 <- 0.4*(t-0.15)^2 + 0.1
datavect <- rep(0, Ti+1)
datavect[1] <- rnorm(1, 0, sqrt(mut0[1]))
datavect[2] <- rnorm(1, 0, sqrt(mut0[2]+at10[2]*datavect[1]^2))
for(i in 3:Ti){
  datavect[i] <- rnorm(1, 0, sqrt(mut0[i]+at10[i]*datavect[i-1]^2+at20[i]*datavect[i-2]^2))
}

data <- datavect

fitS <- fit.tvARMCMCunbd(data, order = 1, knot = 5, norder = 4, Total_itr = 10000, burn = 5000)
