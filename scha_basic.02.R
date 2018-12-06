library(dlm)


build_dlm_scha <- function(order_trend,order_cycle,freq_names,rho) {
hours_per_cycle <- c(
4*2.*12.420612,
23.93447213,
12.4206012,
25.81933871,
12.,
26.868350,    #O1
12.65834751,
12.19162085,  #0.0820235525 S2
6.210300601, 
8.177140247,
4.140200401,
4.930880215,
12.4206012*2,
12.4206012*2/3.,
12.4206012*2/5.,
12.4206012*2/8.
) #0.2028035475

names(hours_per_cycle) <- c("LOW","K1","M2","O1","S2","Q1","N2","L2","M4","MK3","M6","MK5","M1","M3","M5","M8")
radians_per_hour <- 2.*pi/hours_per_cycle

freq_used <- 0.25*radians_per_hour[freq_names]
nfreq <- length(freq_used)

print("nfreq")
print(nfreq)
print("order_trend")
print(order_trend)
print("order_cycle")
print(order_cycle)

nstate <- order_trend + 2*order_cycle*nfreq

# Construct observation matrix
obs_trend <- c(1,rep(0,order_trend-1))
obs_per_freq <- c(1,0,rep(0,2*(order_cycle-1)))
FF <- t(as.matrix(c(obs_trend, rep(obs_per_freq,nfreq))))

# Construct state transition matrix
state_trend <- diag(order_trend)
if (order_trend>1) {
    for (itrend in 1:(order_trend-1)){state_trend[itrend,itrend+1] = 1}
}
# cycles 
rhovec = rep(rho,nfreq)

print("FF")
print(freq_used)
print(rho)

#cycles_matrices <- matrix(c(cos(x),sin(x),-sin(x),cos(x)),2)
#print(cycles_matrices)

single_freq <- function(f,rho){
    #cycle_identity <- rep(c(rep(1,(order_cycle-1)*2),0,0),nfreq)[1:(2*order_cycle*nfreq-2)]
    #cycle_lag_part <-lapply(1:4,function(x) diag(min(1,x%%order_cycle),2))

    cycle_trig_part <- rep(list(rho*matrix(c(cos(f),-sin(f),sin(f),cos(f)),2)),order_cycle)
    single <-bdiag(cycle_trig_part)

    for(i in 1:(2*(order_cycle-1))){
       single[i,i+2] = 1
    }

    single
}

#blah <-lapply(freq_used,function(x) bdiag(rep(list(matrix(c(cos(x),sin(x),-sin(x),cos(x)),2)),order_cycle)))
cycle_matrices <-mapply(single_freq,freq_used,rhovec,SIMPLIFY=FALSE)
#cycle_identity <- rep(c(rep(1,(order_cycle-1)*2),0,0),nfreq)[1:(2*order_cycle*nfreq-2)]


state_cycle <- bdiag(cycle_matrices)
print("state_cycle dimensions")
print(dim(state_cycle))


# combine trend and cycle components
GG <- bdiag(state_trend,state_cycle)
list(FF=FF,GG=GG)
}

##########################

hours_per_cycle <- c(
4*2.*12.420612,
23.93447213,
12.4206012,
25.81933871,
12.,
26.868350,    #O1
12.65834751,
12.19162085,  #0.0820235525 S2
6.210300601, 
8.177140247,
4.140200401,
4.930880215,
12.4206012*2,
12.4206012*2/3.,
12.4206012*2/5.,
12.4206012*2/8.
) #0.2028035475

freq_names <- c("LOW","K1","M2","O1","S2","Q1","N2","L2","M4","MK3","M6","MK5","M1","M3","M5","M8")
names(hours_per_cycle)  <- freq_names
radians_per_hour <- 2.*pi/hours_per_cycle
freq_used <- 0.25*radians_per_hour[freq_names]


###################

fname <- "F:/projects/stochastic_cycle/freeport_flow.csv"
scale<-10000.

y <- read.csv(fname, header = FALSE,stringsAsFactors=FALSE)
y <- as.matrix(y)
y <- y - mean(y)

ysave <- y
select <- 1200:3500

#remove <- 1850:2000
#remove <- 1250:1470
#remove <- 1650:1720
#remove <- 2050:2180
#remove <- 2250:2320
#remove <- 1650:1720
#remove <- 2450:2520
#remove <- 2600:2725

remove <- c(1290:1410,2230:2290,2600:2690)
remove <- c(1290:1410,2230:2430,2600:2690)
y[remove] <- NA

regress_names <- c("K1","M2","O1","S2","N2","L2","Q1","M4","MK3","M6","MK5")
regress_freq = freq_used[regress_names]

regress = TRUE
if (regress){
  t <- 1:length(y)
  X <- lapply(regress_freq, function(x) cbind(cos(x*t),sin((x*t))))
  X <- as.data.frame(X)
  lm1 <- lm(y~ K1.1+K1.2+M2.1+M2.2+O1.1+O1.2+S2.1+S2.2+Q1.1+Q1.2+N2.1+N2.2+L2.1+L2.2+M4.1+M4.2+MK3.1+MK3.2,na.action=na.exclude,data=X)
  yres <- residuals(lm1)
  y[t,] <- yres
  ylm <- predict(lm1,as.data.frame(X))
  rho <- 0.99   # 0.99, trend = 2 cycle = 2 works
} else {
  rho <- 0.97   # 0.99, trend = 2 cycle = 2 works
  ylm <- y     # dummy for reconstrution
  ylm[] <- 0
}

# 0.99, trend = 2 cycle = 2 works
# 0.97, trend = 3 cycle = 3 works

ordertrend = 2
ordercycle = 2
freqs <- c("K1","M2","MK3","M4","MK5","M6")
freqs <- c("LOW","M1","M2","M3","M4","M5","M6")
freqs <- c("M1","M2","M3","M4","M5")
freqs <- c("M1","M2","M3","M4","M5")
#freqs <- c("K1","O1","M2","MK3","M4")
freqs <- c("M1","M2","M3","M4","M6")

numfreq <- length(freqs)
determ <- build_dlm_scha( ordertrend, ordercycle, freqs, rho=rho)
V <- diag(c(1))/100.
Wtrend <- c(rep(0.,ordertrend-1),1.)
Wcycle <- rep( c(rep(0,2*(ordercycle-1)),1,1),numfreq)
W <- diag(c(Wtrend,Wcycle))/200.
W[3:6,3:6] <- W[3:6,3:6]/2
W[11:nrow(W),11:ncol(W)] <- W[11:nrow(W),11:ncol(W)]/4

m0 <- rep(0., dim(W)[1])
C0 <- diag(10000,dim(W)[1])

model = dlm(FF=determ$FF,V=V,GG=determ$GG,W=W,m0=m0,C0=C0)
modFilt <- dlmFilter(y[,1],model)
modSmooth <- dlmSmooth(modFilt)

ycompare <- ysave
ycompare[is.na(y)] <- NA

abline(h=0)
D0 = 1
D1 = 1+ordertrend
D2 = 1+ordertrend+ordercycle*2
D3 = 1+ordertrend+ordercycle*4
D4 = 1+ordertrend+ordercycle*6
D5 = 1+ordertrend+ordercycle*8

sumf <- modSmooth$s[-1,D0] + modSmooth$s[-1,D1] + modSmooth$s[-1,D2] + modSmooth$s[-1,D3] +modSmooth$s[-1,D4] +modSmooth$s[-1,D5]
yreconstruct <- sumf + ylm

plot.ts(ysave[select]*1,ylim=c(-15000,17000.),col="black",xlim=c(100,2400))
lines(ycompare[select],lwd=2,col="red")
lines(yreconstruct[select],col="green")
lines(sumf[select],col="brown")

lines(modSmooth$s[select,D0],col="brown")
lines(modSmooth$s[select,D1],col="pink")
lines(modSmooth$s[select,D2],col="orange")
lines(modSmooth$s[select,D3],col="chocolate")
lines(modSmooth$s[select,D4],col="aquamarine1")
lines(modSmooth$s[select,D5],col="aquamarine2")
lines(ylm[select],col="gray")




#####################################################################
#####################################################################


plot.ts(y[select],ylim=c(-4*scale,4.*scale),col="black")
lines(sumf[select],col="red")

ts.plot(modSmooth$s[,D0])
lines(modSmooth$s[,D2])
ts.plot(ylm)








#####################################################################################
#####################################################################################


lines(modSmooth$s[select,1])
lines(modSmooth$s[select,1+ordertrend],col="dark green")
lines(modSmooth$s[select,1+ordertrend+ordercycle*2],col="blue")
lines(modSmooth$s[select,1+ordertrend+ordercycle*4],col="purple")
lines(modSmooth$s[select,1+ordertrend+ordercycle*6],col="purple")
lines(modSmooth$s[select,27],col="cyan")
lines(modSmooth$s[select,35],col="brown")
lines(modSmooth$s[select,43],col="orange")
lines(modSmooth$s[select,51],col="darkslategray")
lines(modSmooth$s[select,1]+modSmooth$s[select,3],col="dark green")
plot.ts(fitted(lm1))
lines(sumf)
plot.ts(modFilt$f[select,3])

fname <- paste0(sta,'_',type,'_fit_data',case,'.csv')
write.csv(dropFirst(modSmooth$s),file=fname)






write_fit_data <- function(type,sta,case='',V,fname='') {
    if (fname == '') {
        fname <- paste0('D:/control_volume2/flow_data/kalman_filter_input/',sta,'_',type,'.dat')}  
  
    y <- read.csv(fname, header = FALSE)
    y <- as.matrix(y)
    fname <- paste0(sta,'fit_',type,case,'.dat')
    fit <- dget(fname)
    V <- V
    model <- buildFunM2O1K1S2Q1N2L2M4MK3M6MK5_v1(V,fit$par)
    modFilt <- dlmFilter(y[,2],model)
    modSmooth <- dlmSmooth(modFilt)
    fname <- paste0(sta,'_',type,'_fit_data',case,'.csv')
    write.csv(dropFirst(modSmooth$s),file=fname)
}


fitData <- function(infile,outfile,parm) {
    y <- read.csv(infile, header = FALSE)
    y = as.matrix(y)
    print(Sys.time())
    flush.console()
    fit <- dlmMLE(y[,2], control=list(maxit = 30000), parm = parm,
                    build = buildFunM2O1K1S2Q1N2L2M4MK3M6MK5)  
    print(warnings())
    flush.console()
    print(Sys.time())
    flush.console()
    dput(fit,file=outfile)
    return(fit)
}



write_fit_data <- function(type,sta,case='',V,fname='') {
    if (fname == '') {
        fname <- paste0('D:/control_volume2/flow_data/kalman_filter_input/',sta,'_',type,'.dat')}    
    y <- read.csv(fname, header = FALSE)
    y <- as.matrix(y)
    fname <- paste0(sta,'fit_',type,case,'.dat')
    fit <- dget(fname)
    V <- V
    model <- buildFunM2O1K1S2Q1N2L2M4MK3M6MK5_v1(V,fit$par)
    modFilt <- dlmFilter(y[,2],model)
    modSmooth <- dlmSmooth(modFilt)
    fname <- paste0(sta,'_',type,'_fit_data',case,'.csv')
    write.csv(dropFirst(modSmooth$s),file=fname)
}














