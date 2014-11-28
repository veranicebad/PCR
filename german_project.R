panel.lmpolyline <- function(x, y, groups = NULL,
                             degree = 1, col.line = par.line$col,
                             lty = par.line$lty, lwd = par.line$lwd,
                             alpha = par.line$alpha, ..., identifier = "lmpolyline") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (!is.null(groups)) {
    par.line <- trellis.par.get("superpose.line")
    panel.superpose(x = x, y = y, groups = groups,
                    degree = degree, col.line = col.line,
                    lty = lty, lwd = lwd, alpha = alpha,
                    panel.groups = sys.function(),
                    ...)
  } else {
    if (length(x) > degree) {
      l.par.get("plot.line")
      panel.curve(predict(l, list(x = x)),
                  from = min(x), to = max(x),
                  col.line = col.line, lty = lty,
                  lwd = lwd, alpha = alpha,
                  ..., identifier = identifier)
    }
  }
}


xyplot(y ~ x, data = d) + layer_(panel.lmpolyline(...,col.line = 134, degree = 20))



fit_realdata<-function()
{
  mydata <- read.table("c:/Users/Vera/Documents/������/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
  den<-density(as.numeric(mydata[mydata$Target == "AMPL661570029",]))
  d<-data.frame(x=den$x,y=den$y)
  print(d)
  #plot(density(as.numeric(mydata[mydata$Target == "AMPL661570029",])), xlab = "Number of amplicons", ylab = "density", main = "Density plot of real data")
  fit<-lm(formula = d$x ~ poly(d$y, 4))
  #plot(data$x, data$y, type="p", lwd=3)
  pol <- function(x)
  {
    a<-0
    for(i in 1:4)
    {
      a<-a+fit$coefficient[i+1]*x^(i)
    }
    return(a)
  }
  curve(pol,-500, 500, ylim=c(0, 0.001), col="red", lwd=2)
  points(data$x, data$y, type="p", lwd=3)
  #abline(fit, col="red")
  #abline(fit)
}

plot_real_data<-function()
{
  #jpeg('rplot1.jpg')
  mydata <- read.table("c:/Users/Vera/Documents/������/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
  #plot(mydata$IonXpress_14_015, mydata$IonXpress_14_016)
  #plot(as.numeric(mydata[mydata$Target == "AMPL1380948436",]), as.numeric(mydata[mydata$Target == "AMPL655136916",]))
  plot(density(as.numeric(mydata[mydata$Target == "AMPL661570029",])), xlab = "Number of amplicons", ylab = "density", main = "Density plot of real data")
  lst1<-as.numeric(mydata[mydata$Target == "AMPL661570029",])
  lst2<-as.numeric(mydata[mydata$Target == "AMPL655136916",])
  YZ<-data.frame(y=lst1,z=lst2)
  sort_yz<-YZ[ order(YZ[1]),]
  print(sort_yz)
  return(sort_yz)
  #get_linear_model(sort_yz)
  #dev.off()
}

remove_deletions_from_real_data<-function(){
  mydata <- read.table("c:/Users/Vera/Documents/������/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
  mydata$IonXpress_17_013<-NULL
  mydata$IonXpress_17_021<-NULL
  mydata$IonXpress_17_045<-NULL
  mydata$IonXpress_18_042<-NULL
  mydata$IonXpress_18_044<-NULL
  mydata$IonXpress_18_046<-NULL
  mydata$IonXpress_18_001<-NULL
  mydata$IonXpress_19_037<-NULL
  mydata$IonXpress_20_012<-NULL
  mydata$IonXpress_20_016<-NULL
  return(mydata)
}

plot_parAB_real_data<-function(){
  mydata<-remove_deletions_from_real_data()
  run_20<-mydata[, grep("IonXpress_20", colnames(mydata))]
  run_17<-mydata[, grep("IonXpress_17", colnames(mydata))]
  run_18<-mydata[, grep("IonXpress_18", colnames(mydata))]
  run_19<-mydata[, grep("IonXpress_19", colnames(mydata))]
  run_16<-mydata[, grep("IonXpress_16", colnames(mydata))]
  run_15<-mydata[, grep("IonXpress_15", colnames(mydata))]
  run_14<-mydata[, grep("IonXpress_14", colnames(mydata))]
  plot(0,0,col='white',xlim=c(0,15),ylim=c(0,15))
  for(i in 1:nrow(mydata)){
    amp_i_run_14 <- as.numeric(run_14[i,])
    amp_i_run_15 <- as.numeric(run_15[i,])
    amp_i_run_16 <- as.numeric(run_16[i,])
    amp_i_run_17 <- as.numeric(run_17[i,])
    amp_i_run_18 <- as.numeric(run_18[i,])
    amp_i_run_19 <- as.numeric(run_19[i,])
    amp_i_run_20 <- as.numeric(run_20[i,])
    ab_amp_i_run_14<-nsamplGibs(1, amp_i_run_14)
    ab_amp_i_run_15<-nsamplGibs(1, amp_i_run_15)
    ab_amp_i_run_16<-nsamplGibs(1, amp_i_run_16)
    ab_amp_i_run_17<-nsamplGibs(1, amp_i_run_17)
    ab_amp_i_run_18<-nsamplGibs(1, amp_i_run_18)
    ab_amp_i_run_19<-nsamplGibs(1, amp_i_run_19)
    ab_amp_i_run_20<-nsamplGibs(1, amp_i_run_20)
    points(ab_amp_i_run_14[1],ab_amp_i_run_14[2], pch=16, col="red")
    points(ab_amp_i_run_15[1],ab_amp_i_run_15[2], pch=16, col="orange")
    points(ab_amp_i_run_16[1],ab_amp_i_run_16[2], pch=16, col="yellow")
    points(ab_amp_i_run_17[1],ab_amp_i_run_17[2], pch=16, col="green")
    points(ab_amp_i_run_18[1],ab_amp_i_run_18[2], pch=16, col="cyan")
    points(ab_amp_i_run_19[1], ab_amp_i_run_19[2],pch=16, col="blue")
    points(ab_amp_i_run_20[1],ab_amp_i_run_20[2], pch=16, col="violet")
  }
}

nsamplGibs<-function(n, data){
  library(e1071) 
  mean_data<-mean(data)
  std_data<-sd(data)  
  skewness_data<-skewness(data)
  kurtosis_data<-kurtosis(data)
  a_min=0
  b_min=0
  a_max=15
  b_max=15
  a_ans=0
  b_ans=0
  best_r=1000
  
  for(i in 1:n){
    a0=runif(1,min=a_min, max=a_max)
    b0=runif(1,min=b_min, max=b_max)
    ab<-samplGibs(a0,b0,data)
    data_sample=get_Y_2(ab[1],ab[2])
    mean_sample<-mean(data_sample)    
    sd_sample<-sd(data_sample)
    skewness_sample<-skewness(data_sample)
    kurtosis_sample<-kurtosis(data_sample)
    r=abs(mean_data-mean_sample)+abs(std_data-sd_sample)+abs(skewness_data-skewness_sample)+abs(kurtosis_data-kurtosis_sample)+as.numeric(quantile(data_sample,0.025))+as.numeric(quantile(data_sample,0.975))
    if(r<=best_r)
    {
      best_r <- r
      a_ans=ab[1]
      b_ans=ab[2]
    }
  }
  ans<-c()
  ans<-c(ans,a_ans)
  ans<-c(ans,b_ans)
}

samplGibs<-function(a0,b0,data)
{
  n=100
  m=100
  mean_data<-mean(data)
  std_data<-sd(data)
  skewness_data<-skewness(data)
  kurtosis_data<-kurtosis(data)
  best_r=1000
  r = 1000
  a_min=0
  b_min=0
  a_max=15
  b_max=15
  a_ans=0
  b_ans=0
  a=a0
  b=b0
    for(j in 1:n){
      print(best_r)
      print(j)
      if(j%%2==0){
        prev_a = a
        a=runif(1,min=a_min, max=a_max)
      }
      else{
        prev_b = b
        b=runif(1,min=b_min, max=b_max)
      }
      mean_sample<-0
      sd_sample<-0
      skewness_sample<-0
      kurtosis_sample<-0
      for(i in 1:m){
        data_sample=get_Y_2(a,b)
        mean_sample<-mean_sample+mean(data_sample)
        sd_sample<-sd_sample+sd(data_sample)
        skewness_sample<-skewness_sample+skewness(data_sample)
        kurtosis_sample<-kurtosis_sample+kurtosis(data_sample)
      }
      mean_sample<-mean_sample/m
      sd_sample<-sd_sample/m
      skewness_sample<-skewness_sample/m      
      kurtosis_sample<-kurtosis_sample/m
      r_pr<-r
      r=abs(mean_data-mean_sample)+abs(std_data-sd_sample)+abs(skewness_data-skewness_sample)+abs(kurtosis_data-kurtosis_sample)+as.numeric(quantile(data_sample,0.025))+as.numeric(quantile(data_sample,0.975))
      if(r<=best_r)
      {
        best_r <- r
        a_ans=a
        b_ans=b
      }
      k=r_pr/(r+r_pr)
      t<-runif(1,0,1)
      if(t<=k)
      {
        if(j%%2==0){
          a = prev_a
        }
        else
        {
          b = prev_b
        }
      }
    
  }
  ans<-c()
  ans<-c(ans,a_ans)
  ans<-c(ans,b_ans)
}

get_Y_2<-function(a,b)
{
  #jpeg('rplot.jpg')
  #p<-0.95
  lst1 <- c()
  y0<-1
  for (i in 1:100)
  {
    scale <- rnorm(1, 0, 0.65) + 2
    y <- y0* scale #+ rnorm(1, 0, 0.005)
    for (j in 1:12)
    { 
      p<-1/(1+exp((j-a)*b))
      t<-runif(1,0,1)
      #print(t)
      if(t<=p)
      {
        eff<-2
      }
      else
      {
        eff<-1
      }
      #print(eff)
      y <- y * eff
    }
    lst1 <- c(lst1, y)
  }
  Y<-data.frame(y=lst1)
  sort_y<-Y[order(Y[1]),]
  d<-density(lst1)
  #plot(d, xlab = "Number of amplicons", ylab = "density", main = "Density plot of simulated data")
  #dev.off()
  return(sort_y)
}

get_eff<-function(Y)
{
  maxCor <- 0
  eff_value <- 1
  for(i in seq(1, 2, 0.01))
  {
    eff<-i
    #print(eff)
    y0<-1
    if(maxCor < cor(Y, get_Y(y0,eff)))
    {
      maxCor <- cor(Y, get_Y(y0,eff))
      eff_value <- eff
    }
    #print(cor(Y, get_Y(y0,eff)))
  }
  ans<-c(maxCor, eff)
  print(ans)
  return(ans)
}
#   varY_fin=var(Y)
#   varY=y0*0.65+0.005
#   for (i in 1:48)
#   {
#     for (j in 1:12)
#     { 
#       varK=Y[i]/7000     
#       varYdeleff=varY*varK/(varY+varK)
#       varY=varYdeleff*eff
#     }
#   }
#   
#   return(eff)
#}

plot_density_xy<- function(lst)
{
  lst3<-lst[[1]]*lst[[2]]
  plot(density(lst3))
}

plot_density_x<- function(lst)
{  
  jpeg('rplot2.jpg')
  plot(density(lst[[1]]), xlab ="Number of amplicon 1", ylab="Density", main = "Simulated data")
  dev.off()
}

get_Y_as_in_lect <- function(y0, p)
{
  lst1 <- c()
  for (i in 1:48)
  {
    #scale <- rnorm(1, 0, 0.65) + 3
    y <- y0
    for (j in 1:12)
    {
      xn <- rnorm(1, 0, p*(1-p)*y)
      y <- (1+p)*y+xn
    }
    lst1 <- c(lst1, y)
  }
  Y<-data.frame(y=lst1)
  print(mean(lst1))
  sort_y<-Y[order(Y[1]),]
  d<-density(lst1)
  plot(d, xlab ="Number of amplicon 1", ylab="Density", main = "Simulated data")
  #print(sqrt(var(sort_y)))
  return(lst1)
}

get_YZ_as_in_lect <- function(x, p1, p2)
{
  lst1 <- c()
  lst2 <- c()
  for (i in 1:480)
  {
    y <- x
    z <- x
    
    for (j in 1:12)
    {
      yn <- rnorm(1, 0, p1*(1-p1)*y)#p1*(1-p1)*
      y <- (1+p1)*y+yn
      zn <- rnorm(1, 0, p2*(1-p2)*z/7000)#p2*(1-p2)*
      z <- (1+p2)*z+zn
    }
      lst1 <- c(lst1, y)
      lst2 <- c(lst2, z)
  }
  d<-density(lst1)
  #plot(d)
  plot(lst1, lst2, xlim=c(0,3000), ylim=c(0,3000))
  YZ<-data.frame(y=lst1,z=lst2)
  sort_yz<-YZ[ order(YZ[1]),]
  return(sort_yz)
}

get_Y <- function(y0, eff1)
{
  lst1 <- c()
  for (i in 1:48)
  {
    scale <- rnorm(1, 0, 0.65) + 3
    y <- y0 * scale + rnorm(1, 0, 0.005)
    for (j in 1:12)
    {
      k <- rnorm(1, 1, y/7000)
      y <- y * (eff1 * k)
    }
    lst1 <- c(lst1, y)
  }
  Y<-data.frame(y=lst1)
  sort_y<-Y[order(Y[1]),]
  d<-density(lst1)
  plot(d, xlab ="Number of amplicon 1", ylab="Density", main = "Simulated data")
  #print(sqrt(var(sort_y)))
  return(lst1)
}

get_YZ <- function(x, eff1, eff2)
{
  #print(x)
  #print(eff1)
  #print(eff2)
  lst1 <- c()
  lst2 <- c()
  for (i in 1:48)
  {
    scale <- rnorm(1, 0, 0.65) + 3
    #print(scale)
    y <- x * scale + rnorm(1, 0, 0.005)
    z <- x * scale + rnorm(1, 0, 0.005)
    
    for (j in 1:12)
    {
      k <- rnorm(1, 1, y/7000)
      y <- y * (eff1 * k) 
      k <- rnorm(1, 1, z/7000)
      z <- z * (eff2 * k) 
    }
    lst1 <- c(lst1, y)
    lst2 <- c(lst2, z)
  }
  #print(min(lst1))
  #print(min(lst2))
  plot(lst1, lst2)
  YZ<-data.frame(y=lst1,z=lst2)
  sort_yz<-YZ[ order(YZ[1]),]
  return(sort_yz)
}

get_linear_model_2 <- function(p1,p2)
{
  lst1=get_Y_2(p1)
  lst2=get_Y_2(p2)
  YZ<-data.frame(y=lst1,z=lst2)
  sort_yz<-YZ[ order(YZ[1]),]
  get_linear_model(sort_yz)
}

get_linear_model <- function(lst)
{
  jpeg('rplot.jpg')
  plot(lst[[1]],lst[[2]], xlab ="Number of amplicon 1", ylab="Number of amplicon 2", main = "Lineal model of Simulated data")
  fit<-lm(lst[[2]]~lst[[1]])
  #Prediction intervals
  pred.int<-predict(fit, interval="prediction")
  abline(fit)
  print(lst)
  #print(pred.int)
  #matlines(sort(lst[[2]]), pred[,c("lwr","upr")],  lty=1, type="l")
  #Confidence intervals
  #conf.int =  predict(fit,interval="confidence")
  fitted.values <- pred.int[,1]
  pred.lower <- pred.int[,2]
  pred.upper <- pred.int[,3]
  lines(lst[[1]],fitted.values[],lwd=2, col="red")
  lines(lst[[1]],pred.lower[],lwd=2, col="blue")
  lines(lst[[1]],pred.upper[],lwd=2, col="blue")

  dev.off()
  return(fit) 
} 

plot_log<-function(YZ){
  lz<-log(YZ$z)
  ly<-log(YZ$y)
  fit.log=lm(lz~ly)
  print(fit.log)
  plot(ly,lz)
  abline(fit.log)
}