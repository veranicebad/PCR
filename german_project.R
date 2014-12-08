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


#xyplot(y ~ x, data = d) + layer_(panel.lmpolyline(...,col.line = 134, degree = 20))



fit_realdata<-function()
{
  mydata <- read.table("c:/Users/Vera/Documents/научкa/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
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

plot_reaSl_data<-function()
{
  #jpeg('rplot1.jpg')
  mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
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
  mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
  mydata<-normalization(mydata)
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

get_deletions_from_real_data<-function(){
  mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
  mydata<-normalization(mydata)
  data_without_del<-remove_deletions_from_real_data()
  data <- mydata
  data<-data[,-which(names(data) %in% names(data_without_del))]
  return(data)
}

normalization<-function(data){
  cs<-colSums(data[,3:ncol(data)])
  newdata <-data[,3:ncol(data)]
  for(i in 1:ncol(newdata)){
      newdata[,i]<-newdata[,i]*100/cs[i]
  }
  data[,3:ncol(data)] <- newdata
  return(data)
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
  plot(0,0,col='white',xlim=c(0,2.5),ylim=c(0,150))
  run_20_ab<-c()
  run_19_ab<-c()
  run_18_ab<-c()
  run_17_ab<-c()
  run_16_ab<-c()
  run_15_ab<-c()
  run_14_ab<-c()
  for(i in 1:nrow(mydata)){
    amp_i_run_14 <- as.numeric(run_14[i,])
    amp_i_run_15 <- as.numeric(run_15[i,])
    amp_i_run_16 <- as.numeric(run_16[i,])
    amp_i_run_17 <- as.numeric(run_17[i,])
    amp_i_run_18 <- as.numeric(run_18[i,])
    amp_i_run_19 <- as.numeric(run_19[i,])
    amp_i_run_20 <- as.numeric(run_20[i,])
    ab_amp_i_run_14<-nsamplGibs_a_b(1, amp_i_run_14,1)
    ab_amp_i_run_15<-nsamplGibs_a_b(1, amp_i_run_15,1)
    ab_amp_i_run_16<-nsamplGibs_a_b(1, amp_i_run_16,1)
    ab_amp_i_run_17<-nsamplGibs_a_b(1, amp_i_run_17,1)
    ab_amp_i_run_18<-nsamplGibs_a_b(1, amp_i_run_18,1)
    ab_amp_i_run_19<-nsamplGibs_a_b(1, amp_i_run_19,1)
    ab_amp_i_run_20<-nsamplGibs_a_b(1, amp_i_run_20,1)
    run_20_ab<-c(run_20_ab, ab_amp_i_run_20)
    run_19_ab<-c(run_19_ab, ab_amp_i_run_19)
    run_18_ab<-c(run_18_ab, ab_amp_i_run_18)
    run_17_ab<-c(run_17_ab, ab_amp_i_run_17)
    run_16_ab<-c(run_16_ab, ab_amp_i_run_16)
    run_15_ab<-c(run_15_ab, ab_amp_i_run_15)
    run_14_ab<-c(run_14_ab, ab_amp_i_run_14)
    points(ab_amp_i_run_14[1],ab_amp_i_run_14[2], pch=16, col="red")
    points(ab_amp_i_run_15[1],ab_amp_i_run_15[2], pch=16, col="orange")
    points(ab_amp_i_run_16[1],ab_amp_i_run_16[2], pch=16, col="yellow")
    points(ab_amp_i_run_17[1],ab_amp_i_run_17[2], pch=16, col="green")
    points(ab_amp_i_run_18[1],ab_amp_i_run_18[2], pch=16, col="cyan")
    points(ab_amp_i_run_19[1], ab_amp_i_run_19[2],pch=16, col="blue")
    points(ab_amp_i_run_20[1],ab_amp_i_run_20[2], pch=16, col="violet")
  }
  run_20_19_18_17_16_15_14_ab<-c(run_20_ab, run_19_ab, run_18_ab, run_17_ab,
                                 run_16_ab, run_15_ab, run_14_ab)
  return(run_20_19_18_17_16_15_14_ab)
}

nsamplGibs_a_b<-function(n, data,N0){
  library(e1071) 
  mean_data<-mean(data)
  std_data<-sd(data)  
  skewness_data<-skewness(data)
  kurtosis_data<-kurtosis(data)
  a_min=0
  b_min=0
  a_max=2.5
  b_max=150
  a_ans=0
  b_ans=0
  best_r=1000
  
  for(i in 1:n){
    a0=runif(1,min=a_min, max=a_max)
    b0=runif(1,min=b_min, max=b_max)
    ab<-samplGibs_a_b(a0,b0,data,N0)
    data_sample=get_Y_2(ab[1],ab[2],N0)
    mean_sample<-mean(data_sample)    
    sd_sample<-sd(data_sample)
    skewness_sample<-skewness(data_sample)
    kurtosis_sample<-kurtosis(data_sample)
    r=abs(mean_data-mean_sample)+abs(std_data-sd_sample)+
      abs(skewness_data-skewness_sample)+abs(kurtosis_data-kurtosis_sample)++
      abs(as.numeric(quantile(data_sample,0.025))-as.numeric(quantile(data,0.025)))
    +abs(as.numeric(quantile(data_sample,0.975))-as.numeric(quantile(data,0.975)))
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

samplGibs_a_b<-function(a0,b0,data,N0)
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
  a_max=2.5
  b_max=150
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
        data_sample=get_Y_2(a,b,N0)
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
      r=abs(mean_data-mean_sample)+abs(std_data-sd_sample)+
        abs(skewness_data-skewness_sample)+abs(kurtosis_data-kurtosis_sample)+
        abs(as.numeric(quantile(data_sample,0.025))-as.numeric(quantile(data,0.025)))
      +abs(as.numeric(quantile(data_sample,0.975))-as.numeric(quantile(data,0.975)))
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

nsamplGibs_N0<-function(n, data,a,b){
  library(e1071) 
  mean_data<-mean(data)
  std_data<-sd(data)  
  skewness_data<-skewness(data)
  kurtosis_data<-kurtosis(data)
  N0_min=0
  N0_max=10
  N0_ans=0
  best_r=1000
  
  for(i in 1:n){
    N0_0=runif(1,min=N0_min, max=N0_max)
    N0<-samplGibs_N0(N0_0,data,a,b)
    data_sample=get_Y_2(a,b,N0)
    mean_sample<-mean(data_sample)    
    sd_sample<-sd(data_sample)
    skewness_sample<-skewness(data_sample)
    kurtosis_sample<-kurtosis(data_sample)
    r=abs(mean_data-mean_sample)+abs(std_data-sd_sample)+
      abs(skewness_data-skewness_sample)+abs(kurtosis_data-kurtosis_sample)+
      abs(as.numeric(quantile(data_sample,0.025))-as.numeric(quantile(data,0.025)))
    +abs(as.numeric(quantile(data_sample,0.975))-as.numeric(quantile(data,0.975)))
    if(r<=best_r)
    {
      best_r <- r
      N0_ans=N0
    }
  }
  ans<-N0_ans
}

samplGibs_N0<-function(N0_0,data,a,b)
{
  n=50
  m=50
  mean_data<-mean(data)
  std_data<-sd(data)
  skewness_data<-skewness(data)
  kurtosis_data<-kurtosis(data)
  best_r=1000
  r = 1000
  N0_min=0
  N0_max=10
  N0_ans=0
  N0=N0_0
  prev_N0<-N0
  for(j in 1:n){
    print(best_r)
    print(j)
    mean_sample<-0
    sd_sample<-0
    skewness_sample<-0
    kurtosis_sample<-0
    for(i in 1:m){
      data_sample=get_Y_2(a,b,N0)
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
    r=abs(mean_data-mean_sample)+abs(std_data-sd_sample)+
      abs(skewness_data-skewness_sample)+abs(kurtosis_data-kurtosis_sample)+
      abs(as.numeric(quantile(data_sample,0.025))-as.numeric(quantile(data,0.025)))
    +abs(as.numeric(quantile(data_sample,0.975))-as.numeric(quantile(data,0.975)))
    if(r<=best_r)
    {
      best_r <- r
      N0_ans=N0
    }
    k=r_pr/(r+r_pr)
    t<-runif(1,0,1)
    if(t<=k)
    {
      N0 = prev_N0
    }
    prev_N0 = N0
    N0=runif(1,min=N0_min, max=N0_max)
    
  }
  ans<-N0_ans
}

get_N0<-function(a,b,N)
{
  N0_min<-1
  N0_max=100
  rmin<-1000
  N0<-N0_min
  for(i in N0_min:N0_max){
    val<-get_Y_2(a,b,i)
    mean_val<-mean(val)
    r<-abs(N-mean_val)
    if(r<rmin){
      N0<-i
      rmin<-r
    }
    print(rmin)
  }
  return(N0)
}

get_N0_for_deletion<-function(){
  mydata<-remove_deletions_from_real_data()
  run_20_pat_016<-mydata[, grep("IonXpress_20_016", colnames(mydata))]
  #run_20_pat_012<-mydata[, grep("IonXpress_20_012", colnames(mydata))]
  #run_19_pat_037<-mydata[, grep("IonXpress_19_037", colnames(mydata))]
  #run_18_pat_001<-mydata[, grep("IonXpress_18_001", colnames(mydata))]
  #run_18_pat_046<-mydata[, grep("IonXpress_18_046", colnames(mydata))]
  #run_18_pat_044<-mydata[, grep("IonXpress_18_044", colnames(mydata))]
  #run_18_pat_042<-mydata[, grep("IonXpress_18_042", colnames(mydata))]
  #run_17_pat_045<-mydata[, grep("IonXpress_17_045", colnames(mydata))]
  #run_17_pat_021<-mydata[, grep("IonXpress_17_021", colnames(mydata))]
  #run_17_pat_013<-mydata[, grep("IonXpress_17_013", colnames(mydata))]
  N_run_20_pat_016_AMPL1316862546 <- as.numeric(run_20_pat_016[grep("AMPL1316862546", row.names(mydata)),])
  N_run_20_pat_016_AMPL655136916 <- as.numeric(run_20_pat_016[grep("AMPL655136916", row.names(mydata)),])
  N_run_20_pat_016_AMPL478031510 <- as.numeric(run_20_pat_016[grep("AMPL478031510", row.names(mydata)),])
  N_run_20_pat_016_AMPL468281303 <- as.numeric(run_20_pat_016[grep("AMPL468281303", row.names(mydata)),])
  N_run_20_pat_016_AMPL612960426 <- as.numeric(run_20_pat_016[grep("AMPL612960426", row.names(mydata)),])
  N_run_20_pat_016_AMPL612959905 <- as.numeric(run_20_pat_016[grep("AMPL612959905", row.names(mydata)),])
  
  run_20<-mydata[, grep("IonXpress_20", colnames(mydata))]
  
  run_20_AMPL1316862546 <- as.numeric(run_20[grep("AMPL1316862546", row.names(mydata)),])
  ab_run_20_AMPL1316862546<-nsamplGibs_a_b(1, run_20_AMPL1316862546,1)
  #run_20_AMPL655136916 <- as.numeric(run_20[grep("AMPL655136916", row.names(mydata)),])
  #ab_run_20_AMPL655136916<-nsamplGibs_a_b(1, run_20_AMPL655136916,1) 
  #run_20_AMPL478031510 <- as.numeric(run_20[grep("AMPL478031510", row.names(mydata)),])
  #ab_run_20_AMPL478031510<-nsamplGibs_a_b(1, run_20_AMPL478031510,1)
  #run_20_AMPL468281303 <- as.numeric(run_20[grep("AMPL468281303", row.names(mydata)),])
  #ab_run_20_AMPL468281303<-nsamplGibs_a_b(1, run_20_AMPL468281303,1)
  #run_20_AMPL612960426 <- as.numeric(run_20[grep("AMPL612960426", row.names(mydata)),])
  #ab_run_20_AMPL612960426<-nsamplGibs_a_b(1, run_20_AMPL612960426,1)
  #run_20_AMPL612959905 <- as.numeric(run_20[grep("AMPL612959905", row.names(mydata)),])
  #ab_run_20_AMPL612959905<-nsamplGibs_a_b(1, run_20_AMPL612959905,1)
  N0_run_20_pat_016_AMPL1316862546<-c()
  for(j in 1:100){
    N0_run_20_pat_016_AMPL1316862546<-c(N0_run_20_pat_016_AMPL1316862546,
                                        get_N0(ab_run_20_AMPL1316862546[1],
                                               ab_run_20_AMPL1316862546[2],
                                               N_run_20_pat_016_AMPL1316862546))
  }
  plot(density(N0_run_20_pat_016_AMPL1316862546))
 }


get_Y_2<-function(a,b,N0)
{
  #jpeg('rplot.jpg')
  #p<-0.95
  lst1 <- c()
  y0<-N0
  for (i in 1:50)
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