N0_run_20_pat_016_AMPL655136916 <- c()
N0_run_20_pat_016_AMPL478031510 <- c()
N0_run_20_pat_016_AMPL468281303 <- c()
N0_run_20_pat_016_AMPL612960426 <- c()
N0_run_20_pat_016_AMPL612959905 <- c()
for(j in 1:100){
N0_run_20_pat_016_AMPL1316862546<-c(N0_run_20_pat_016_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_016_AMPL1316862546))
N0_run_20_pat_016_AMPL655136916<-c(N0_run_20_pat_016_AMPL655136916,
get_N0(ab_run_20_AMPL655136916[1],
ab_run_20_AMPL655136916[2],
N_run_20_pat_016_AMPL655136916))
N0_run_20_pat_016_AMPL478031510<-c(N0_run_20_pat_016_AMPL478031510,
get_N0(ab_run_20_AMPL478031510[1],
ab_run_20_AMPL478031510[2],
N_run_20_pat_016_AMPL478031510))
N0_run_20_pat_016_AMPL468281303<-c(N0_run_20_pat_016_AMPL468281303,
get_N0(ab_run_20_AMPL468281303[1],
ab_run_20_AMPL468281303[2],
N_run_20_pat_016_AMPL468281303))
N0_run_20_pat_016_AMPL612960426<-c(N0_run_20_pat_016_AMPL612960426,
get_N0(ab_run_20_AMPL612960426[1],
ab_run_20_AMPL612960426[2],
N_run_20_pat_016_AMPL612960426))
N0_run_20_pat_016_AMPL612959905<-c(N0_run_20_pat_016_AMPL612959905,
get_N0(ab_run_20_AMPL612959905[1],
ab_run_20_AMPL612959905[2],
N_run_20_pat_016_AMPL612959905))
}
N_run_20_pat_016_AMPL1316862546
plot(density(N0_run_20_pat_016_AMPL1316862546))
plot(density(N0_run_20_pat_016_AMPL655136916))
plot(density(N0_run_20_pat_016_AMPL478031510))
plot(density(N0_run_20_pat_016_AMPL468281303))
plot(density(N0_run_20_pat_016_AMPL468281303))
plot(density(N0_run_20_pat_016_AMPL612960426))
plot(density(N0_run_20_pat_016_AMPL612959905))
mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
mydata<-normalization(mydata)
run_20_pat_016<-mydata[, c(2, grep("IonXpress_20_016", colnames(mydata)))]
run_20_pat_012<-mydata[, c(2,grep("IonXpress_20_012", colnames(mydata)))]
run_19_pat_037<-mydata[, c(2,grep("IonXpress_19_037", colnames(mydata)))]
run_18_pat_001<-mydata[, c(2,grep("IonXpress_18_001", colnames(mydata)))]
run_18_pat_046<-mydata[, c(2,grep("IonXpress_18_046", colnames(mydata)))]
run_18_pat_044<-mydata[, c(2,grep("IonXpress_18_044", colnames(mydata)))]
run_18_pat_042<-mydata[, c(2,grep("IonXpress_18_042", colnames(mydata)))]
run_17_pat_045<-mydata[, c(2,grep("IonXpress_17_045", colnames(mydata)))]
run_17_pat_021<-mydata[, c(2,grep("IonXpress_17_021", colnames(mydata)))]
run_17_pat_013<-mydata[, c(2,grep("IonXpress_17_013", colnames(mydata)))]
N_run_20_pat_016_AMPL1316862546 <- as.numeric(run_20_pat_016[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL655136916 <- as.numeric(run_20_pat_016[grep("AMPL655136916", mydata$Target),2])
N_run_20_pat_016_AMPL478031510 <- as.numeric(run_20_pat_016[grep("AMPL478031510", mydata$Target),2])
N_run_20_pat_016_AMPL468281303 <- as.numeric(run_20_pat_016[grep("AMPL468281303", mydata$Target),2])
N_run_20_pat_016_AMPL612960426 <- as.numeric(run_20_pat_016[grep("AMPL612960426", mydata$Target),2])
N_run_20_pat_016_AMPL612959905 <- as.numeric(run_20_pat_016[grep("AMPL612959905", mydata$Target),2])
mydata<-remove_deletions_from_real_data()
run_20<-mydata[, grep("IonXpress_20", colnames(mydata))]
run_20_AMPL1316862546 <- as.numeric(run_20[grep("AMPL1316862546", mydata$Target),])
ab_run_20_AMPL1316862546<-nsamplGibs_a_b(1, run_20_AMPL1316862546,2)
run_20_AMPL655136916 <- as.numeric(run_20[grep("AMPL655136916", mydata$Target),])
ab_run_20_AMPL655136916<-nsamplGibs_a_b(1, run_20_AMPL655136916,1)
run_20_AMPL478031510 <- as.numeric(run_20[grep("AMPL478031510", mydata$Target),])
ab_run_20_AMPL478031510<-nsamplGibs_a_b(1, run_20_AMPL478031510,1)
run_20_AMPL468281303 <- as.numeric(run_20[grep("AMPL468281303", mydata$Target),])
ab_run_20_AMPL468281303<-nsamplGibs_a_b(1, run_20_AMPL468281303,1)
run_20_AMPL612960426 <- as.numeric(run_20[grep("AMPL612960426", mydata$Target),])
ab_run_20_AMPL612960426<-nsamplGibs_a_b(1, run_20_AMPL612960426,1)
run_20_AMPL612959905 <- as.numeric(run_20[grep("AMPL612959905", mydata$Target),])
ab_run_20_AMPL612959905<-nsamplGibs_a_b(1, run_20_AMPL612959905,1)
N0_run_20_pat_016_AMPL1316862546 <- c()
N0_run_20_pat_016_AMPL655136916 <- c()
N0_run_20_pat_016_AMPL478031510 <- c()
N0_run_20_pat_016_AMPL468281303 <- c()
N0_run_20_pat_016_AMPL612960426 <- c()
N0_run_20_pat_016_AMPL612959905 <- c()
for(j in 1:100){
N0_run_20_pat_016_AMPL1316862546<-c(N0_run_20_pat_016_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_016_AMPL1316862546))
N0_run_20_pat_016_AMPL655136916<-c(N0_run_20_pat_016_AMPL655136916,
get_N0(ab_run_20_AMPL655136916[1],
ab_run_20_AMPL655136916[2],
N_run_20_pat_016_AMPL655136916))
N0_run_20_pat_016_AMPL478031510<-c(N0_run_20_pat_016_AMPL478031510,
get_N0(ab_run_20_AMPL478031510[1],
ab_run_20_AMPL478031510[2],
N_run_20_pat_016_AMPL478031510))
N0_run_20_pat_016_AMPL468281303<-c(N0_run_20_pat_016_AMPL468281303,
get_N0(ab_run_20_AMPL468281303[1],
ab_run_20_AMPL468281303[2],
N_run_20_pat_016_AMPL468281303))
N0_run_20_pat_016_AMPL612960426<-c(N0_run_20_pat_016_AMPL612960426,
get_N0(ab_run_20_AMPL612960426[1],
ab_run_20_AMPL612960426[2],
N_run_20_pat_016_AMPL612960426))
N0_run_20_pat_016_AMPL612959905<-c(N0_run_20_pat_016_AMPL612959905,
get_N0(ab_run_20_AMPL612959905[1],
ab_run_20_AMPL612959905[2],
N_run_20_pat_016_AMPL612959905))
}
plot(density(N0_run_20_pat_016_AMPL1316862546))
plot(density(N0_run_20_pat_045_AMPL1316862546))
mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
mydata<-normalization(mydata)
run_20_pat_045<-mydata[, c(2, grep("IonXpress_20_045", colnames(mydata)))]
run_20_pat_016<-mydata[, c(2, grep("IonXpress_20_016", colnames(mydata)))]
run_20_pat_012<-mydata[, c(2,grep("IonXpress_20_012", colnames(mydata)))]
run_19_pat_037<-mydata[, c(2,grep("IonXpress_19_037", colnames(mydata)))]
run_18_pat_001<-mydata[, c(2,grep("IonXpress_18_001", colnames(mydata)))]
run_18_pat_046<-mydata[, c(2,grep("IonXpress_18_046", colnames(mydata)))]
run_18_pat_044<-mydata[, c(2,grep("IonXpress_18_044", colnames(mydata)))]
run_18_pat_042<-mydata[, c(2,grep("IonXpress_18_042", colnames(mydata)))]
run_17_pat_045<-mydata[, c(2,grep("IonXpress_17_045", colnames(mydata)))]
run_17_pat_021<-mydata[, c(2,grep("IonXpress_17_021", colnames(mydata)))]
run_17_pat_013<-mydata[, c(2,grep("IonXpress_17_013", colnames(mydata)))]
N_run_20_pat_045_AMPL1316862546 <- as.numeric(run_20_pat_045[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL1316862546 <- as.numeric(run_20_pat_016[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL655136916 <- as.numeric(run_20_pat_016[grep("AMPL655136916", mydata$Target),2])
N_run_20_pat_016_AMPL478031510 <- as.numeric(run_20_pat_016[grep("AMPL478031510", mydata$Target),2])
N_run_20_pat_016_AMPL468281303 <- as.numeric(run_20_pat_016[grep("AMPL468281303", mydata$Target),2])
N_run_20_pat_016_AMPL612960426 <- as.numeric(run_20_pat_016[grep("AMPL612960426", mydata$Target),2])
N_run_20_pat_016_AMPL612959905 <- as.numeric(run_20_pat_016[grep("AMPL612959905", mydata$Target),2])
mydata<-remove_deletions_from_real_data()
run_20<-mydata[, grep("IonXpress_20", colnames(mydata))]
run_20_AMPL1316862546 <- as.numeric(run_20[grep("AMPL1316862546", mydata$Target),])
ab_run_20_AMPL1316862546<-nsamplGibs_a_b(1, run_20_AMPL1316862546,2)
run_20_AMPL655136916 <- as.numeric(run_20[grep("AMPL655136916", mydata$Target),])
ab_run_20_AMPL655136916<-nsamplGibs_a_b(1, run_20_AMPL655136916,1)
run_20_AMPL478031510 <- as.numeric(run_20[grep("AMPL478031510", mydata$Target),])
ab_run_20_AMPL478031510<-nsamplGibs_a_b(1, run_20_AMPL478031510,1)
run_20_AMPL468281303 <- as.numeric(run_20[grep("AMPL468281303", mydata$Target),])
ab_run_20_AMPL468281303<-nsamplGibs_a_b(1, run_20_AMPL468281303,1)
run_20_AMPL612960426 <- as.numeric(run_20[grep("AMPL612960426", mydata$Target),])
ab_run_20_AMPL612960426<-nsamplGibs_a_b(1, run_20_AMPL612960426,1)
run_20_AMPL612959905 <- as.numeric(run_20[grep("AMPL612959905", mydata$Target),])
ab_run_20_AMPL612959905<-nsamplGibs_a_b(1, run_20_AMPL612959905,1)
N0_run_20_pat_045_AMPL1316862546 <- c()
N0_run_20_pat_016_AMPL1316862546 <- c()
N0_run_20_pat_016_AMPL655136916 <- c()
N0_run_20_pat_016_AMPL478031510 <- c()
N0_run_20_pat_016_AMPL468281303 <- c()
N0_run_20_pat_016_AMPL612960426 <- c()
N0_run_20_pat_016_AMPL612959905 <- c()
for(j in 1:100){
N0_run_20_pat_045_AMPL1316862546<-c(N0_run_20_pat_045_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_045_AMPL1316862546))
N0_run_20_pat_016_AMPL1316862546<-c(N0_run_20_pat_016_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_016_AMPL1316862546))
N0_run_20_pat_016_AMPL655136916<-c(N0_run_20_pat_016_AMPL655136916,
get_N0(ab_run_20_AMPL655136916[1],
ab_run_20_AMPL655136916[2],
N_run_20_pat_016_AMPL655136916))
N0_run_20_pat_016_AMPL478031510<-c(N0_run_20_pat_016_AMPL478031510,
get_N0(ab_run_20_AMPL478031510[1],
ab_run_20_AMPL478031510[2],
N_run_20_pat_016_AMPL478031510))
N0_run_20_pat_016_AMPL468281303<-c(N0_run_20_pat_016_AMPL468281303,
get_N0(ab_run_20_AMPL468281303[1],
ab_run_20_AMPL468281303[2],
N_run_20_pat_016_AMPL468281303))
N0_run_20_pat_016_AMPL612960426<-c(N0_run_20_pat_016_AMPL612960426,
get_N0(ab_run_20_AMPL612960426[1],
ab_run_20_AMPL612960426[2],
N_run_20_pat_016_AMPL612960426))
N0_run_20_pat_016_AMPL612959905<-c(N0_run_20_pat_016_AMPL612959905,
get_N0(ab_run_20_AMPL612959905[1],
ab_run_20_AMPL612959905[2],
N_run_20_pat_016_AMPL612959905))
}
plot(density(N0_run_20_pat_045_AMPL1316862546))
plot(density(N0_run_20_pat_016_AMPL1316862546))
plot(density(N0_run_20_pat_016_AMPL655136916))
plot(density(N0_run_20_pat_016_AMPL478031510))
plot(density(N0_run_20_pat_016_AMPL468281303))
plot(density(N0_run_20_pat_016_AMPL612960426))
mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
mydata<-normalization(mydata)
run_20_pat_045<-mydata[, c(2, grep("IonXpress_20_045", colnames(mydata)))]
run_20_pat_016<-mydata[, c(2, grep("IonXpress_20_016", colnames(mydata)))]
run_20_pat_012<-mydata[, c(2,grep("IonXpress_20_012", colnames(mydata)))]
run_19_pat_037<-mydata[, c(2,grep("IonXpress_19_037", colnames(mydata)))]
run_18_pat_001<-mydata[, c(2,grep("IonXpress_18_001", colnames(mydata)))]
run_18_pat_046<-mydata[, c(2,grep("IonXpress_18_046", colnames(mydata)))]
run_18_pat_044<-mydata[, c(2,grep("IonXpress_18_044", colnames(mydata)))]
run_18_pat_042<-mydata[, c(2,grep("IonXpress_18_042", colnames(mydata)))]
run_17_pat_045<-mydata[, c(2,grep("IonXpress_17_045", colnames(mydata)))]
run_17_pat_021<-mydata[, c(2,grep("IonXpress_17_021", colnames(mydata)))]
run_17_pat_013<-mydata[, c(2,grep("IonXpress_17_013", colnames(mydata)))]
N_run_20_pat_045_AMPL1316862546 <- as.numeric(run_20_pat_045[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL1316862546 <- as.numeric(run_20_pat_016[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL655136916 <- as.numeric(run_20_pat_016[grep("AMPL655136916", mydata$Target),2])
N_run_20_pat_016_AMPL478031510 <- as.numeric(run_20_pat_016[grep("AMPL478031510", mydata$Target),2])
N_run_20_pat_016_AMPL468281303 <- as.numeric(run_20_pat_016[grep("AMPL468281303", mydata$Target),2])
N_run_20_pat_016_AMPL612960426 <- as.numeric(run_20_pat_016[grep("AMPL612960426", mydata$Target),2])
N_run_20_pat_016_AMPL612959905 <- as.numeric(run_20_pat_016[grep("AMPL612959905", mydata$Target),2])
mydata<-remove_deletions_from_real_data()
run_20<-mydata[, grep("IonXpress_20", colnames(mydata))]
run_20_AMPL1316862546 <- as.numeric(run_20[grep("AMPL1316862546", mydata$Target),])
ab_run_20_AMPL1316862546<-nsamplGibs_a_b(1, run_20_AMPL1316862546,2)
run_20_AMPL655136916 <- as.numeric(run_20[grep("AMPL655136916", mydata$Target),])
ab_run_20_AMPL655136916<-nsamplGibs_a_b(1, run_20_AMPL655136916,1)
run_20_AMPL478031510 <- as.numeric(run_20[grep("AMPL478031510", mydata$Target),])
ab_run_20_AMPL478031510<-nsamplGibs_a_b(1, run_20_AMPL478031510,1)
run_20_AMPL468281303 <- as.numeric(run_20[grep("AMPL468281303", mydata$Target),])
ab_run_20_AMPL468281303<-nsamplGibs_a_b(1, run_20_AMPL468281303,1)
run_20_AMPL612960426 <- as.numeric(run_20[grep("AMPL612960426", mydata$Target),])
ab_run_20_AMPL612960426<-nsamplGibs_a_b(1, run_20_AMPL612960426,1)
run_20_AMPL612959905 <- as.numeric(run_20[grep("AMPL612959905", mydata$Target),])
ab_run_20_AMPL612959905<-nsamplGibs_a_b(1, run_20_AMPL612959905,1)
N0_run_20_pat_045_AMPL1316862546 <- c()
N0_run_20_pat_016_AMPL1316862546 <- c()
N0_run_20_pat_016_AMPL655136916 <- c()
N0_run_20_pat_016_AMPL478031510 <- c()
N0_run_20_pat_016_AMPL468281303 <- c()
N0_run_20_pat_016_AMPL612960426 <- c()
N0_run_20_pat_016_AMPL612959905 <- c()
for(j in 1:100){
N0_run_20_pat_045_AMPL1316862546<-c(N0_run_20_pat_045_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_045_AMPL1316862546))
N0_run_20_pat_016_AMPL1316862546<-c(N0_run_20_pat_016_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_016_AMPL1316862546))
N0_run_20_pat_016_AMPL655136916<-c(N0_run_20_pat_016_AMPL655136916,
get_N0(ab_run_20_AMPL655136916[1],
ab_run_20_AMPL655136916[2],
N_run_20_pat_016_AMPL655136916))
N0_run_20_pat_016_AMPL478031510<-c(N0_run_20_pat_016_AMPL478031510,
get_N0(ab_run_20_AMPL478031510[1],
ab_run_20_AMPL478031510[2],
N_run_20_pat_016_AMPL478031510))
N0_run_20_pat_016_AMPL468281303<-c(N0_run_20_pat_016_AMPL468281303,
get_N0(ab_run_20_AMPL468281303[1],
ab_run_20_AMPL468281303[2],
N_run_20_pat_016_AMPL468281303))
N0_run_20_pat_016_AMPL612960426<-c(N0_run_20_pat_016_AMPL612960426,
get_N0(ab_run_20_AMPL612960426[1],
ab_run_20_AMPL612960426[2],
N_run_20_pat_016_AMPL612960426))
N0_run_20_pat_016_AMPL612959905<-c(N0_run_20_pat_016_AMPL612959905,
get_N0(ab_run_20_AMPL612959905[1],
ab_run_20_AMPL612959905[2],
N_run_20_pat_016_AMPL612959905))
}
plot(density(N0_run_20_pat_016_AMPL612960426))
plot(density(N0_run_20_pat_016_AMPL612959905))
plot(density(N0_run_20_pat_016_AMPL612960426))
plot(density(N0_run_20_pat_016_AMPL468281303))
x<-c(1,2,3,4,6,5,7,8,7,8,9,11,11,15,20,22,30,40,44,50)
y<-c(1,2,3,4,6,5,7,8,7,8,10,11,11,15,20,22,30,43,44,54)
plot(lm(y~x))
lm(y~x)
fit<-lm(y~x)
fit
plot(fit)
abline(fit)
abline(fit)
plot(x,y)
abline(fit)
fit
plot(fit)
mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
mydata<-normalization(mydata)
run_20_pat_045<-mydata[, c(2, grep("IonXpress_20_045", colnames(mydata)))]
run_20_pat_016<-mydata[, c(2, grep("IonXpress_20_016", colnames(mydata)))]
run_20_pat_012<-mydata[, c(2,grep("IonXpress_20_012", colnames(mydata)))]
run_19_pat_037<-mydata[, c(2,grep("IonXpress_19_037", colnames(mydata)))]
run_18_pat_001<-mydata[, c(2,grep("IonXpress_18_001", colnames(mydata)))]
run_18_pat_046<-mydata[, c(2,grep("IonXpress_18_046", colnames(mydata)))]
run_18_pat_044<-mydata[, c(2,grep("IonXpress_18_044", colnames(mydata)))]
run_18_pat_042<-mydata[, c(2,grep("IonXpress_18_042", colnames(mydata)))]
run_17_pat_045<-mydata[, c(2,grep("IonXpress_17_045", colnames(mydata)))]
run_17_pat_021<-mydata[, c(2,grep("IonXpress_17_021", colnames(mydata)))]
run_17_pat_013<-mydata[, c(2,grep("IonXpress_17_013", colnames(mydata)))]
N_run_20_pat_045_AMPL1316862546 <- as.numeric(run_20_pat_045[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL1316862546 <- as.numeric(run_20_pat_016[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL655136916 <- as.numeric(run_20_pat_016[grep("AMPL655136916", mydata$Target),2])
N_run_20_pat_016_AMPL478031510 <- as.numeric(run_20_pat_016[grep("AMPL478031510", mydata$Target),2])
N_run_20_pat_016_AMPL468281303 <- as.numeric(run_20_pat_016[grep("AMPL468281303", mydata$Target),2])
N_run_20_pat_016_AMPL612960426 <- as.numeric(run_20_pat_016[grep("AMPL612960426", mydata$Target),2])
N_run_20_pat_016_AMPL612959905 <- as.numeric(run_20_pat_016[grep("AMPL612959905", mydata$Target),2])
mydata<-remove_deletions_from_real_data()
run_20<-mydata[, grep("IonXpress_20", colnames(mydata))]
run_20_AMPL1316862546 <- as.numeric(run_20[grep("AMPL1316862546", mydata$Target),])
ab_run_20_AMPL1316862546<-nsamplGibs_a_b(1, run_20_AMPL1316862546,2)
run_20_AMPL655136916 <- as.numeric(run_20[grep("AMPL655136916", mydata$Target),])
ab_run_20_AMPL655136916<-nsamplGibs_a_b(1, run_20_AMPL655136916,1)
run_20_AMPL478031510 <- as.numeric(run_20[grep("AMPL478031510", mydata$Target),])
ab_run_20_AMPL478031510<-nsamplGibs_a_b(1, run_20_AMPL478031510,1)
run_20_AMPL468281303 <- as.numeric(run_20[grep("AMPL468281303", mydata$Target),])
ab_run_20_AMPL468281303<-nsamplGibs_a_b(1, run_20_AMPL468281303,1)
run_20_AMPL612960426 <- as.numeric(run_20[grep("AMPL612960426", mydata$Target),])
ab_run_20_AMPL612960426<-nsamplGibs_a_b(1, run_20_AMPL612960426,1)
run_20_AMPL612959905 <- as.numeric(run_20[grep("AMPL612959905", mydata$Target),])
ab_run_20_AMPL612959905<-nsamplGibs_a_b(1, run_20_AMPL612959905,1)
N0_run_20_pat_045_AMPL1316862546 <- c()
N0_run_20_pat_016_AMPL1316862546 <- c()
N0_run_20_pat_016_AMPL655136916 <- c()
N0_run_20_pat_016_AMPL478031510 <- c()
N0_run_20_pat_016_AMPL468281303 <- c()
N0_run_20_pat_016_AMPL612960426 <- c()
N0_run_20_pat_016_AMPL612959905 <- c()
for(j in 1:100){
N0_run_20_pat_045_AMPL1316862546<-c(N0_run_20_pat_045_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_045_AMPL1316862546))
N0_run_20_pat_016_AMPL1316862546<-c(N0_run_20_pat_016_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_016_AMPL1316862546))
N0_run_20_pat_016_AMPL655136916<-c(N0_run_20_pat_016_AMPL655136916,
get_N0(ab_run_20_AMPL655136916[1],
ab_run_20_AMPL655136916[2],
N_run_20_pat_016_AMPL655136916))
N0_run_20_pat_016_AMPL478031510<-c(N0_run_20_pat_016_AMPL478031510,
get_N0(ab_run_20_AMPL478031510[1],
ab_run_20_AMPL478031510[2],
N_run_20_pat_016_AMPL478031510))
N0_run_20_pat_016_AMPL468281303<-c(N0_run_20_pat_016_AMPL468281303,
get_N0(ab_run_20_AMPL468281303[1],
ab_run_20_AMPL468281303[2],
N_run_20_pat_016_AMPL468281303))
N0_run_20_pat_016_AMPL612960426<-c(N0_run_20_pat_016_AMPL612960426,
get_N0(ab_run_20_AMPL612960426[1],
ab_run_20_AMPL612960426[2],
N_run_20_pat_016_AMPL612960426))
N0_run_20_pat_016_AMPL612959905<-c(N0_run_20_pat_016_AMPL612959905,
get_N0(ab_run_20_AMPL612959905[1],
ab_run_20_AMPL612959905[2],
N_run_20_pat_016_AMPL612959905))
}
N0_run_20_pat_016_AMPL1316862546
N0_run_20_pat_016_AMPL1316862546 <- as.data.frame(N0_run_20_pat_016_AMPL1316862546)
N0_run_20_pat_016_AMPL1316862546
N0_run_20_pat_016_AMPL1316862546$ind <- seq(1)
N0_run_20_pat_016_AMPL1316862546
colnames(N0_run_20_pat_016_AMPL1316862546) <- c("X", "ind")
head(N0_run_20_pat_016_AMPL1316862546)
library(ggplot2)
ggplot(aes(x = X), data = N0_run_20_pat_016_AMPL1316862546)
ggplot(aes(x = X), data = N0_run_20_pat_016_AMPL1316862546) + geom_shapes()
ggplot(aes(x = X), data = N0_run_20_pat_016_AMPL1316862546) + geom_point()
ggplot(aes(x = X), data = N0_run_20_pat_016_AMPL1316862546) + geom_density()
N0_run_20_pat_016_AMPL468281303 <- as.data.frame(N0_run_20_pat_016_AMPL468281303)
N0_run_20_pat_016_AMPL468281303$ind <- seq(1)
colnames(N0_run_20_pat_016_AMPL468281303) <- c("X", "ind")
N0_run_20_pat_016_AMPL468281303$ind <- seq(2)
totaldf <- rbind(N0_run_20_pat_016_AMPL1316862546, N0_run_20_pat_016_AMPL468281303)
totaldf
ggplot(totaldf, aes(x = X, group = ind)) + geom_density()
ggplot(aes(x = X), data = N0_run_20_pat_016_AMPL1316862546) + geom_density()
ggplot(totaldf, aes(x = X, group = ind)) + geom_density()
ggplot(aes(x = X), data = N0_run_20_pat_016_AMPL1316862546) + geom_density()
ggplot(totaldf, aes(x = X, group = ind)) + geom_density()
N0_run_20_pat_016_AMPL468281303$ind
N0_run_20_pat_016_AMPL468281303$ind <- seq(2,2)
N0_run_20_pat_016_AMPL468281303$ind
totaldf <- rbind(N0_run_20_pat_016_AMPL1316862546, N0_run_20_pat_016_AMPL468281303)
ggplot(aes(x = X), data = N0_run_20_pat_016_AMPL1316862546) + geom_density()
ggplot(totaldf, aes(x = X, group = ind)) + geom_density()
mydata <- read.table("c:/Users/Vera/Documents/научка/Run_22_hg19_v3.bcmatrix.xls", header=TRUE)
mydata<-normalization(mydata)
run_20_pat_043<-mydata[, c(2, grep("IonXpress_20_043", colnames(mydata)))]
run_20_pat_016<-mydata[, c(2, grep("IonXpress_20_016", colnames(mydata)))]
run_20_pat_012<-mydata[, c(2,grep("IonXpress_20_012", colnames(mydata)))]
run_19_pat_037<-mydata[, c(2,grep("IonXpress_19_037", colnames(mydata)))]
run_18_pat_001<-mydata[, c(2,grep("IonXpress_18_001", colnames(mydata)))]
run_18_pat_046<-mydata[, c(2,grep("IonXpress_18_046", colnames(mydata)))]
run_18_pat_044<-mydata[, c(2,grep("IonXpress_18_044", colnames(mydata)))]
run_18_pat_042<-mydata[, c(2,grep("IonXpress_18_042", colnames(mydata)))]
run_17_pat_045<-mydata[, c(2,grep("IonXpress_17_045", colnames(mydata)))]
run_17_pat_021<-mydata[, c(2,grep("IonXpress_17_021", colnames(mydata)))]
run_17_pat_013<-mydata[, c(2,grep("IonXpress_17_013", colnames(mydata)))]
N_run_20_pat_043_AMPL1316862546 <- as.numeric(run_20_pat_043[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL1316862546 <- as.numeric(run_20_pat_016[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_016_AMPL655136916 <- as.numeric(run_20_pat_016[grep("AMPL655136916", mydata$Target),2])
N_run_20_pat_016_AMPL478031510 <- as.numeric(run_20_pat_016[grep("AMPL478031510", mydata$Target),2])
N_run_20_pat_016_AMPL468281303 <- as.numeric(run_20_pat_016[grep("AMPL468281303", mydata$Target),2])
N_run_20_pat_016_AMPL612960426 <- as.numeric(run_20_pat_016[grep("AMPL612960426", mydata$Target),2])
N_run_20_pat_016_AMPL612959905 <- as.numeric(run_20_pat_016[grep("AMPL612959905", mydata$Target),2])
mydata<-remove_deletions_from_real_data()
run_20<-mydata[, grep("IonXpress_20", colnames(mydata))]
run_20_AMPL1316862546 <- as.numeric(run_20[grep("AMPL1316862546", mydata$Target),])
ab_run_20_AMPL1316862546<-nsamplGibs_a_b(1, run_20_AMPL1316862546,2)
N0_run_20_pat_016_AMPL1316862546
df1 <- as.data.frame(N0_run_20_pat_016_AMPL1316862546)
N0_run_20_pat_016_AMPL1316862546$ind <- seq(1,1)
colnames(N0_run_20_pat_016_AMPL1316862546) <- c("X", "ind")
totaldf <- rbind(totaldf, N0_run_20_pat_016_AMPL1316862546)
totaldf <- data.frame()
N0_run_20_pat_016_AMPL1316862546
df1 <- as.data.frame(N0_run_20_pat_016_AMPL1316862546)
N0_run_20_pat_016_AMPL1316862546$ind <- seq(1,1)
colnames(N0_run_20_pat_016_AMPL1316862546) <- c("X", "ind")
totaldf <- rbind(totaldf, N0_run_20_pat_016_AMPL1316862546)
N0_run_20_pat_016_AMPL1316862546
df1 <- as.data.frame(N0_run_20_pat_016_AMPL1316862546)
N0_run_20_pat_016_AMPL1316862546$ind <- seq(1,1)
colnames(N0_run_20_pat_016_AMPL1316862546) <- c("X", "ind")
totaldf <- rbind(totaldf, N0_run_20_pat_016_AMPL1316862546)
appenddf <- function(totaldf, dftoappend, index) {
dftoappend <- as.data.frame(dftoappend)
dftoappend$ind <- seq(index,index)
colnames(dftoappend) <- c("X", "ind")
totaldf <- rbind(totaldf, dftoappend)
}
totaldf <- data.frame()
colnames(totaldf) <- c("X", "ind")
totaldf <- data.frame(2)
colnames(totaldf) <- c("X", "ind")
totaldf <- data.frame(0,2)
colnames(totaldf) <- c("X", "ind")
totaldf <- appenddf(totaldf, N0_run_20_pat_016_AMPL1316862546, 1)
totaldf
totaldf <- appenddf(totaldf, N0_run_20_pat_016_AMPL468281303, 2)
totaldf <- appenddf(totaldf, N0_run_20_pat_016_AMPL478031510, 3)
totaldf <- appenddf(totaldf, N0_run_20_pat_016_AMPL612959905, 4)
totaldf <- appenddf(totaldf, N0_run_20_pat_016_AMPL612960426, 5)
totaldf <- appenddf(totaldf, N0_run_20_pat_016_AMPL655136916, 6)
totaldf <- appenddf(totaldf, N0_run_20_pat_045_AMPL1316862546, 7)
totaldf
as.factor(totaldf$ind)
labels(totaldf$ind)
totaldf$ind <- as.factor(totaldf$ind)
totaldf
labels(totaldf$ind)
totaldf$ind
levels(totaldf$ind)
levels(totaldf$ind) <- c("N0_run_20_pat_016_AMPL1316862546", "N0_run_20_pat_016_AMPL468281303", "N0_run_20_pat_016_AMPL478031510", "N0_run_20_pat_016_AMPL612959905", "N0_run_20_pat_016_AMPL612960426", "N0_run_20_pat_016_AMPL655136916", "N_run_20_pat_043_AMPL1316862546")
levels(totaldf$ind)
totaldf$ind
ggplot(aes(x = X, groups = ind), data = totaldf) + geom_density()
head(totaldf)
totaldf
ggplot(aes(x = X, groups = ind), data = totaldf) + geom_density()
ggplot(aes(x = X, group = ind), data = totaldf) + geom_density()
plot(density(N_run_20_pat_016_AMPL1316862546))
plot(density(N0_run_20_pat_016_AMPL1316862546))
plot(density(N0_run_20_pat_016_AMPL1316862546$X))
plot(density(N0_run_20_pat_016_AMPL468281303$X))
plot(density(N0_run_20_pat_016_AMPL478031510$X))
plot(density(N0_run_20_pat_016_AMPL478031510))
plot(density(N0_run_20_pat_016_AMPL612959905))
plot(density(N0_run_20_pat_016_AMPL612960426))
plot(density(N0_run_20_pat_016_AMPL655136916))
plot(density(N0_run_20_pat_045_AMPL1316862546))
ggplot(aes(x = X, group = ind), data = totaldf) + geom_density()
ggplot(aes(x = X, group = ind, col = ind), data = totaldf) + geom_density()
ggplot(aes(x = X, group = ind, col = ind, fill = ind), data = totaldf) + geom_density()
ggplot(aes(x = X, group = ind, col = c("blue", "grey"), fill = ind), data = totaldf) + geom_density()
ggplot(aes(x = X, group = ind, col = c("blue", "grey", "blue", "grey", "blue", "grey","grey"), fill = ind), data = totaldf) + geom_density()
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density(col = c("blue", "grey", "blue", "grey", "blue", "grey","grey"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density(col = ind)
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density(data = totaldfcol = ind)
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density(data = totaldf, col = ind)
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density()
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_hue(l=40)
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_brewer(palette="Set1")
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_brewer(palette="Set2")
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_manual(values=c("red", "blue", "green"))
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_brewer()
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_manual()
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_manual(c(, "#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_manual(c("#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind, col = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228"))
ggplot(aes(x = X, group = ind, col = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#FF131E", "#FF121E", "#FF111E", "#FF1G1E", "#FF1F1E", "#FF1E1E", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#FF131E", "#FF121E", "#FF111E", "#FF141E", "#FF1F1E", "#FF1E1E", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#FB0000", "#FB2800", "#FB5000", "#FB7800", "#FBA000", "#FBC800", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#FBA000", "#FB2800", "#FB5000", "#FB7800", "#FB0000", "#FBC800", "#41E228"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#FBA000",  "#FB7800","#FB5000", "#FB2800", "#FB0000", "#FBC800", "#41E228"))
run_20_pat_044<-mydata[, c(2, grep("IonXpress_20_044", colnames(mydata)))]
run_20_pat_045<-mydata[, c(2, grep("IonXpress_20_045", colnames(mydata)))]
N_run_20_pat_044_AMPL1316862546 <- as.numeric(run_20_pat_044[grep("AMPL1316862546", mydata$Target),2])
N_run_20_pat_045_AMPL1316862546 <- as.numeric(run_20_pat_045[grep("AMPL1316862546", mydata$Target),2])
N0_run_20_pat_044_AMPL1316862546 <- c()
N0_run_20_pat_045_AMPL1316862546 <- c()
for(j in 1:100){
N0_run_20_pat_044_AMPL1316862546<-c(N0_run_20_pat_044_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_044_AMPL1316862546))
N0_run_20_pat_045_AMPL1316862546<-c(N0_run_20_pat_045_AMPL1316862546,
get_N0(ab_run_20_AMPL1316862546[1],
ab_run_20_AMPL1316862546[2],
N_run_20_pat_045_AMPL1316862546))
}
labels(totaldf$ind) = c(1,2,3,4,5,6,7)
labels(totaldf$ind)
levels(totaldf$ind) <- c(1,2,3,4,5,6,7)
totaldf$ind <- as.numeric(totaldf$ind)
appenddf(totaldf, N0_run_20_pat_044_AMPL1316862546, 8)
totaldf <- appenddf(totaldf, N0_run_20_pat_044_AMPL1316862546, 8)
totaldf <- appenddf(totaldf, N0_run_20_pat_045_AMPL1316862546, 9)
totaldf$ind <- as.factor(totaldf$ind)
levels(totaldf$ind) <- c("N0_run_20_pat_016_AMPL1316862546", "N0_run_20_pat_016_AMPL468281303", "N0_run_20_pat_016_AMPL478031510", "N0_run_20_pat_016_AMPL612959905", "N0_run_20_pat_016_AMPL612960426", "N0_run_20_pat_016_AMPL655136916", "N_run_20_pat_043_AMPL1316862546", "N0_run_20_pat_044_AMPL1316862546", "N0_run_20_pat_045_AMPL1316862546")
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#FBA000",  "#FB7800","#FB5000", "#FB2800", "#FB0000", "#FBC800", "#41E228", "#41F128"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf) + geom_density() + scale_fill_manual(values = c("#FBA000",  "#FB7800","#FB5000", "#FB2800", "#FB0000", "#FBC800", "#41E228", "#41F128", "#41D328"))
ggplot(aes(x = X, group = ind, fill = ind), data = totaldf, ) + geom_density() + scale_fill_manual(values = c("#FBA000",  "#FB7800","#FB5000", "#FB2800", "#FB0000", "#FBC800", "#41E228", "#41F128", "#41D328")) + coord_cartesian(xlim = c(0, 3))
totaldf
savehistory("~/R/PCR/history.R")
