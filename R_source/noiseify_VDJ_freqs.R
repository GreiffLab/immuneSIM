
.noiseify_VDJ_freqs<-function(vdj_list,vdj_noise_level=0.2){
  #introduces noise into the VDJ frequencies.

  #from https://stackoverflow.com/questions/19343133/setting-upper-and-lower-limits-in-rnorm
  rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
    stats::qnorm(stats::runif(n, stats::pnorm(a, mean, sd), stats::pnorm(b, mean, sd)), mean, sd)
  }

  #create rtnorm distributions of needed length with a mean=0 and
  #sd=vdj_noise_level bounded by -1,1.
  V_noise<-rtnorm(length(vdj_list$V$frequency),0,vdj_noise_level,a=-1,b=1)
  D_noise<-rtnorm(length(vdj_list$D$frequency),0,vdj_noise_level,a=-1,b=1)
  J_noise<-rtnorm(length(vdj_list$J$frequency),0,vdj_noise_level,a=-1,b=1)

  #add noise to the frequencies and readjust to sum up to 1
  vdj_list$V$frequency<-vdj_list$V$frequency+V_noise*vdj_list$V$frequency
  vdj_list$V$frequency<-vdj_list$V$frequency/sum(vdj_list$V$frequency)

  vdj_list$D$frequency<-vdj_list$D$frequency+D_noise*vdj_list$D$frequency
  vdj_list$D$frequency<-vdj_list$D$frequency/sum(vdj_list$D$frequency)

  vdj_list$J$frequency<-vdj_list$J$frequency+J_noise*vdj_list$J$frequency
  vdj_list$J$frequency<-vdj_list$J$frequency/sum(vdj_list$J$frequency)


  return(vdj_list)

}

