rob_kurt_skew <- function(x) {
  q = as.numeric(quantile(x, c(0.01,0.25,0.50,0.75,0.99), na.rm = T))
  skew = (q[4]+q[2]-2*q[3])/(q[4]-q[2])
  kurt = (q[5]-q[1])/(q[4]-q[2])
  
  return(c(skew,kurt))
}


# Summarise fft on each electrode.
fft_summary <- function(y) {
  Nyq.Freq = 400/2
  FFTFreqs <- getFFTFreqs(Nyq.Freq, y)
  FFT <- fft(y)
  modFFT <- Mod(FFT)
  FFTdata <- cbind(FFTFreqs, modFFT)
  FFTdata = data.table(FFTdata[1:nrow(FFTdata)/2,])
  wave_1 = log(sum(FFTdata[FFTFreqs>=0 & FFTFreqs<4,modFFT]))
  wave_2 = log(sum(FFTdata[FFTFreqs>=4 & FFTFreqs<8,modFFT]))
  wave_3 = log(sum(FFTdata[FFTFreqs>=8 & FFTFreqs<12,modFFT]))
  wave_4 = log(sum(FFTdata[FFTFreqs>=12 & FFTFreqs<30,modFFT]))
  wave_5 = log(sum(FFTdata[FFTFreqs>=30 & FFTFreqs<70,modFFT]))
  wave_6 = log(sum(FFTdata[FFTFreqs>=70 & FFTFreqs<200,modFFT]))  
  return(list(wave_1,wave_2,wave_3,wave_4,wave_5,wave_6))
}

# Summarise coherence between waves. Produces 5*15*16 values.
coh_summary <- function(dt_N) {
  u=0
  #N = sum(seq(1,15))
  #N=120
  N=56
  wave_1 = rep(NA,N)
  wave_2 = rep(NA,N)
  wave_3 = rep(NA,N)
  wave_4 = rep(NA,N)
  wave_5 = rep(NA,N)
  wave_6 = rep(NA,N)
  #for ( i in seq(1,15) ) {
  #  for ( j in seq(i+1,16) ) {
  #for ( i in seq(1,14,1) ) {
  #  for ( j in seq(i+2,16,2) ) {
  pairs = list(list(1,5),list(9,13),list(10,14),list(11,15),list(12,16),list(3,7),list(2,6),list(4,8))
  for ( pair in pairs ) {
      i = pair[[1]][[1]]
      j = pair[[2]][[1]]
      u = u + 1
      #if ( u%%2==0 ) cat('\r',i,',',j,',',u)
      #c=data.table(coh(dt_N[,i,with=FALSE],dt_N[,j,with=FALSE], f=400, plot=F))
      #c=data.table(coh(dt[id==idi,V1],dt[id==idi,V2], f=400))
      c=data.table(my_coh(dt_N[,i,with=FALSE],dt_N[,j,with=FALSE], f=400))
      #c=data.table(coh(dt[,V1],dt[,V2], f=400))
      
      c[,X:=X*1000]
      wave_1[u] = log(sum(c[X>=0 & X<4,V2]))
      wave_2[u] = log(sum(c[X>=4 & X<8,V2]))
      wave_3[u] = log(sum(c[X>=8 & X<12,V2]))
      wave_4[u] = log(sum(c[X>=12 & X<30,V2]))
      wave_5[u] = log(sum(c[X>=30 & X<70,V2]))
      wave_6[u] = log(sum(c[X>=70 & X<200,V2]))
    }
  }
  avg_wave1 = mean(wave_1, na.rm=T)
  avg_wave2 = mean(wave_2, na.rm=T)
  avg_wave3 = mean(wave_3, na.rm=T)
  avg_wave4 = mean(wave_4, na.rm=T)
  avg_wave5 = mean(wave_5, na.rm=T)
  avg_wave6 = mean(wave_6, na.rm=T)
  sum_wave1 = sum(wave_1, na.rm=T)
  sum_wave2 = sum(wave_2, na.rm=T)
  sum_wave3 = sum(wave_3, na.rm=T)
  sum_wave4 = sum(wave_4, na.rm=T)
  sum_wave5 = sum(wave_5, na.rm=T)
  sum_wave6 = sum(wave_6, na.rm=T)  
  #return(c(wave_1,wave_2,wave_3,wave_4,wave_5,
  #         avg_wave1,avg_wave2,avg_wave3,avg_wave4,avg_wave5,
  #         sd_wave1,sd_wave2,sd_wave3,sd_wave4,sd_wave5))
  return(c(avg_wave1,avg_wave2,avg_wave3,avg_wave4,avg_wave5,avg_wave6,
           sum_wave1,sum_wave2,sum_wave3,sum_wave4,sum_wave5,sum_wave6))
}
 
my_coh <- function(wave1, wave2, f) {
  input1 <- inputw(wave = wave1, f = f)
  wave1 <- input1$w
  f <- input1$f
  rm(input1)
  wave2 <- inputw(wave = wave2, f = f)$w
  n1 <- nrow(wave1)
  n2 <- nrow(wave2)
  
  #Y <- my_specpgram(cbind(wave1, wave2), fast = FALSE, taper = FALSE, spans = c(3, 3), plot = FALSE)$coh
  Y <-  my_specpgram(cbind(wave1, wave2), fast = FALSE, taper = FALSE, spans = c(3, 3))
  
  X <- seq(0, f/2000, length.out = nrow(Y))
  return(cbind(X, Y))
}

my_specpgram <- function(x, spans = NULL, kernel = NULL, taper = 0.1, pad = 0, 
          fast = TRUE, detrend = TRUE, na.action = na.fail) {
  #series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  #cat(x)
  xfreq <- frequency(x)
  x <- as.matrix(x)
  N <- N0 <- nrow(x)
  nser <- ncol(x)
  #cat(nser)
  #if (!is.null(spans)) 
  #  kernel <- {
  #    if (is.tskernel(spans)) 
  #      spans
  #    else kernel("modified.daniell", spans%/%2)
  #  }
  kernel <- kernel("modified.daniell", spans%/%2)
  #if (!is.null(kernel) && !is.tskernel(kernel)) 
  #  stop("must specify 'spans' or a valid kernel")
  # detrend also removes the mean.
  if (detrend) {
    t <- 1L:N - (N + 1)/2
    sumt2 <- N * (N^2 - 1)/12
    for ( i in 1L:ncol(x) ) {
      x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2
    }
  }
  #else if (demean) {
  #  x <- sweep(x, 2, colMeans(x), check.margin = FALSE)
  #}
  x <- spec.taper(x, taper)
  #u2 <- (1 - (5/8) * taper * 2)
  #u4 <- (1 - (93/128) * taper * 2)
  #if (pad > 0) {
  #  x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
  #  N <- nrow(x)
  #}
  #NewN <- if (fast) 
  #  nextn(N)
  #else N
  #x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
  N <- nrow(x)
  Nspec <- floor(N/2)
  #freq <- seq.int(from = xfreq/N, by = xfreq/N, length.out = Nspec)
  xfft <- mvfft(x)
  pgram <- array(NA, dim = c(N, ncol(x), ncol(x)))
  for (i in 1L:ncol(x)) {
    for (j in 1L:ncol(x)) {
      pgram[, i, j] <- xfft[, i] * Conj(xfft[, j])/(N0 * xfreq)
      pgram[1, i, j] <- 0.5 * (pgram[2, i, j] + pgram[N, i, j])
    }
  }
  #if (!is.null(kernel)) {
  for (i in 1L:ncol(x) ) {
    for (j in 1L:ncol(x)) {
      pgram[, i, j] <- kernapply(pgram[, i, j], kernel, circular = TRUE)
      #df <- df.kernel(kernel)
      #bandwidth <- bandwidth.kernel(kernel)
    }
  }
  #}
  #else {
  #  df <- 2
  #  bandwidth <- sqrt(1/12)
  #}
  #df <- df/(u4/u2^2)
  #df <- df * (N0/N)
  #bandwidth <- bandwidth * xfreq/N
  pgram <- pgram[2:(Nspec + 1), , , drop = FALSE]
  spec <- matrix(NA, nrow = Nspec, ncol = nser)
  for (i in 1L:nser) spec[, i] <- Re(pgram[1L:Nspec, i, i])
  #if (nser == 1) {
  #  coh <- phase <- NULL
  #}
  #else {
    coh <- phase <- matrix(NA, nrow = Nspec, ncol = nser * (nser - 1)/2)
    for (i in 1L:(nser - 1)) {
      for (j in (i + 1):nser) {
        coh[, i + (j - 1) * (j - 2)/2] <- Mod(pgram[, i, j])^2/(spec[, i] * spec[, j])
        #phase[, i + (j - 1) * (j - 2)/2] <- Arg(pgram[, 
        #                                              i, j])
      }
    }
  #}
  #for (i in 1L:nser) spec[, i] <- spec[, i]/u2
  #spec <- drop(spec)
  #spg.out <- list(freq = freq, spec = spec, coh = coh, phase = phase, 
  #                kernel = kernel, df = df, bandwidth = bandwidth, n.used = N, 
  #                orig.n = N0, series = series, snames = colnames(x), method = ifelse(!is.null(kernel), 
  #                                                                                    "Smoothed Periodogram", "Raw Periodogram"), taper = taper, 
  #                pad = pad, detrend = detrend, demean = demean)
  #class(spg.out) <- "spec"
  return(coh)
}

# u=0
# for ( i in seq(1,14,1) ) {
#   for ( j in seq(i+2,16,2) ) {
#    cat(i,' ',j,'\n')
#      u = u + 1
#   }
# }
cacf_calculator <- function(ts) {
   ACF = acf(ts, lag.max = 2000, plot = F)
   #sum_ACF_100  = sum(ACF$acf[50:100])
   #sum_ACF_50   = sum(ACF$acf[2:50])
   #sum_ACF_150  = sum(ACF$acf[100:150])
   #sum_ACF_200  = sum(ACF$acf[150:200])
   #sum_ACF_200p = sum(ACF$acf[200:500])
   cacf = sum(ACF$acf[0:2000])
   return(cacf)
}


