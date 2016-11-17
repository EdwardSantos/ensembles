# Gets the frequencies returned by the FFT function
getFFTFreqs <- function(Nyq.Freq, data)
{
  if ((length(data) %% 2) == 1) # Odd number of samples
  {
    FFTFreqs <- c(seq(0, Nyq.Freq, length.out=(length(data)+1)/2), 
                  seq(-Nyq.Freq, 0, length.out=(length(data)-1)/2))
  }
  else # Even number
  {
    FFTFreqs <- c(seq(0, Nyq.Freq, length.out=length(data)/2), 
                  seq(-Nyq.Freq, 0, length.out=length(data)/2))
  }
  
  return (FFTFreqs)
}

# FFT plot
# Params:
# x,y -> the data for which we want to plot the FFT 
# samplingFreq -> the sampling frequency
# shadeNyq -> if true the region in [0;Nyquist frequency] will be shaded
# showPeriod -> if true the period will be shown on the top
# Returns a list with:
# freq -> the frequencies
# FFT -> the FFT values
# modFFT -> the modulus of the FFT
plotFFT <- function(x, y, samplingFreq, shadeNyq=TRUE, showPeriod = TRUE)
{
  Nyq.Freq <- samplingFreq/2
  FFTFreqs <- getFFTFreqs(Nyq.Freq, y)
  
  FFT <- fft(y)
  modFFT <- Mod(FFT)
  FFTdata <- cbind(FFTFreqs, modFFT)
  plot(FFTdata[1:nrow(FFTdata)/2,], t="l", pch=20, lwd=2, cex=0.8, main="",
       xlab="Frequency (Hz)", ylab="Power")
  if (showPeriod == TRUE)
  {
    # Period axis on top        
    a <- axis(3, lty=0, labels=FALSE)
    axis(3, cex.axis=0.6, labels=format(1/a, digits=2), at=a)
  }
  if (shadeNyq == TRUE)
  {
    # Gray out lower frequencies
    rect(0, 0, 2/max(x), max(FFTdata[,2])*2, col="gray", density=30)
  }
  
  ret <- list("freq"=FFTFreqs, "FFT"=FFT, "modFFT"=modFFT)
  return (ret)
}

