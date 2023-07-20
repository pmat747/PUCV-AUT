library ("stats")
library ("signal")
library ("seewave")
library ("psd")
library ("pracma")

########################################################################
data_generator = function (fs, duration, wvf, ampl=0, detector, filter=FALSE, setseed=0,
                            actPlot=TRUE, verbose=TRUE){
########################################################################
  # Input:
    # fs: sampling frequency
    # duration: duration (in second) of the output time serie
    # wvf: dataframe (time=time, hoft=h(t)) that contains the signal waveform sampled at fs
    # ampl: multiplication factor to change the source distance
    # detector: detector name
    # filter: name of the method
      #   "HP" : The fcut parameter is fixed internally (15 Hz)
      #   "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
      #   "AR" : AR model
      #   "prewhiten": use the R prewhiten function
    # setseed: if 0, random seed. Otherwise set to the value
    # The wvf is centered in the noise vector
  # 
  # Output: d$t: time vector
  #         d$x: noise+signal
  #         d$y: filtered (noise+signal)
    
  # check that duration is an integer > 0
  if ((duration%%1)>0){
    print("data_generator:duration must be a positive integer")
    return()
  }
  
  n=duration*fs
  wvf_size=length(wvf$hoft)
  
  if (n<wvf_size){
    print(sprintf("data_generator:the signal waveform duration is larger than %d", duration))
    return()
  }

  # The output vector will be 2 times larger than n
  factor=2
  data=noise_generator(factor,fs, duration, detector, setseed=setseed, filter=FALSE,
                       actPlot=FALSE, verbose)
  
  Y=data$x
  psd=data$psd         # 2 sided PSD
  n_data=length(Y)     # factor x n
  
  if (verbose==TRUE){
    print(sprintf("data_generator:size of the output: %d", n))
    print(sprintf("data_generator:size of the noise: %d", n_data))
    print(sprintf("data_generator:size of the signal: %d", wvf_size))
    print(sprintf("data_generator:amplitude of the signal: %f", ampl))
  }
  
  # Signal addition (centered at the middle of the data vector to avoid filtering leakage
  # at the beggining and end).
  ind1=floor((n_data-wvf_size)/2)
  
  for (i in 1:wvf_size){
    Y[ind1+i]=Y[ind1+i]+ampl*wvf$hoft[i]
  }
    
  # filter the time series if requested
  if (filter != FALSE){
    YY=filtering(Y, fs, filter, psd, verbose)
  }else{
    YY=Y
  }
  
  # generate a time series
  T = seq(1, n_data, by = 1)
  
  # select the original data size
  Tf = seq(1, n, by = 1)
    
  for (i in 1:n){
    Tf[i]=wvf$time[1]+i/fs
  }

  T_wvf=seq(1,wvf_size,by=1)
  for (i in 1:wvf_size){
    T_wvf[i]=wvf$time[1]+i/fs
  }
  
  Yf = seq(1, n, by = 1)
  YYf = seq(1, n, by = 1)
    
  for (i in 1:n){
    Yf[i]=Y[ind1+i]
    YYf[i]=YY[ind1+i]
  }
  
  if (actPlot==TRUE){
    if (filter == "HP" || filter == "spectrum" || filter == "prewhiten" || filter == "AR"){
      plot(T, Y, col="black", type="l", pch=1, panel.first = grid())
      points(T, YY, col="red", type="l", pch=2);              # (noise + signal) filtered 
      leg = c("noise+signal", "(noise+signal) filtered")
      col = c("black","red")
      legend (x=T[1]*1.1,y=max(Y)*.9,legend=leg,cex=.8,col=col,pch=c(1,2))
      
      plot(Tf, Yf, col="black", type="l", pch=1, panel.first = grid())
      points(Tf, YYf, col="red", type="l", pch=2);              # (noise + signal) filtered
      
      points(T_wvf,(wvf$hoft)*ampl,col="green",type="l",pch=3);  # signal only
        
      leg = c("noise", "(noise+signal) filtered", "signal only")
      col = c("black","red","green")
      legend (x=Tf[1]*1.1,y=max(Yf)*.9,legend=leg,cex=.8,col=col,pch=c(1,3))

      #freq2=fs*fftfreq(n)          # two-sided frequency vector
      #freq2[1]=0.0001
      ## Fourier transform
      #YYFT = sqrt(fs) * fft(YYf) / (n*sqrt(n)); # FFT computing and normalization
      #plot (freq2[1:int(n/2)], abs(YYFT)[1:int(n/2)], log="xy", type="l", xlab="frequency",
      #      ylab="ASD", col="grey",xlim=c(1, 4*10^3),panel.first = grid())
      #legend (x=10, y=min(abs(YYFT))*1.5, legend="signal+noise FT")
    }else{
      plot (Tf, Yf, type="l", col="black")
      legend (x=Tf[1]*1.1, y=max(Yf)*.9, legend="noise+signal")
    }
  }
    
  return(list(t=Tf,x=Yf,y=YYf))
  
}

########################################################################
signal_generator = function (fs, signal, pbOff=TRUE,
                             actPlot=TRUE, verbose=FALSE){
########################################################################
  # Assumptions : 
  #   all waveforms are in folder "inputs/New_simulations"
  #   simualted waveforms are sampled 16384 Hz
  #
  # Input:
  #   fs: sampling frequency of the output signal time series
  #   signal: name of the waveform
  #   pbOFF: if TRUE 100ms after the bounce will be discarded in the simulation
  # 
  # Output: 
  #   signal time series
  #   signal duration
  #   ratio for the g2 mode derived from the M and R values produced by the
  #   simulation code.
  
  # name: s11.2--LS220
  folder="inputs/New_simulations/"
  fs_orig=16384
  
  # Metadata  
  metadata_filename = paste(folder,"metadata.csv", sep="")
  meda = read.csv(metadata_filename, stringsAsFactors=FALSE)
  colnames(meda) = c("name","wvf_name","truedata_name", "tb")
  index=which(meda$name == signal)

  gw_filename=paste(folder,meda$wvf_name[index],sep="")
  truedata_filename=paste(folder,"Ratios/",meda$truedata_name[index],sep="")
  t_bounce=meda$tb[index]
  
  if (verbose==TRUE){
    print(gw_filename) 
    print(truedata_filename)
  }
  
  # Signal
  sXX = read.table(gw_filename); # V1 time, V2 signal
  colnames(sXX) = c ("time","hoft");
  duration=length(sXX$time)/fs_orig
 
  # True data to define the ratio
  # For g2 modes, ratio (x variable in TF19) is x = Mpns / Rpns^2 
  true_data = read.table(truedata_filename, sep = ",", comment.char = "#",header = TRUE);
  
  if (signal != "s20.0--SFHo"){
    true_data = cbind(true_data$time, true_data$mass_pns / true_data$r_pns^2);
  }

  colnames(true_data) = c ("time", "ratio");
  true_data = as.data.frame(true_data);
  
  if (verbose==TRUE){
    print(c("Number of samples at 16384 Hz:", length(sXX$time)))
    print(paste(signal, "duration=",duration, "s"))
    print(paste(signal,"tbounce=", t_bounce," s"))
  }
  
  # shift time such that t_bounce=0
  sXX$time=sXX$time - t_bounce  
  true_data$time=true_data$time - t_bounce 
  
  # signal sampled at 16384 Hz. Resampling at fs
  if (fs != fs_orig){
    resamp_factor=fs_orig/fs
    signaly=resample(sXX$hoft, 1./fs_orig, 1./fs)
    # decimate time vector
    signalx=seq(1,(fs*duration),by=1)
    for (i in 1:(fs*duration)) {
      signalx[i]=mean(sXX$time[((i-1)*resamp_factor+1):((i-1)*resamp_factor+resamp_factor)])
    }
  }else{
    resamp_factor=1
    signaly=sXX$hoft
    signalx=sXX$time}
    
  sYY = data.frame("time"=signalx,"hoft"=signaly)
  
  # remove times corresponding to the post-bounce period that is removed (100 ms) for all wvfs
  if (pbOff==TRUE){
    t_start=t_bounce+0.100
    true_data=subset(true_data, true_data$time >= 0.100)
    sYY=subset(sYY, sYY$time >= 0.100)  
  }
  
  if (actPlot == TRUE){
    plot(sXX$time, sXX$hoft, type="l", xlab="Time after bounce [s]", ylab="h(t)", xlim=c(0,duration),pch=1)
    points(sYY$time, sYY$hoft, col="red",pch=2)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
    title(signal)
    
    leg = c("wvf @16384 Hz", "resampled wvf (first 100ms after bounce removed)")
    col = c("black","red")
    legend (x=duration*.1,y=max(sXX$hoft)*.9,legend=leg,col=col,pch=c(1,2))
  }
  
  if (verbose == TRUE){
    print(c("Initial time of the output wvf", sYY$time[1]))
  }
  
  return (list(wvf=sYY, true_data=true_data, duration=duration))
  
}

########################################################################
noise_generator = function (factor,fs, duration, detector, setseed=0,
                            filter=FALSE, actPlot=TRUE, verbose=FALSE){
########################################################################
  # fs: sampling frequency
  # duration: duration (in second) of the output time serie
  # detector: name of the detector whose ASD will be used to generate colored noise
  # setseed: if 0, random seed. Otherwise set to the value.
  # filter method:
  #   "HP" : The fcut parameter is fixed internally (15 Hz)
  #   "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
  #   "AR" : AR model
  #   "prewhiten": use the R prewhiten function
  #
  # output: d$t: time vector
  #         d$x: noise
  ########################################################################  
  # check that duration is an integer > 0

  if ((duration%%1)>0){
    print("noise_generator:duration must be a positive integer")
    return()
  }
  
  if (duration < 10){
    # For 3G detectors we need to use a frequency resolution smaller than 0.1 Hz
    n=factor*duration*fs 
    }else{
      n=duration*fs
    }

  if (verbose==TRUE){
    print(sprintf("noise_generator:size of the noise output vector:%d", n))
  }
  
  # Noise generation
  freq2=fs*fftfreq(n)          # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:int(n/2)]        # one-sided frequency vector
  
  # Get the 2 sided PSD
  if (detector == "ALIGO"){
    psd=aLIGO_PSD_new(freq2, 2)
    }else{
    psd=PSD_fromfiles(freq2, 2, detector, verbose)
  }
    
  if (setseed >0){
    set.seed(setseed)
  }

  X = rnorm(n, mean=0, sd=1);           # Gaussian white noise
  XX = fft(X);                          # FFT computing
  XXX = XX*sqrt(psd)*sqrt(fs);          # Coloring
  Y = fft(XXX, inverse = TRUE);         # FFT inverse
  Y = Re(Y)/n;                          # noise in time domain
  # Note on the normalisation factor: 
  #  - n comes from the FFT and FFT inverse (sqrt(n) each)
  #  - to color properly the noise and keep the amplitude right 
  #    one needs to multiply by sqrt(psd) x sqrt(fs) 
  
  # filter the time series if requested
  if (filter != FALSE){
    YY=filtering(Y, fs, filter, psd, verbose)
  }else{
    YY=Y
  }
  
  if (verbose==TRUE){
    ss <- std(Y)
    print(sprintf("noise_generator:noise time serie sigma:%g", ss))
    ss <- std(YY)
    print(sprintf("noise_generator:filtered noise time serie sigma:%g", ss))
  }
  
  # generate a time series vector sampled at fs
  Tf = seq(1, n, by = 1)
  
  for (i in 1:n){
    Tf[i]=i/fs
  }
  
  if (actPlot==TRUE){
    # Time series
    T = seq(1, n, by = 1)
    ts.plot(Y); # noise only
    points(T, Y, col="black", type="l", pch=1, panel.first = grid())
    points(T, YY, col="red", type="l",pch=2)
    
    legend_str=c("simulated noise", "filtered noise")
    legend (x=0, y=abs(max(Y)), legend=legend_str, col=c("black","red"), pch=c(1,2))   
    
    # spectrum estimated
    psdest <- pspectrum(Y, Y.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, plot=FALSE,verbose=FALSE)
    psdest_filtered <- pspectrum(YY, YY.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, plot=FALSE,verbose=FALSE)
 
    # Fourier transform
    YFT = sqrt(2)*fft(Y)/sqrt(n);
    WFT = sqrt(2)*fft(YY)/sqrt(n);
    ymin=10^(ceiling(log10(min(abs(YFT)[1:int(n/2)])/sqrt(fs))))
    ymax=10^(ceiling(log10(max(abs(YFT)[1:int(n/2)])/sqrt(fs))))
    #ymin=1e-24
    #ymax=2e-21
    plot (freq1, abs(YFT)[1:int(n/2)]/sqrt(fs), log="xy", type="l", xlab="Frequency", ylab="ASD", 
          col="grey", xlim=c(1, fs/2), ylim=c(ymin,ymax), pch=1, panel.first = grid())

    lines(fs*psdest$freq, sqrt(psdest$spec)/sqrt(fs), col="blue", pch=2)
   
    lines(freq1, abs(WFT)[1:int(n/2)]/sqrt(fs), col="black", pch=4)        # factor 2 because FT is 2 sided
    lines(fs*psdest_filtered$freq[1:int(n/2)],                      # pspectrum is 1 sided
          sqrt(psdest_filtered$spec[1:int(n/2)])/sqrt(fs), col="green", pch=5)

    lines(freq1, sqrt(2*psd[1:int(n/2)]), col="red", pch=3)         # PSD is 2 sided PSD
    
    legend_str=c("col noise FT", "col noise spectrun", "ASD model", "filtered FT", "filtered spectrum")
    legend (x=100, y=min(abs(tail(YFT,-1)))*50000, legend=legend_str, col=c("grey","blue","red","black","green"), pch=c(1,2,3,4,5))   

    if (verbose==TRUE){
      s1 <- sqrt(2*trapz(fs*psdest$freq[1:int(n/2)], psdest$spec[1:int(n/2)]/fs))
      print(sprintf("noise_generator:colored noise rms:%g", s1))
    
      s2 <- sqrt(2*trapz(fs*psdest_filtered$freq[1:int(n/2)], psdest_filtered$spec[1:int(n/2)]/fs))
      print(sprintf("noise_generator:filtered noise rms:%g", s2))
    
      Sn_min=sqrt(2*min(psd))
      print(sprintf("minimal asd value:%g",Sn_min))
    }
    
  }
  return(list(t=Tf,x=Y,y=YY,psd=psd))
}

########################################################################
filtering = function(X, fs, method, psd=0, verbose=FALSE){
########################################################################
  # data processing of the input vector according to different methods
  # X: input data
  # fs: sampling frequency of X
  # method: filtering method
  #   "HP" : The fcut parameter is fixed internally (15 Hz)
  #   "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
  #   "AR" : AR model
  #   "prewhiten": use the R prewhiten function
  # psd: PSD required by the AR filering method

  # warning: the psd must be the 2 sided PSD. The size of the psd and data vectors must be equal
  if (length(X) != length(psd)){
    print(length(X))
    print(length(psd))
    warning("filtering:filtering::the data and psd vectors must have the same size. Please check")
  }
 
  n=length(X)
  duration=n/fs

  # compute noise sigma 
  freq2=fs*fftfreq(n)          # two-sided frequency vector

  s0 <- sqrt(trapz(freq2, psd))
  if (verbose==TRUE){
    print(sprintf("filtering:%s ASD noise rms: %g", detector, s0))
  }
  
  if (method == "HP"){
    fcut=10
    # filtfilt : zero phase filter (forward& backward)
    myfilter=butter(n=5, W=fcut/(fs/2), type="high")
    Y=filtfilt(filt=myfilter, x=X)} 
  else if (method == "AR"){
    if (psd==0){
      print("Filtering with AR method cannot be performed because noise psd has not been provided")
    }else{
      # generate another noise TS
      X1 = rnorm(n, mean=0, sd=1);            # Gaussian white noise
      XX1 = fft(X1);                          # FFT computing
      XXX1 = XX1*sqrt(psd);                   # Coloring
      Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
      Y1 = Re(Y1)*sqrt(fs)/n;                   # noise in time domain
    
      ar_model <- stats::ar(Y1,order.max=10, aic=FALSE ,method=c("yule-walker"), demean=TRUE);
      b <- stats::filter(x=X, filt=c(1, -ar_model$ar[1], -ar_model$ar[2], -ar_model$ar[3], 
                                   -ar_model$ar[4], -ar_model$ar[5], -ar_model$ar[6], 
                                   -ar_model$ar[7], -ar_model$ar[8], -ar_model$ar[9], 
                                   -ar_model$ar[10]), method="convolution", sides = 1);
      b[1]=b[2]=b[3]=b[4]=b[5]=b[6]=b[7]=b[8]=b[9]=b[10]=b[11]
      Y=b}
    }
  else if (method == "spectrum"){
    if (psd[1]==0){
      print("Filtering with specrum method cannot be performed because noise psd has not been provided")
    }else{
      # generate another noise TS
      X1 = rnorm(n, mean=0, sd=1);           # Gaussian white noise
      XX1 = fft(X1);                          # FFT computing
      XXX1 = XX1*sqrt(psd);                   # Coloring
      Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
      Y1 = Re(Y1);                   # noise in time domain
    
      # compute the PSD
      #myts <- ts(Y1, start=0, end=duration, frequency=fs)
      psdest <- pspectrum(Y1, Y1.frqsamp=fs, ntap.init=6, Nyquist.normalize = TRUE, 
                          plot=FALSE,verbose = FALSE)
      psdwhitening=rep(0, n);
      for(i in 1:(int(n/2))){
        psdwhitening[i]=psdest$spec[i]
        psdwhitening[n+1-i]=psdest$spec[i]
      }
    
      a = fft(X)                        # FFT computing and normalization
      b = a/sqrt(psdwhitening)          # whitening
      c = fft(b, inverse = TRUE);       # FFT inverse
      Y = s0*Re(c);                     # Normalisation factor of the 2 FFTs

      #    myfilter=butter(n=4,W=10/(fs/2),type="high")
      #    YY=filtfilt(filt=myfilter,x=YY)
    }
  }
  else if (method == "prewhiten"){
    # prewhiten
    myts <- ts(X, start=0, end=duration, frequency=fs)
    myts <- prewhiten(myts, AR.max=100, zero.pad="rear", plot=FALSE, verbose=FALSE)
    Y <- myts[['prew_ar']][1:n]}
  else{
    print("No filtering method specify")
    Y=X
  }
  
  return (Y)
}



########################################################################
fftfreq = function(n, d = 1){
########################################################################
  # surogate for the numpy fft.fftfreq function that generates the two sided 
  # frequency vector. Defaults d=1 means sampling frequency is 1. 
  # https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
  #
  # n: samples number
  # d: sample spacing (inverse of the sampling rate). Defaults to 1
  
  if(n%%2 == 0){# n is even
    out = c(seq(0, n/2-1, by = 1), seq(-n/2, -1, by=1)) / (d*n);
  }else{ # n is odd
    out = c(seq(0, (n-1)/2, by = 1), seq(-(n-1)/2, -1, by=1)) / (d*n);
  }
  
  return(out);
}

########################################################################
int = function(n){
########################################################################
  # https://stackoverflow.com/questions/31036098/what-is-the-difference-between-int-and-floor-in-python-3
  if(n < 0 ){
    return(-floor(abs(n)))
  }else{
    return(floor(n))
  }
}

########################################################################
aLIGO_PSD = function(f,type){
########################################################################
  # Original aLIGO PSD function used by Patricio
  cutoff = -109.35 + log(2e10);
  fn     = length(f);
  logpsd = rep(0, fn);
  if(f[1]==0){
    f[1]=f[2]
  }
  if(type == 1){
    for(i in 1:fn){
      x  = f[i]/215;
      x2 = x*x;
      logpsd[i] = log(1e-49) + log(x^(-4.14) -5/x2 + 111*(1-x2+0.5*x2*x2)/(1.+0.5*x2));
      
      if(logpsd[i]>cutoff){
        logpsd[i]=cutoff;
      }
    }
    output=exp(logpsd);
    
  }else{
    for(i in 1:(int(fn/2)+1)){
      x  = abs(f[i]/215);
      x2 = x*x;
      logpsd[i] = log(1e-49) + log(x^(-4.14) -5./x2 + 111.*(1-x2+0.5*x2*x2)/(1.+0.5*x2));
      
      if(logpsd[i]>cutoff){
        logpsd[i]=cutoff;
      }
      if(i>0){
        logpsd[fn-i]=logpsd[i];
      }
    }
    output=exp(logpsd)/2;            # Two sided PSD
  }
  return (output)
}

########################################################################
aLIGO_PSD_new = function(f,type){
########################################################################
  # aLIGO sensitivity curve: fit the data point from https://dcc.ligo.org/LIGO-T1800044/public
  # Type=1 --> one-sided PSD.
  # Type=2 --> two-sided PSD. 
  
  S1 = 5.0e-26;
  S2 = 1.0e-40;
  S3 = 1.4e-46;
  S4 = 2.7e-51;
  fcut = 10;
  cutoff = 1e-42;
  fn = length(f);
  output = rep(0, fn); #np.zeros(len(f))
  
  # to avoid issue with f=0
  if(f[1]==0){
    f[1]=f[2];
  }
  if(type == 1){
    for(i in 1:fn){
      x = abs(f[i]);
      output[i] =  S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }
    }
  }else{
    for(i in 1:(int(fn/2)+1)){ # range(int(len(f)/2)+1)
      x = abs(f[i]);
      output[i] = S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }

      # Wraparound frequency: f=0 is the first element (i=1), 
      # and all elements are symetric around index fn/2+1
      if(i>1 && i<fn/2+1){
        output[fn+2-i]=output[i]
      }
      
      
    }  
    output=output/2;            # Two sided PSD
#    output=shifter(output,-1)  # No more needed after correction 4 lines above
  }
  return(output);
}

shifter = function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

########################################################################
PSD_fromfiles=function(f, type, detector, verbose=FALSE){
########################################################################
  # Sensitivity curves for advnaced LIGO, advanced Virgo and KAGRA.
  # [Add refs here]
  # f: frequency vector
  # type=1 --> one-sided PSD.
  # type=2 --> two-sided PSD.
  # detector: name of the detector
  
  cutoff=1e-42            # For 2nd generator detectors

  # depending on the detector, the ASD is located in different columns of the data vector
  if (detector=="CE1"){
    psd_filename="PSD/curves_Jan_2020/ce1.txt"
    data=read.table(psd_filename);
    sens=data$V2        
    cutoff=1e-44}   
  
  if (detector=="CE2"){
    psd_filename="PSD/curves_Jan_2020/ce2.txt"
    data=read.table(psd_filename);
    sens=data$V2        
    cutoff=1e-44}   

  if (detector=="ET_B"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V2
    cutoff=1e-44}
  
  if (detector=="ET_C"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V2
    cutoff=1e-44}
  
  if (detector=="ET_D"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V4   # HF + LF
    cutoff=1e-44}
  
  if (detector=="ADV"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V7}   # Design
  
  if (detector=="aLIGO"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", toupper(detector))
    data=read.table(psd_filename);
    sens=data$V6}   # Design
  
  if (detector=="aLIGO2"){
    psd_filename="PSD/aLIGODesign.txt"
    data=read.table(psd_filename);
    sens=data$V2}   # Design  
  
  
  if (detector=="KAGRA"){
    psd_filename=sprintf("PSD/%s_sensitivity.txt", detector)
    data=read.table(psd_filename);
    sens=data$V6}
    
  if (exists("sens")==FALSE){
    stop(sprintf("Detector %s is not implemented in this code. You may want to use CE1, CE2, ET_B, ET_C, ET_D, aLIGO, ADV, KAGRA or ALIGO",detector))
  }
  
  n = length(f)
  fmin=f[1]
  if (type==1){
    fmax=f[n]
  } else{
    fmax=abs(f[n/2+1])}

  yl=sens[1]
  yr=sens[length(data$V1)]
  
  asd_func = approxfun(x = data$V1, y = sens, method = "linear",
                yleft=yl, yright=yr, rule = 1, f = 0, ties = "mean");
  
  if (type==1){
    asd = asd_func(f)
    psd = asd*asd
  }else{
    asd = rep(0, n);
    asd_1sided = asd_func(abs(f[1:int(n/2)]))
  
    if (length(asd_1sided) != int(n/2)){
      print (sprintf("Warning: ASD vector size %d is different from the frequency vector size %d", 
                   length(asd_1sided), n/2))
    }
    
    for(i in 1:(int(n/2))){
      asd[i]=asd_1sided[i];
      
      # Wraparound frequency: f=0 is the first element (i=1), 
      # and all elements are symetric around index n/2+1
      if(i>1){
        asd[n+2-i]=asd[i]
      }
    }  
    asd[n/2+1]=asd_func(abs(f[int(n/2)+1]))

    # Two sided psd
    asd=asd/sqrt(2);
    psd=asd*asd;
  }
  
  for (i in 1:n){
    if (psd[i]>cutoff){
      psd[i]=cutoff
    }
  }
  
  if (verbose==TRUE){
    fN=4096
    if (type==1){
      plot(f,psd,log="xy",col="blue",xlim=c(1, fN/2),pch=2)
      points(data$V1,sens*sens,col="red",type="l",pch=1)
    }else{
      plot(f,psd,log="y",col="blue",xlim=c(1, fN/2),pch=2)
      points(data$V1,0.5*sens*sens,col="red",type="l",pch=1)
    }
    leg = c(detector,"interpolated")
    col = c("red","blue")
    legend (x=500,y=psd[1]*0.8,legend=leg,cex=.8,col=col,pch=c(1,2))
  }
  
  return(psd)
}

########################################################################
compute_SNR = function(name, detector, fcut=0, dist=10, pbOff=FALSE){
########################################################################
  fs=16384
  signal=signal_generator(fs, name, pbOff=FALSE, actPlot=TRUE, verbose=FALSE)
  waveform=signal$wvf
  
  n=length(waveform$hoft)
  a = nextpow2(10*n)         #zero padding and rounding to the next power of 2
  n2=2^a

  # Remove or not 0.100s after the bounce (set hoft values to 0)
  if (pbOff == TRUE){
    ext=which(waveform$time<0.1)
    waveform$hoft[ext]=0
  }
  
  
  freq2 = fs*fftfreq(n2)         # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:int(n2/2)]       # one-sided frequency vector

  # Get the 1 sided PSD
  if (detector == "ALIGO"){
    psd=aLIGO_PSD_new(freq1, 1)
  }else{
    psd=PSD_fromfiles(freq1, 1, detector)
  }
  
  vec=rep(0,n2)
  for (i in 1:n){
    vec[n2/4+i]=vec[n2/4+i]+waveform$hoft[i]*10./dist
  }  
  
  hf=fft(vec)/sqrt(fs);                   # normalisation wrt the sampling

  hf=hf[1:(n2/2)]                # The integral is performed over positive freqs
  
  hf=subset(hf,freq1-fcut>0)
  psd=subset(psd,freq1-fcut>0)
  freq1=subset(freq1, freq1-fcut>0)

  integrand=abs(hf*Conj(hf))
  p=integrand/psd/fs
  
  snr=sqrt(4*trapz(freq1,p))
  
  #print(c(name,"SNR:",snr))
  
  plot (freq1, sqrt(freq1)*abs(hf), log="xy", type="l", xlab="Frequency", ylab="hchar", 
        col="grey", xlim=c(1, fs/2), ylim=c(1e-24,1e-20), pch=1, panel.first = grid())
  points(freq1,sqrt(psd), type="l", col="black",pch=2)
  leg = c("sqrt(f) x h~(f)", "ASD")
  col = c("grey","black")
  legend (x=1,y=6e-22,legend=leg,cex=.8,col=col,pch=c(1,2))
  title(c(name,"SNR:",snr))
  return(snr)  
  }

########################################################################
aLIGO_PSD_new2 = function(f,type){
########################################################################
# aLIGO sensitivity curve: fit the data point from https://dcc.ligo.org/LIGO-T1800044/public
# Type=1 --> one-sided PSD.
# Type=2 --> two-sided PSD. 
  
  S1 = 5.0e-26;
  S2 = 1.0e-40;
  S3 = 1.4e-46;
  S4 = 2.7e-51;
  fcut = 10;
  cutoff=1e-42;
  fn     = length(f);
  output = rep(0, fn); #np.zeros(len(f))
  
  # to avoid issue with f=0
  if(f[1]==0){
    f[1]=f[2];
  }
  if(type == 1){
    for(i in 1:fn){
      x = abs(f[i]);
      output[i] =  S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }
    }
  }else{
    for(i in 1:(int(fn/2)+1)){ # range(int(len(f)/2)+1)
      x = abs(f[i]);
      output[i] = S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }
      if(i>0){
        output[fn-i]=output[i]
      }
    }  
    output=output/2;            # Two sided PSD
    output=shifter(output,-1)   # Wraparound frequency: f=0 must be the last element
  }
  return(output);
}

