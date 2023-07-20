source("data_simulation.R")

N        = 2000; # Number of simulations x signal
M        = 2000; # Number of noise simulations

dist_s   = c(0.1, 5.05, 10); # Distances

duration = 2;    # Noise duration in seconds
fs       = 2^14; # Sampling frequency for signal and noise

#####################################
### GENERATE MULTIPLE SIMULATIONS ###
#####################################

# One signal simulation
signal_s = c("s11.2--LS220"); # signals to be used 
detector = c("aLIGO");

# We can also simulate multiple signals using the following
# signal_s=c("s20.0--SFHo", "s11.2--LS220", "s15.0--GShen", "s15.0--LS220", 
#          "s15.0--SFHo", "s20.0--LS220", "s25.0--LS220", "s40.0--LS220");

#detectors = c("aLIGO", "CE1", "CE2", "ET_B", "ET_C", "ET_D");
#detectors = c("CE1", "CE2", "ET_B", "ET_C", "ET_D", "aLIGO", "ADV", "KAGRA", "ALIGO");

filtering_method="prewhiten";

for(s in signal_s){
  
  signal_name = s;
  signal      = signal_generator(fs, signal_name, actPlot=FALSE, verbose=FALSE);
  wvf         = signal$wvf; # Pure signal
  duration    = signal$duration;
  
  for(dist in dist_s){
    
    for(i in 1:N){
      
      aux = data_generator(fs, duration, wvf, ampl=10/dist, detector, 
                           filter = filtering_method, 
                           setseed=0, actPlot=FALSE, verbose=FALSE);
      
      aux = data.frame("V1"=aux$t,"V2"=aux$y);
      
      name = paste("simulations/",signal_name,"_", dist, "kpc_sim", i, ".txt", sep = "");
      
      write.table(aux, name, row.names=FALSE, col.names = FALSE);
      
    }# end signal simulation
  }# end distances
}

########################
### Simulating noise ###
########################

signal_name = c("s20.0--SFHo"); # No relevant for generating noise
signal      = signal_generator(fs, signal_name, actPlot=FALSE, verbose=FALSE);
wvf         = signal$wvf; # Pure signal
duration    = signal$duration;

for(i in 1:M){
  
  aux = data_generator(fs, duration, wvf, ampl=0, detector, 
                       filter = filtering_method, 
                       setseed=0, actPlot=FALSE, verbose=FALSE);
  
  aux = data.frame("V1"=aux$t,"V2"=aux$y);
  
  name = paste("simulations/",detector,"_noise_", duration, "sec_sim", i, ".txt", sep = "");
  
  write.table(aux, file = name, row.names=FALSE, col.names = FALSE);
  
}



