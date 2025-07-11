# this script includes the functions and helper settings used within DBEN:

# when run from within scripts - relative paths for download from github:
# base.path <- "../"

models_avail = c("JULES-RED","FATES","ORCHIDEE","LPJ-GUESS","CABLE-POP","BiomeEP","BiomeE-Standalone","SEIB-DGVM","EDv3")

#knitr::opts_chunk$set(echo = TRUE)
library(ncdf4)
#devtools::install_github("MagicForrest/DGVMTools", ref = "master", dependencies = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)
library(DGVMTools)
library(dplyr)
library(maps)
library(maptools) # no available after R 4.2.0. Alternative: sf
library(zoo) # datetime  handling of observations
library(gridExtra) # for plotting multiple ggplot objects beside each other
library(psych) # for scatter-history plotting for mort and wbgrowth rates
library(here)
#source the format metadata for DBEN project to integrate with DGVMTools:
#source("Format-DBEN_paper1.R")
# conversion helpers:
year_to_seconds = 60*60*24*365
seconds_to_year = 1/(60*60*24*365)
create_figs = TRUE

# vector of all variables simulated for DBEN 
variables <-  c("cveg","cwood","cwood_size","nstem_size","lai","CA","BA","height","WBgrowth","BAgrowth","cmort","stemmort","gpp","npp","nbp")
# vector of all sites simulated for DBEN
sites = c("FIN","BIA","BCI")

basepath = "../model_outputs/"
# set all file directories here:
file.dir.julesred <- paste0(basepath,"JULES-RED/")
file.dir.fates    <- paste0(basepath,"FATES/")# we are at v5
file.dir.lpjguess <- paste0(basepath,"LPJGUESS/2_processed/")#long_runs/")
#file.dir.lpjguess <- "/Users/annemarie/OneDrive - Lund University/1_TreeMort_onedrive/2_Analysis/3_analysis_demographic_model_intercomparison/D-BEN-site-sims/Paper_1/Outputs/LPJGUESS/2_processed_paper0/"
file.dir.cablepop <- paste0(basepath,"CABLE-POP/")
#file.dir.julesred <- paste0(base.path,"Outputs/JULES-RED/v6/")
file.dir.biomeEP  <- paste0(basepath,"BiomeEP/")
file.dir.biomeE_Standalone <- paste0(basepath,"BiomeE-Standalone/")
file.dir.seibdgvm <- paste0(basepath,"SEIB-DGVM/")
file.dir.EDv3     <- paste0(basepath,"EDv3/")
file.dir.orchidee <- paste0(basepath,"ORCHIDEE/")

abcdfeghijkl <- function(model_name){
  idx <- which(models_avail==model_name)
  abcdfeghijkl <- c("a)", "b)" ,"c)","d)", "e)" ,"f)","g)", "h)" ,"i)","j)", "k)" ,"l)")
  return(abcdfeghijkl[idx])
}

get_output_LPJGUESS <- function(site,run,var,file.dir,co2_levels){
  model_name = "LPJGUESS"
  
  # to make backwards compatible:
  if(co2_levels == "412ppm" | co2_levels == "PS_412ppm/" ){
    co2_levels = "PS_412ppm/"
  }else if(co2_levels == "562ppm" | co2_levels == "PS_562ppm/" ){
    co2_levels = "PS_562ppm/"
  }else(# enable read-in for runs where the CO2 levels do not matter ( multi and single PFT runs for succession/turnover part of the paper)
    co2_levels <- co2_levels
  )
  
  file.dir = paste0(file.dir,co2_levels)
  #print(file.dir)
  
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste("LPJ-GUESS", site,'-', run),
                            forcing.data = "cru_jra2.2")
  
  var.to.plot <- getField_DBEN(source.in, 
                               quant = get_quantity(var),
                               file.name = paste0(file.dir,
                                                  model_name,"_",var,"_",run,"_",site,".nc"),
                               model_name = model_name)
  
  # change metadata
  var.to.plot@first.year <- 1 # first year of simulation
  var.to.plot@last.year  <- dim(var.to.plot@data)[1]
  # change real data:
  var.to.plot@data$Year  <- seq(1,dim(var.to.plot@data)[1])
  
  #[TODO] probably change unit later in actual model output and for dben standard..
  if(var== "cmort" |  var == "cmort_age" |
     var == "cmort_other" | var == "cmort_greff" | var == "cmort_thin"| var =="nbp") {
    #convert to -yr:
    # all units are in year now.
    #var.to.plot <- change_unit_stoy(var.to.plot) 
    # LPJG outputs in units of kgC m2 yr, so just changing the unit name here.
    var.to.plot@quant@units <- "kgC m-2 yr-1"
  }
  
  # if(var =="gpp" ){
  #  var.to.plot <- change_unit_stoy(var.to.plot)
  # }
  
  return(var.to.plot)
}

get_output_BiomeEP <- function(site,run,var,file.dir,co2_levels){
  model_name = "BiomeEP"
  #site = c("Fi1")
  #run = c("p0")
  # set metadata 
  
    if(co2_levels =="412ppm"){
      co2_levels ="aCO2"
    }
    if(co2_levels =="562ppm"){
      co2_levels ="eCO2"
    }
    cmort_var_quant ="cmort_size"
    stemmort_var_quant = "stemmort_size"
 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste(model_name, site,'-', run),
                            forcing.data = "cru_jra2.2")
  
  # address the fact that for filename with var = "cmort_size", the variable name within = "cmort"
  if(var== "cmort_size"){
    var.to.plot <- getField_DBEN(source.in, 
                                 quant = get_quantity(cmort_var_quant),
                                 file.name = paste0(file.dir,co2_levels,"/",site,"/",
                                                    model_name,"_",var,"_",run,"_",site,"_",co2_levels,".nc"),
                                 model_name = model_name)
  }else if(var == "cmort_pft"){ #BiomeEP has a cmort by pft. Making this importable here:
    var.to.plot <- getField_DBEN(source.in, 
                                 quant = get_quantity("cmort"),
                                 file.name = paste0(file.dir,co2_levels,"/",site,"/",
                                                    model_name,"_",var,"_",run,"_",site,"_",co2_levels,".nc"),
                                 model_name = model_name)
  }else if(var == "stemmort_size"){ # see above.
    var.to.plot <- getField_DBEN(source.in, 
                                 quant = get_quantity(stemmort_var_quant),
                                 file.name = paste0(file.dir,co2_levels,"/",site,"/",
                                                    model_name,"_",var,"_",run,"_",site,"_",co2_levels,".nc"),
                                 model_name = model_name)
  }else{ # normal procedure:
    var.to.plot <- getField_DBEN(source.in, 
                                 quant = get_quantity(var),
                                 file.name = paste0(file.dir,co2_levels,"/",site,"/",
                                                    model_name,"_",var,"_",run,"_",site,"_",co2_levels,".nc"),
                                 model_name = model_name)
  }
  
  #[TODO] probably change unit later..
  if(var == "cmort" | var == "cmort_age" |
     var == "cmort_other" | var == "cmort_greff") {
    # convert to -yr:
    #var.to.plot <- change_unit_stoy(var.to.plot)
  }
  
  if(var == "cmort_size" ){
    var.to.plot@quant@units <- "kgC m-2 yr-1"
  }
  
  if(var=="stemmort_size"|var =="stemmort"){ #Ensheng. 17.05.2023: The unit of "Stemmort" is n_stems per year per m2 (how many trees died in a year in one square meter of land. So it is very small.
    # convert to -ha:
    end <- dim(var.to.plot@data)[2]
    var.to.plot@data[,4:end] <- var.to.plot@data[,4:end] *10000 # m to ha
  }
  #  if(var =="npp"){
  #  var.to.plot <- change_npp_unit(var.to.plot)
  #}
  return(var.to.plot)
}

get_output_BiomeE_standalone <- function(site,run,var,file.dir,co2_levels){
  model_name = "BiomeE"
  if(co2_levels=="562ppm"){
    co2_levels ="eCO2"
  }
  if(co2_levels =="412ppm"){
    co2_levels ="aCO2"
  }
  #run = c("p0")
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir.biomeE_Standalone ,
                            format = DBEN,
                            name = paste(model_name, site,'-', run),
                            forcing.data = "cru_jra2.2")
  
  var.to.plot <- getField_DBEN(source.in, 
                               quant = get_quantity(var),
                               file.name = paste0(file.dir,model_name,"_PS_",site,"_",co2_levels,"_",run,"_",var,".nc") ,
                               model_name = model_name)
  
  #paste0(file.dir,model_name,"_",run,"_",site,"_",co2_levels,"_00_",var,".nc")
  # /Users/annemarie/Documents/1_TreeMort/2_Analysis/3_analysis_demographic_model_intercomparison/Outputs/BiomeE-Standalone/810yrsRun/NetCDF/BiomeE_PS_FIN_aCO2_00_AGB.nc
  #[TODO] probably change unit later..
  if(var == "cmort" | var == "cmort_age" |
     var == "cmort_other" | var == "cmort_greff") {
    #convert to -yr:
    var.to.plot <-change_unit_stoy(var.to.plot)
  }
  
  if(var=="stemmort_size"){ #Ensheng. 17.05.2023: The unit of "Stemmort" is n_stems per year per m2 (how many trees died in a year in one square meter of land. So it is very small.
    #convert to -ha:
    end <- dim(var.to.plot@data)[2]
    var.to.plot@data[,4:end] <- var.to.plot@data[,4:end] *10000 # m to ha
  }
  #  if(var =="npp"){
  #  var.to.plot <- change_npp_unit(var.to.plot)
  #}
  
  #update quantities:
  if(var =="WDgrowTot"){
    #WDgrow+ WDrepr  (new seedlings from reproduction)+ WDmgrt (migrated plants, once in year 31), just for traceability, that will be updated as WBgrowth, as that is what "Total wood growth" means for DBEN. 
    var.to.plot@quant = get_quantity("WBgrowth")
  }
  
  if(var =="WDmortTot"){
    #WDmort+WDkill (killed trees because of low density)+WDdstb (removed by disturbance, year 31)
    var.to.plot@quant = get_quantity("cmort")
  }
  
  return(var.to.plot)
}

get_output_CABLEPOP <- function(site,run,var,file.dir,co2_levels){
  model_name="CABLE-POP"
  # flexibly account for file naming:
  if(run == "0"){
    #benchmark run
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var,"_P0_",site,".nc")
  }
  else{ # sensitivity runs:
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var,"_PS_",site,"_",run,".nc")
  }
  
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste(model_name, site,'-', run),
                            forcing.data =  "cru_jra2.2")
  
  var.to.plot <- getField_DBEN(source.in, 
                               quant = get_quantity(var),
                               file.name = file.name,
                               model_name="CABLE-POP")
  
  # c("cmort_crowd", "cmort_dist", "cmort_res")
  #[TODO] Juergen unit seems to be off by magnitudes, but if I multiply with year_to_seconds multiplier, t is fine. 
  # has this maybe been "converted" to year twice, and I have to roll back one conversion here?
  if(var == "cmort"| var == "cmort_crowd"| var == "cmort_dist" | var == "cmort_res" | var == "gpp" | var == "npp"){
    var.to.plot <- change_unit_stoy(var.to.plot)
  }
  
  return(var.to.plot)
  
}

get_output_JULESRED <- function(site,run,var,file.dir,co2_levels){
  model_name = "JULES-RED"
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste(model_name, site,'-', run),
                            forcing.data = "cru_jra2.2")
  
  var.to.plot <- getField_DBEN(source.in, 
                               quant = get_quantity(var),
                               file.name = paste0(file.dir,
                                                  model_name,"_",var,"_PS_",co2_levels,"_",site,"_",run,".nc"),
                               model_name = model_name)
  
  # split up size classes graph to make more readable.
  # create shared maximum y axis to make graphs more interpretable. maybe change back later  
  ymax <- max(var.to.plot@data$Total)
  
  #[TODO] probably change unit later..
  if(var == "cmort" | var == "cmort_age" |
     var == "cmort_other" | var == "cmort_greff" ) {
    #Arthur must have changed these to -yr in his post-processing now (20.07.2023)
    var.to.plot@quant@units <- "kgC m-2 yr-1"
  }
  
  if(var =="npp"| var == "gpp"){
    var.to.plot <- change_unit_stoy(var.to.plot)
  }
  return(var.to.plot)
  
}

get_output_FATES <- function(site,run,var,file.dir,model_name,co2_levels ="412ppm"){
  model_name ="FATES"
  if(site == "FIN"){siteFates="fi"}
  if(site == "BIA"){siteFates ="bia"}
  if(site == "BCI"){siteFates ="bci"}
  if(run  == "P0"){
    runFates ="p0"
  }else{
    runFates="p1"
    if(as.numeric(strsplit(run,split=NULL)[[1]][2])>=1){ # evaluate for all run inputs > 1
      srun = as.numeric(strsplit(run,split=NULL)[[1]][2]) # retrieve sensitivity run value from string
    }
  }
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste(model_name, site,'-', run),
                            forcing.data = "cru_jra2.1")# different to the other models.
  
  if(runFates =="p0"){
    if(var =="cveg_total"){
      #name in netcdf is cveg, but it differs from the cveg file. so updating the quantity after reading in :
      var.to.plot <- getField_DBEN(source.in, 
                                   quant = get_quantity("cveg"),
                                   file.name = paste0(file.dir,co2_levels,"/",site,"/",runFates,"/",
                                                      model_name,"_",var,"_",runFates,"_",siteFates,".nc"),
                                   model_name = model_name)
      var.to.plot@quant <- get_quantity("cveg_total")
    } else if( var =="stemmort"){
      var.to.plot_under <- getField_DBEN(source.in, 
                                         quant = get_quantity(var),
                                         file.name = paste0(file.dir,co2_levels,"/",site,"/",runFates,"/",
                                                            model_name,"_",var,"_",runFates,"_",siteFates,"_understory.nc"),
                                         model_name = model_name)
      
      var.to.plot_over <- getField_DBEN(source.in, 
                                        quant = get_quantity(var),
                                        file.name = paste0(file.dir,co2_levels,"/",site,"/",runFates,"/",
                                                           model_name,"_",var,"_",runFates,"_",siteFates,"_overstory.nc"),
                                        model_name = model_name)
      
      var.to.plot  <- var.to.plot_under # "initalise" to keep metadata
      var.to.plot@data[,4:19]  <- var.to.plot_under@data[,4:19] + var.to.plot_over@data[,4:19] #combine # 06.06.2023 sizecalsses lost
      
    } else{
        var.to.plot <- getField_DBEN(source.in, 
                                     quant = get_quantity(var),
                                     file.name = paste0(file.dir,co2_levels,"/",site,"/",runFates,"/",
                                                        model_name,"_",var,"_",runFates,"_",siteFates,".nc"),
                                     model_name = model_name)
    }
 }
    
    if(runFates =="p1"   ){ 
      # run is P1 as in - sensitivity runs:
      if(var == "stemmort"){ # 
        var.to.plot_under <- getField_DBEN(source.in, 
                                           quant = get_quantity(var),
                                           file.name = paste0(file.dir,co2_levels,"/",site,"/",runFates,"/",
                                                              model_name,"_",var,"_",runFates,"_",siteFates,"_understory_",srun,".nc"),
                                           model_name = model_name)
        
        var.to.plot_over <- getField_DBEN(source.in, 
                                          quant = get_quantity(var),
                                          file.name = paste0(file.dir,co2_levels,"/",site,"/",runFates,"/",
                                                             model_name,"_",var,"_",runFates,"_",siteFates,"_overstory_",srun,".nc"),
                                          model_name = model_name)
        
        var.to.plot  <- var.to.plot_under # "initalise" to keep metadata
        var.to.plot@data[,4:19]  <- var.to.plot_under@data[,4:19] + var.to.plot_over@data[,4:19] #combine # 06.06.2023 sizecalsses lost
        
      } else{
        var.to.plot <- getField_DBEN(source.in, 
                                     quant = get_quantity(var),
                                     file.name = paste0(file.dir,co2_levels,"/",site,"/",runFates,"/",
                                                        model_name,"_",var,"_",runFates,"_",siteFates,"_",srun,".nc"),
                                     model_name = model_name)
        
      } # done reading in variables
   }
  
  # some unit adjustments
  if(var =="npp"| var == "gpp"){
    var.to.plot <- change_unit_stoy(var.to.plot)
  }
  
  #change units for nstem_size for 562ppm runs ( email Jessie 12/07/2023, 17:14)
  if(co2_levels=="562ppm" & var=="nstem_size" ){
    #must scale up to ha:
    var.to.plot@data[,4:20] <- var.to.plot@data[,4:20]*10000
  } 
  #AHES Sept 2024: seems like they also have to be changed for 412pp p0 runs in v3
  if(var =="nstem_size" | var =="stemmort"){
    #must scale up to ha:
    var.to.plot@data[,4:20] <- var.to.plot@data[,4:20]*10000
  }
  
  return(var.to.plot)
}

get_output_EDv3 <- function(site,run,var_in,file.dir,co2_levels = "412ppm"){
  model_name = "EDv3"
  
  if(co2_levels == "562ppm"){
    co2_levels = "PS_562ppm"
    #print(file.dir)
  }
  if(co2_levels == "412ppm"){
    co2_levels = "PS_412ppm"
    #print(file.dir)
  }
  
  if(run == "P0"){
    # benchmark run
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var_in,"_P0_",site,".nc")
  }else{ # sensitivity runs:
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var_in,"_PS_",site,"_",run,".nc")
  }
  
  print(file.name)
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste(model_name, site,'-', run),
                            forcing.data = "cru_jra2.2")
  
  var.to.plot <- getField_DBEN(source.in, 
                               quant = get_quantity(var_in),
                               file.name = file.name,
                               model_name = model_name)

  return(var.to.plot)
  
}
get_output_SEIBDGVM <- function(site,run,var,file.dir,co2_levels = "412ppm"){
  model_name = "seib"
  
  if(run == "P0"){
    #benchmark run
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var,"_P0_",site,".csv")
  }else{ # sensitivity runs:
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var,"_",run,"_",site,".csv")
  }
  
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste("SEIB-DGVM", site,'-', run),
                            forcing.data = "cru_jra2.2")
  
  var.to.plot <- getField_DBEN(source.in, 
                               quant = get_quantity(var),
                               file.name = file.name,
                               model_name = model_name)
  
  #reupdate metadata:
  
  if(var == "cmort" | var == "cmort_age" |
     var == "cmort_other" | var == "cmort_greff" ) {
    var.to.plot@quant@units <- "kgC m-2 yr-1"
  }
  
  if(var =="gpp"| var =="npp"){
    var.to.plot <- change_unit_stoy(var.to.plot)
  }
  return(var.to.plot)
  
}

get_output_ORCHIDEE <- function(site,run,var,file.dir,co2_levels = "412ppm"){
  model_name="ORCHIDEE"
  
  if(run == "P0"){
    #benchmark run
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var,"_P0_",site,".nc")
  }else{ # sensitivity runs:
    file.name = paste0(file.dir,co2_levels,"/",model_name,"_",var,"_",run,"_",site,".nc")
  }
  
  # set metadata 
  source.in <- defineSource(id = site,
                            dir = file.dir ,
                            format = DBEN,
                            name = paste(model_name, site,'-', run),
                            forcing.data = "crujra.2.2")
  
  var.to.plot <- getField_DBEN(source.in, 
                               quant = get_quantity(var),
                               file.name = file.name,
                               model_name = model_name)
  
  #[TODO] probably change unit later..
  #if(var == "cmort"){
  #  #convert to -yr:
  #  var.to.plot <- change_unit_stoy(var.to.plot)
  #}
  
  return(var.to.plot)
}

#temporary function, to change units for some fluxes where needed
change_unit_stoy <- function(var_to_plot){
  
  seconds_to_year = 60*60*24*365
  #convert to -yr:
  end <- dim(var_to_plot@data)[2]
  var_to_plot@data[,4:end] <- var_to_plot@data[,4:end]  * seconds_to_year
  var_to_plot@quant@units <- "kgC m-2 yr-1"
  return(var_to_plot)
}

### helper function for Sensitivity analysis ( disturbance analysis) plotting:
#collect_ambient_x =  data frame containing xaxis-values for ambient simulation output for variable of choice
#collect_elevated_x = xaxis elevated simulation output for variable of choice
#collect_ambient_y = yaxis ambient simulation output for variable of choice
#collect_elevated_y = yaxis elevated simulation output for variable of choice
#model_name name of model to plot "FATES", ""BiomeE-Standalone" etc..
#ymax = upper yaxis limit
#site name "FIN","BIA","BCI"
# xmax = upper xaxis limit

plot_SA_output <- function(collect_ambient_x,collect_ambient_y,collect_elevated_x,collect_elevated_y,model_name,ymax,site,xmax,col="site"){
  
  #more flexible comparison between cmort and WBgrowth
  if(col=="grey"){
    col="grey"
  }else{
    col= get_biome_colour(site)
  }
  plot(collect_ambient_x[which(collect_ambient_x$model_name == model_name),site],
       collect_ambient_y[which(collect_ambient_y$model_name == model_name),site], 
       col=col,ylim=c(0,ymax), ylab="", xlab="", xaxt="n", yaxt="n", xlim=c(0,xmax), pch=get_model_pch(model_name), lty=1)
  lines(collect_ambient_x[which(collect_ambient_x$model_name == model_name),site],
        collect_ambient_y[which(collect_ambient_y$model_name == model_name),site], 
        col=col,  lty=1)
  lines(collect_elevated_x[which(collect_elevated_x$model_name == model_name),site],
        collect_elevated_y[which(collect_elevated_y$model_name == model_name),site],
        col=col,  pch=8, type="p")
  lines(collect_elevated_x[which(collect_elevated_x$model_name == model_name),site],
        collect_elevated_y[which(collect_elevated_y$model_name == model_name),site],
        col=col,lty=2)
}

# not all naming conventions were adhered to. 
# this function is used to adress this:
# map individual model outputs to D-BEN run number, 
# function is used by get_model_output to automatise 
# the calling of the same runs across models.
get_model_run_number <- function(model_name,protocol_run_name){
  
  # create dataframe that collects all models and the run number/name used:
  df = data.frame(Mnames = models_avail)
  tmp <- as.data.frame(matrix(NA,ncol=8,nrow= length(models_avail)))
  names(tmp) <-  c("P0","P1","P2","P3","P4","P5","P6","P7")
  df <- cbind(df, tmp)
  
  df[df$Mnames=="FATES",2:9]             <- c("P0","P1","P2","P3","P4","P5","P6","P7") # P0 and P1 are equivalent in parameterisation, but P0 is not run for eCO2, so no output exists, but the identical
  # output exists for P1. Some variables only exist for P0, which right now cannot be accessed thugh this setup.
  # first index P0 chenaged from P1 AHES, Sept 2024 in v5
  
  df[df$Mnames=="LPJ-GUESS",2:9]         <- c("P0","P1","P2","P3","P4","P5","P6","P7")
  df[df$Mnames=="BiomeEP",2:9]           <- c("P0","P0","PS1","PS2","PS3","PS4","PS5","PS6") #no P0 run, equivalent to PS1
  df[df$Mnames=="BiomeE-Standalone",2:9] <- c("00","00","01", "02", "04", "08", "20", "40") #no P0 run, equivalent to PS1"
  df[df$Mnames=="JULES-RED",2:9]         <- c("1","1","2","3","4","5","6","7") # no P0 run, because that is equivalent to PS1
  df[df$Mnames=="CABLE-POP",2:9]         <- c("0","1","2","3","4","5","6","7")
  df[df$Mnames =="EDv3",2:9]             <- c("P0","1","2","3","4","5","6","7")
  df[df$Mnames =="SEIB-DGVM",2:9]        <- c("P0","P1","P2","P3","P4","P5","P6","P7") # [TODO] 17.07. checked with H.Sato on PS runs ( only go up to 6), waiting for reply.
  df[df$Mnames =="ORCHIDEE",2:9]        <- c("P0","P0","P0","P0","P0","P0","P0","P0") # [TODO] 11.09.2023 SEnsitivity runs not yet incorporates in output post-processing, only allow P0
  
  # FATES' model output calling function already is set up in a way that it works with default run names.
  # [TODO] re-naming of biomeEP scenarios once PS1 is run.
  
  # to deal with different file naming conventions:
  runs = c("P0","P1","P2","P3","P4","P5","P6","P7")
  runs_biomeEP <- c("P0","P0","PS1","PS2","PS3","PS4","PS5","PS6")
  
  # disturbance: stand-replacing (resetting to initial conditions) and stochastic with mean frequency
  # of:0.01, 0.02, 0.04, 0.08, 0.20, 0.40 (corresponding to the file names: _01, _02, _04, _08, _20, _40)
  runs_biomeES <- c("00","01", "02", "04", "08", "20", "40","08")
  runs_julesred = c("1","2","3","4","5","6",NA)
  runs_cablepop = c("00","1","2","3","4","5","6","7")
  if(model_name =="FATES" & protocol_run_name =="P0"){
    print("warning, only versions >v3 of FATES have P0 runs. for earlier versions, change the first index in this function for FATES to P1")
  }
  return(df[df$Mnames==model_name,protocol_run_name])
}

#load selected model's output
get_model_output <- function(var,site,run, model_name,co2_levels){
  
  ########################LPJGUESS####################################
  if(model_name =="LPJ-GUESS"){
    if(co2_levels == "412ppm"){
      co2_levels = "PS_412ppm/"
    } 
    if(co2_levels == "562ppm"){
      co2_levels = "PS_562ppm/"
      # print(paste("testing",co2_levels))
    }
    
    model_out <- get_output_LPJGUESS(site, run = get_model_run_number(model_name,run), 
                                     var = var, file.dir = file.dir.lpjguess, co2_levels = co2_levels)
    
    ## print(model_name)
  }
  
  ########################FATES####################################
  if(model_name =="FATES"){
    #print(co2_levels)
    model_out    <- get_output_FATES(site,run = get_model_run_number(model_name,run),
                                     var = var,file.dir.fates,co2_levels = co2_levels)
    # print(model_name)
  }
  
  ########################BIOME_EP####################################
  if(model_name == "BiomeEP"){
    # print(model_name)
    if(var =="cmort"){ # filename called cmort_size in BIomeEP, but cmort in DBEN, to make sure I can loop across DBEN variables easily, I added this switch:
      var_loc = "cmort"
    } else{
      var_loc = var
    }
    model_out <- get_output_BiomeEP(site,run = get_model_run_number(model_name,run),
                                    var = var_loc,file.dir = file.dir.biomeEP,co2_levels = co2_levels)
    
  }
  ########################BIOME_ES####################################
  if(model_name == "BiomeE-Standalone"){
    
    if(var=="AGcwood"){
      var_loc = "AGB" #It is above ground wood biomass. The biomass of leaves is low, 0.2~0.4 kg C m-2.
    }else if(var =="cmort"){
      var_loc = "cmort_size"
    }else{
      var_loc=var
    }
    #print(var_loc)
    model_out <- get_output_BiomeE_standalone(site,run = get_model_run_number(model_name,run),
                                              var = var_loc, file.dir.biomeE_Standalone,co2_levels = co2_levels)
    
    # print(model_name)
  }
  ########################CABLEPOP####################################
  if(model_name == "CABLE-POP"){
    if(co2_levels == "412ppm"){
      co2_levels = "PS_412ppm"
    }
    if(co2_levels == "562ppm"){
      co2_levels = "PS_562ppm"
    }
    
    if(var=="cmort"){ # account for the fact the CABLE does not provide total mortality:
      cmortdist_out <- get_output_CABLEPOP(site,run = get_model_run_number(model_name,run),
                                           var = "cmort_dist",file.dir.cablepop,co2_levels = co2_levels)
      cmortcrowd_out <- get_output_CABLEPOP(site,run = get_model_run_number(model_name,run),
                                            var = "cmort_crowd",file.dir.cablepop,co2_levels = co2_levels)
      cmortres_out <- get_output_CABLEPOP(site,run = get_model_run_number(model_name,run),
                                          var = "cmort_res",file.dir.cablepop,co2_levels = co2_levels)
      
      #recycle one of the above objects, but change metadata accordingly
      model_out             <- cmortcrowd_out
      model_out@data$Total  <- cmortres_out@data$Total + cmortdist_out@data$Total + cmortcrowd_out@data$Total
      model_out@quant       <- get_quantity("cmort")
      model_out@quant@units <- "kgC m-2 yr-1"
      
    }else{
      
      model_out <- get_output_CABLEPOP(site,run = get_model_run_number(model_name,run),
                                       var = var,file.dir.cablepop,co2_levels = co2_levels)
    }
    
    # print(model_name)
  }
  
  ########################JULESRED####################################
  if(model_name =="JULES-RED"){
    
    if(var=="AGcwood"){
      var_loc= "AGB" #leaves negligible at BIA and BCI, and for FIN, made so leaves are excluded ( email 4.3.2023)
    }else{
      var_loc=var
    }
    model_out <- get_output_JULESRED(site,run = get_model_run_number(model_name,run),
                                     var = var_loc,file.dir.julesred,co2_levels = co2_levels)
    
    # print(model_name)
  }
  
  ########################EDv3####################################
  if(model_name =="EDv3"){
    
    # account for EDv3 not reporting a cmort variable:
    if(var == "cmort"){
      model_out_disb  <- get_output_EDv3(site,run = get_model_run_number(model_name,run),
                                         var = "cmort_disb",file.dir.EDv3,co2_levels = co2_levels)
      model_out_mort  <- get_output_EDv3(site,run = get_model_run_number(model_name,run),
                                         var = "cmort_mort",file.dir.EDv3,co2_levels = co2_levels)
      
      model_out <- model_out_mort
      model_out@data[,4:20] <- model_out_mort@data[,4:20] + model_out_disb@data[,4:20]
      model_out@quant <- get_quantity(var)
      
    }else if(var == "stemmort"){
      model_out_disb  <- get_output_EDv3(site,run = get_model_run_number(model_name,run),
                                         var = "stemmort_disb",file.dir.EDv3,co2_levels = co2_levels)
      model_out_mort  <- get_output_EDv3(site,run = get_model_run_number(model_name,run),
                                         var = "stemmort_mort",file.dir.EDv3,co2_levels = co2_levels)
      
      
      model_out   <- model_out_mort
      model_out@data[,4:20]  <- model_out_mort@data[,4:20] + model_out_disb@data[,4:20]
      model_out@quant   <- get_quantity(var)
      
    }else{
      
      model_out    <- get_output_EDv3(site,run = get_model_run_number(model_name,run),
                                      var = var,file.dir = file.dir.EDv3,co2_levels = co2_levels)
      
    }
    # print(model_name)
  }
  
  ########################SEIB-DGVM####################################
  if(model_name =="SEIB-DGVM"){
    # print(co2_levels)
    model_out    <- get_output_SEIBDGVM(site,run = get_model_run_number(model_name,run),
                                        var = var,file.dir.seibdgvm,co2_levels = co2_levels)
    # print(model_name)
  }
  
  ########################ORCHIDEE####################################
  if(model_name == "ORCHIDEE"){
    #print(var)
    if(var =="cwood_size"){
      var_loc = "cwood_size"# in v5: cwood variable now available. moved away from "AGcwood_size"
      
    } else{
      var_loc = var
    }
    #print(var_loc)
    model_out <- get_output_ORCHIDEE(site,run = get_model_run_number(model_name,run),
                                     var = var_loc,file.dir = file.dir.orchidee,co2_levels = co2_levels)
    
  }
  
  ###########################################################
  return(model_out)
  
}

# aesthetics for consistent plotting
# x allowed to be biome name or site name, for more convenient use
get_biome_colour <- function(x, black_only=FALSE){
  
  # preparation:
  biomes  <- c("Boreal","Boreal","Temperate","Tropics")
  sites    <- c("FIN","FI", "BIA" ,"BCI")
  cols   <- c("dark green","dark green","light green","brown")
  
  # to make compatible for all use-cases in this script:
  if(black_only){
    cols   <- c(1,1,1,1)
  } 
  
  #automatic colour coding for biome based on site name:
  if(x %in% sites){
    #retrieval:
    ret_col <- cols[sites==x]
  }
  
  if(x %in% c("Boreal","Temperate","Tropics")){
    #retrieval:
    ret_col <- cols[biomes==x]
  }
  
  return(ret_col)
  
}

#turn into object with function that returns colour when giving it model-name later..
get_model_colour <- function(model_in){
  # PREPARATION:
  models  <- c( "JULES-RED"  ,       "FATES"  ,           "ORCHIDEE"  ,        "LPJ-GUESS" ,        "CABLE-POP" ,  
                "BiomeEP"  ,         "BiomeE-Standalone", "SEIB-DGVM" ,"EDv3")
  
  
  # adress colour-blindness:
  #palette.colors(8)
  model_cols <- c(palette.colors()[c(7,5,6,3,9,8,2)],"#BCE9C5",palette.colors()[c(4)])
  
  #model_cols <- c("red","orange","blue","light blue","grey","purple","magenta","green","dark green")
  plot_mod_cols <- data.frame(models,model_cols)
  
  #RETRIEVAL:
  ret_col <-  plot_mod_cols[ models %in% model_in , ]$model_cols
  return(ret_col)
}

get_model_pch <- function(model){
  # PREPARATION:
  models <- c("FATES","JULES-RED","ORCHIDEE","LPJ-GUESS","CABLE-POP","BiomeEP","BiomeE-Standalone")
  model_pch <- c(3,11,1,18,12,6,8)
  plot_mod_pch <- data.frame(models,model_pch)
  
  #RETRIEVAL:
  ret_pch <-   plot_mod_pch[which(plot_mod_pch$models == model), ]$model_pch
  return(ret_pch)
}

#stand-structure related shared objects:
sc <- c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90" ,  "<100",  "<150" , "<200",  ">=200")

#to deal with different file naming conventions:
runs = c("P1","P2","P3","P4","P5","P6","P7")
runs_biomeEP <- c("P0","PS1","PS2","PS3","PS4","PS5","PS6")

#Disturbance: stand-replacing (resetting to initial conditions) and stochastic with mean frequency
#of:0.01, 0.02, 0.04, 0.08, 0.20, 0.40 (corresponding to the file names: _01, _02, _04, _08, _20, _40)
runs_biomeES <- c("00","01", "02", "04", "08", "20", "40")
runs_julesred = c("1","2","3","4","5","6","7")
runs_cablepop = c("1","2","3","4","5","6","7")
#[TODO] re-naming of biomeEP scenarios once PS1 is run.

#biome-specific points on graphs:
Tropics_pch = 16
Temperate_pch = 16
Boreal_pch = 16

biomes <- c("Boreal" , "Temperate" ,"Tropics"  )

# function to plot a model's  stand sturcutre output alongside observations.
# nstem_size can be plotted on log yaxis (log = TRUE) or not
# ylim_ext can override automatically generated ylim  values.
# model name must match the model name used for get_model_colour
# var the DBEN variable that should be plotted.
# plot_site_name, plots site name in each plot. Optional to make backwards compatible, and to have more flexibility with plotting in publication.
# plot_year introduce the option to plot any simulation year's standstucture, # if no year is given, it uses the last-simulation year's stand structure.
stand_structure_benchmarks <- function(model_name = "LPJ-GUESS",site_in = site,takelog =TRUE, 
                                       var_in = var, ylim_ext = NULL, stand_structure_obs_in = stand_structure_obs,  
                                       run="P1", plot_site_name=TRUE, plot_year = NULL){
  
  # define benchmark run:  must be done here, and then sent through get_model_run_number, 
  # because the models call their runs differently. the get_model_run_number function maps it back to the relevant D-BEN run.
  
  # account for some data not having dbh-classes, but models must be plotted for it. here, retrieve complete set of dbh-classes:
  dbh_classes_plot <- stand_structure_obs_in[which(stand_structure_obs_in$site == "BCI"),]$dbh_classes_num
  
  site=site_in
  
  # call each model: 
  model_out <- get_model_output(site=site, run= run, var = var_in, model_name=model_name, co2_levels = "412ppm")
  
  # if no year  is given for plotting of the stand structure,
  # use the last simulation year:
  if(is.null(plot_year)){
    plot_year = max(model_out@data$Year)
  }
  
  if(model_name=="JULES-RED"){
    # JULES-RED doesn't have a sizeclass <1, so adressing this here, that the plotting works:
    # add sizeclass colummn
    model_out@data$`<1` <- NA
    # reorder dataframe:
    model_out@data <-  model_out@data[,c( "Year" , "Lat" ,  "Lon", "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90", "<100", "<150", "<200", ">=200" ,"Total")]
  }
  
  ############### ############### ############### ############### ############### ############### 
  if(var_in == "nstem_size"){
    # plot models - nstem_size
    
    # subset of dbh- classes for plotting, because in the observations there are no more than that.
    # [TODO]would ideally subset with this vector of column names, but it somehow doesnt, work. doing this manually for now, but ther emust be a better solution..
    dbh_classes_sel <-  c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90", "<100", "<150", "<200", ">=200" )
    
    # make ylims flexible, or prescribe from external:
    max_model_output <- max(model_out@data[,
                                           c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                              "<90" ,  "<100",  "<150" , "<200",  ">=200")])
    min_model_output <-  max(model_out@data[,
                                            c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                               "<90" ,  "<100",  "<150" , "<200",  ">=200")])
    max_obs <- max( na.omit(stand_structure_obs_in[,c("AGB_size_kgCm.2","AGB_size_upper_kgCm.2","AGB_size_lower_kgCm.2")]) )
    min_obs <- min( na.omit(stand_structure_obs_in[,c("AGB_size_kgCm.2","AGB_size_upper_kgCm.2","AGB_size_lower_kgCm.2")]) )
    
    
    #check if external ylims are present or not:
    if(is.null(ylim_ext)){ # not present, obtain ylim_set from within
      ylim_set = c(min(min_obs,min_model_output),max(max_obs,max_model_output))
    }else{
      ylim_set = ylim_ext # present. pass external ylims to ylim_set
    }
    
    #if nstem should be logged:
    if(takelog ==TRUE){
      
      if(is.null(ylim_ext)){ # not present, obtain ylim_set from within
        ylim_set = log(ylim_set)
        ylim_set[ylim_set ==-Inf] <- -2
      }
      
      plot(dbh_classes_plot,
           log(model_out@data[model_out@data$Year == plot_year, 
                              c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                 "<90" ,  "<100",  "<150" , "<200",  ">=200")]), 
           ylim = ylim_set ,col= get_model_colour(model_name),pch=16,cex=1,ylab="",type="p") 
      
      points(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num, 
             log(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$nstem_size_ha.1),  
             col= get_biome_colour(site))
      arrows(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num,
             log(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$nstem_size_lower_ha.1),
             stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num,
             log(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$nstem_size_upper_ha.1),
             length = 0.00, angle = 90, code = 3,
             col= get_biome_colour(site) )
      
      if(plot_site_name){
        mtext(site,adj=0.95,side=3,line=-1.3)
      }
      mtext(side=2,"log(nstem_size) (nstem)",line=2,outer=FALSE)
      
    }else {  # nstem not logged
      
      plot(dbh_classes_plot,
           model_out@data[model_out@data$Year == plot_year, c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                                               "<90" ,  "<100",  "<150" , "<200",  ">=200")], 
           ylim = ylim_set ,col= get_model_colour(model_name),pch=16,cex=1,type="p") 
      
      points(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num, 
             stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$nstem_size_ha.1,  
             col= get_biome_colour(site))
      arrows(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num,
             stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$nstem_size_lower_ha.1,
             stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num,
             stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$nstem_size_upper_ha.1,
             length = 0.00, angle = 90, code = 3,
             col= get_biome_colour(site) )
      
      if(plot_site_name){
        mtext(site,adj=0.95,side=3,line=-1.3)
      }
      mtext(side=2,"nstem_size (nstem)",line=2,outer=FALSE)
    }#stem not logged
    
  }
  
  ############### ############### ############### ############### 
  if(var_in =="cwood_size"){
    # !Note, applying a 0.75 factor to arrive at AGcwood_size; except for ORCHIDEE
    if(model_name=="ORCHIDEE"){roots_off=1.0}else{roots_off = 0.75}
    
    # cwood plots
    max_model_output <- max(model_out@data[, c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                                "<90" ,  "<100",  "<150" , "<200",  ">=200")]*roots_off)
    max_obs <- max( na.omit(stand_structure_obs_in[,c("AGB_size_kgCm.2","AGB_size_upper_kgCm.2","AGB_size_lower_kgCm.2")]) )
    
    # convert to AGcwood:
    
    # check if external ylims are present or not:
    if(is.null(ylim_ext)){ # not present, obtain ylim_set from within
      ylim_set = c(0,max(max_obs,max_model_output) + max(max_obs,max_model_output)*0.1) 
    }else{
      ylim_set = ylim_ext # present. pass external ylims to ylim_set
    }
    
    plot(dbh_classes_plot,
         model_out@data[model_out@data$Year == plot_year, c( "<1" ,   "<5" , "<10", "<15",  "<20" , 
                                                             "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90" ,  "<100",  "<150" , "<200",  ">=200")]*roots_off, 
         ylim = ylim_set ,col= get_model_colour(model_name),pch=16,cex=1,ylab="",type="p") 
    
    points(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num, 
           stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$AGB_size_kgCm.2,  col= get_biome_colour(site))
    arrows(stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num,
           stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$AGB_size_lower_kgCm.2,
           stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$dbh_classes_num,
           stand_structure_obs_in[which(stand_structure_obs_in$site == site),]$AGB_size_upper_kgCm.2,
           length = 0.00, angle = 90, code = 3,
           col= get_biome_colour(site) )
    
    if(plot_site_name){
      mtext(site,adj=0.95,side=3,line=-1.3)
    }
    mtext(side=2,"AGcwood (kgC m-2)",line=2,outer=FALSE)
  }
  
}

# function to plot a model's  stand sturcutre output alongside observations.
# nstem_size can be plotted on log yaxis (log = TRUE) or not
# ylim_ext can override automatically generated ylim  values.
# model name must match the model name used for get_model_colour
# var the DBEN variable that should be plotted.
stand_structure <- function(model_name = "LPJ-GUESS",site_in = site,takelog =TRUE, run ="P1",var_in = "cwood_size", ylim_ext = NULL,co2_levels ="412ppm",add=FALSE,main_in = NULL){
  
  # define benchmark run:  must be done here, and then sent through get_model_run_number, 
  # because the models call their runs differently. the get_model_run_number function maps it back to the relevant D-BEN run.
  
  testing=FALSE
  #call each model: 
  model_out <- get_model_output(site=site_in, run= run, var = var_in, model_name=model_name, co2_levels = co2_levels)
  if(co2_levels == "412ppm"){
    pch = 16
  }else{
    pch = 8
  }
  
  if(add==TRUE){
    par(new=TRUE)
  }
  if(!is.null(main_in)){
    main = run
  }else{
    main=NULL
  }
  
  if(model_name=="JULES-RED"){
    #JULES-RED doesn't seem to have a sizeclass <1, so adressing this here, that the plotting works:
    #add sizeclass colummn
    model_out@data$`<1` <- NA
    #reorder dataframe:
    model_out@data <-  model_out@data[,c( "Year" , "Lat" ,  "Lon", "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90", "<100", "<150", "<200", ">=200" ,"Total")]
  }
  
  if(testing==TRUE){
    if(model_name == "LPJ-GUESS" ){
      model_out <- get_output_LPJGUESS(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.lpjguess,co2_levels = "PS_412ppm/")
    }
    if(model_name == "BiomeEP"){
      model_out <- get_output_BiomeEP(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.biomeEP,co2_levels = "412ppm")
    }
    if(model_name == "SEIB-DGVM"){
      model_out <- get_output_SEIBDGVM(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.seibdgvm,co2_levels = "412ppm")
    }
    if(model_name == "EDv3"){
      model_out <- get_output_EDv3(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.EDv3,co2_levels = "412ppm")
    }
    if(model_name == "FATES"){
      model_out <- get_output_FATES(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.fates,co2_levels = "412ppm")
    }
    if(model_name == "CABLE-POP"){
      model_out <- get_output_CABLEPOP(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.cablepop,co2_levels = "PS_412ppm")
    }
    if(model_name == "BiomeE-Standalone"){
      model_out <- get_output_BiomeE_standalone(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.biomeE_Standalone,co2_levels = "412ppm")
    }
    if(model_name == "ORCHIDEE"){
      # model_out_default <- get_output_ORCHIDEE(site = site_in,run = get_model_run_number(model_name, run) ,var = var_in,file.dir = file.dir.orchidee,co2_levels = "412ppm")
      model_out_testing <- get_model_output(site=site_in, run= run, var = var_in, model_name=model_name, co2_levels = "412ppm")
      
    }
  }
  
  ############### ############### ############### ############### ############### ############### 
  if(var_in == "nstem_size"){
    # plot models - nstem_size
    
    # subset of dbh- classes for plotting, because in the observations there are no more than that.
    # [TODO]would ideally subset with this vector of column names, but it somehow doesnt, work. doing this manually for now, but ther emust be a better solution..
    dbh_classes_sel <-  c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90", "<100", "<150", "<200", ">=200" )
    dbh_classes_num <- c(1,5,10,15,20,30,40,50,60,70,80,90,100,150,200,210)
    
    # make ylims flexible, or prescribe from external:
    max_model_output <- max(model_out@data[,
                                           c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                              "<90" ,  "<100",  "<150" , "<200",  ">=200")],na.rm = TRUE)
    min_model_output <-  max(model_out@data[,
                                            c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                               "<90" ,  "<100",  "<150" , "<200",  ">=200")],na.rm = TRUE)
    
    # check if external ylims are present or not:
    if(is.null(ylim_ext)){ # not present, obtain ylim_set from within
      ylim_set = c(min(min_model_output,na.rm = TRUE),max(max_model_output,na.rm = TRUE))
    }else{
      ylim_set = ylim_ext # present. pass external ylims to ylim_set
    }
    
    # if nstem should be logged:
    if(takelog ==TRUE){
      
      if(is.null(ylim_ext)){ # not present, obtain ylim_set from within
        ylim_set = log(ylim_set)
        ylim_set[ylim_set ==-Inf] <- -2
      }
      
      if(add==TRUE){
        plot(dbh_classes_num,
             log(model_out@data[model_out@data$Year == max(model_out@data$Year), 
                                c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                   "<90" ,  "<100",  "<150" , "<200",  ">=200")]), 
             ylim = ylim_set ,col= get_model_colour(model_name),pch=pch,cex=1,ylab="",yaxt="n")
        
      }else{
        plot(dbh_classes_num,
             log(model_out@data[model_out@data$Year == max(model_out@data$Year), 
                                c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                   "<90" ,  "<100",  "<150" , "<200",  ">=200")]), 
             ylim = ylim_set ,col= get_model_colour(model_name),pch=pch,cex=1,ylab="") 
        
      }
      
      mtext(site,adj=0.95,side=3,line=-1.3)
      mtext(side=2,"log(nstem_size) (nstem)",line=2,outer=FALSE)
      
    }else {  # nstem not logged
      
      plot(dbh_classes_num,
           model_out@data[model_out@data$Year == max(model_out@data$Year), c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                                                              "<90" ,  "<100",  "<150" , "<200",  ">=200")], 
           ylim = ylim_set ,col= get_model_colour(model_name),pch=pch,cex=1)
      
      mtext(site,adj=0.95,side=3,line=-1.3)
      mtext(side=2,"nstem_size (nstem)",line=2,outer=FALSE)
    }#stem not logged
    
  }
  
  ############### ############### ############### ############### 
  if(var_in =="cwood_size"){
    
    dbh_classes_num <- c(1,5,10,15,20,30,40,50,60,70,80,90,100,150,200,210)
    
    # cwood plots
    max_model_output <- max(model_out@data[, c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                                "<90" ,  "<100",  "<150" , "<200",  ">=200")],na.rm = TRUE)
    
    #check if external ylims are present or not:
    if(is.null(ylim_ext)){ # not present, obtain ylim_set from within
      ylim_set = c(0,max(max_model_output,na.rm = TRUE) + max(max_model_output,na.rm = TRUE)*0.1) 
    }else{
      ylim_set = ylim_ext # present. pass external ylims to ylim_set
    }
    
    if(add==TRUE){
      plot(dbh_classes_num,
           model_out@data[model_out@data$Year == max(model_out@data$Year), c( "<1" ,   "<5" , "<10", "<15",  "<20" , 
                                                                              "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90" ,  "<100",  "<150" , "<200",  ">=200")], 
           ylim = ylim_set ,col= get_model_colour(model_name),pch=pch,cex=1.5,ylab="",yaxt="n", main = main) 
      
    }else{
      plot(dbh_classes_num,
           model_out@data[model_out@data$Year == max(model_out@data$Year), c( "<1" ,   "<5" , "<10", "<15",  "<20" , 
                                                                              "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" ,  "<90" ,  "<100",  "<150" , "<200",  ">=200")], 
           ylim = ylim_set ,col= get_model_colour(model_name),pch=pch,cex=1.5,ylab="", main = main) 
      
    }
    #mtext(site,adj=0.95,side=3,line=-1.3)
    
    #mtext(side=2,"AGcwood (kgC m-2)",line=2,outer=FALSE)
  }
  
}

# function to plot a single models' regrowth dynamics alongside regrowth obs and equilibrium dynamics obs.
regrowth_benchmarks <- function(model_name = "LPJ-GUESS",site_in = site, var_in = var,  regrowth_obs_in = regrowth_obs, eq_dyn_obs_in = eq_dyn_obs){
  
  # define benchmark run:  must be done here, and then sent through get_model_run_number, 
  # because the models call their runs differently. the get_model_run_number function (within get_model_output()) maps it back to the relevant D-BEN run.
  run="P0"
  
  # call each model: 
  model_out <- get_model_output(site=site_in, run= run, var = var_in, model_name=model_name, co2_levels = "412ppm")
  
  # plot model:
  last_sim_year = max(model_out@data$Year)
  plot(1:dim(model_out@data)[1], model_out@data$Total, type = "l", ylim = c(0,25),xlim = c(0,last_sim_year+10), col = get_model_colour(model_name),lwd=1.5)
  
  # prepare equilibroum dynamics observations for plotting within regrowth observations, at the end of the regrowth period, where a dynamic equilibrium has instilled by most models.
  eq_dyn_obs_in <- eq_dyn_obs_in %>% group_by(site) %>% mutate(new_plot_loc = last_sim_year + (Year - max(Year)))
  
  # load and plot obs:
  biome <- biomes[which(sites==site)]
  points(regrowth_obs_in[which(regrowth_obs_in$Biome == biome),]$bin_num +30, regrowth_obs_in[which(regrowth_obs_in$Biome == biome),]$AGcwood_kgCm2_med,type="p")
  arrows(regrowth_obs_in[which(regrowth_obs_in$Biome == biome),]$bin_num +30, regrowth_obs_in[which(regrowth_obs_in$Biome == biome),]$AGcwood_kgCm2_10, 
         regrowth_obs_in[which(regrowth_obs_in$Biome == biome),]$bin_num +30 , regrowth_obs_in[which(regrowth_obs_in$Biome == biome),]$AGcwood_kgCm2_90, 
         length = 0.00, angle = 90, code = 3,col = get_biome_colour(biome) )
  mtext(site,adj=0.95,side=3,line=-1.3)
  
  # eq_dyn:
  points(eq_dyn_obs_in[which(eq_dyn_obs_in$site == site),]$new_plot_loc, eq_dyn_obs_in[which(eq_dyn_obs_in$site == site),]$AGB_kgCm2)
  arrows(eq_dyn_obs_in[which(eq_dyn_obs_in$site == site),]$new_plot_loc, eq_dyn_obs_in[which(eq_dyn_obs_in$site == site),]$AGB_lower_kgCm2, 
         eq_dyn_obs_in[which(eq_dyn_obs_in$site == site),]$new_plot_loc, eq_dyn_obs_in[which(eq_dyn_obs_in$site == site),]$AGB_upper_kgCm2, 
         length = 0.00, angle = 90, code = 3,col = get_biome_colour(biome) )
  mtext(site,adj=0.95,side=3,line=-1.3)
  
}

# convenience function that removes the firs 30 or so years from each model simulation, 
# so that only the regrowth phase is considered in further analysis.
# only for P0 runs!!!
# some of these inputs could be automatically extracted from the metadata, but moving on..
omit_equilibrium_phase <- function(data_object,model_name,site,selection_by_year=TRUE){
  
  if(selection_by_year){
    
    # it seems that models have set their first regrowth year to different timings, so I am here manually omitting the equilibrium and disturbance year.
    # for LPJG the first years will inevitably have 0 trees present until they establish, sometimes a 0 value will therefore cause trouble.
    if(model_name=="JULES-RED"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==31)# disturbance at end of 30.
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="FATES"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==34)# disturbance in year 32
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="ORCHIDEE"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==32)# disturbance i year 31
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="LPJ-GUESS"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==32) #disturbance flux at year 31
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="CABLE-POP"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==31)
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="BiomeEP"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==31)
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="BiomeE-Standalone"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==32) # disturbance mortality fluxes recorded in year 31.
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="SEIB-DGVM"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==32) 
      data_object@data <- data_object@data[lower:maxlength,]
    }
    
    if(model_name=="EDv3"){
      maxlength <- dim(data_object@data)[1]
      lower <- which(data_object@data$Year ==31)
      data_object@data <- data_object@data[lower:maxlength,]
    }
  }
  
  if(!selection_by_year){
    # use gpp to extract lowest point ( i.e. identify year from which regrowth commences):
    # This can later be manually adjusted to what year each model actually considers their first regrowth year. 
    # Can be different, depending on the assumptions that were made about the disturbance.
    gpp <- get_model_output(var="gpp",site= site,run = "P0",model_name = model_name,co2_levels = "412ppm") # not a good variable, because it contains grasses gpp.
    
    # Identify disturbance year (year of lowest gpp values) this is the lower index for data
    # selection. The Second (nested) check ([which(Year>2),]) for FATES output, has strange low
    # values at start of simulation
    lower <- which(gpp@data[which(Year>2),]$Total==min(gpp@data[which(Year>2),]$Total))[1]+2
    
    if((model_name =="ORCHIDEE") | (model_name =="FATES")| (model_name =="EDv3")){
      #somehow I cannot use ORCHIDEE output for this meaningfully, so doing by hand:
      #also for FATES and ED:
      maxlength <- dim(data_object@data)[1]
      data_object@data <- data_object@data[33:maxlength,]
    }else{   # all models, except EDv3,ORCHIDEE AND FATES 
      
      if((model_name =="JULES-RED") | (model_name=="BiomeEP")){
        maxlength <- dim(data_object@data)[1]
        data_object@data <- data_object@data[31:maxlength,]
      }else{
        maxlength <- dim(data_object@data)[1]
        data_object@data <- data_object@data[lower:maxlength,]
      }
      
    }
  }
  
  #update first year in timeseries
  data_object@first.year = as.numeric(data_object@data[1,"Year"])
  
  return(data_object)
  
}

# can be used to calculate % Carbon and % individuals lost of total
calc_mort_rate_andrusetal2021 <- function(Td,Tl,Pe=1){
  # https://www-sciencedirect-com.ludwig.lub.lu.se/science/article/pii/S0048969721066808#bb0020
  # In the period 1993–2013, we calculated the mortality rate (M) of each plot following the Eq. (1) as done by Andrus et al. (2021):
  # M = (Td/Tl/Pe)*100
  # where Td represents the number of trees that died (combined or by taxon) during time period Pe and Tl represents the number of living trees in the beginning of each year (or census period Pe).
  
  M = Td/Tl/Pe*100
  #account fo rthe fact that sometimes there are 0 trees, and must not divide by 0
  M[which(is.nan(M))] <- 0
  M[which(is.infinite(M))] <- 0
  
  # if(all(Tl > Td)){
  #    error("dying tree number larger than alive tree number, check input")
  #          }
  return(M)
}

# function to override "Total" column when fewer columns are selected to calculate the total.
# here we update "total" column, now that we have excluded the firest 3 sizeclasses. 
update_total <- function(model_out_df){
  ns <- names(model_out_df)
  if("Total" %in% ns){
    model_out_df$Total <-rowSums(model_out_df[, 4:(dim(model_out_df)[2]-1) ],na.rm = TRUE)
    print("Total column found, updating Total. note, the function assumes that the first 3 columns are: Year   Lat   Lon and omits these from creating the Total")
  }else{
    model_out_df$Total <-rowSums(model_out_df[, 4:(dim(model_out_df)[2]) ],na.rm = TRUE)
    print("No Total column found, creating Total. note, the function assumes that the first 3 columns are: Year   Lat   Lon and omits these from creating the Total")
  }
  return(model_out_df)
}

# to remove the variability caused by the climate:
create_rolling_means <- function(var_to_plot, k = 10, total_only = FALSE,align = "right", rm_firstX_years = NULL){
  # var-to-plot variable that needs plotting, coming from above read-in procedures
  # k = bin width ( default 10)
  # align = rolling means from the start or end of the dataseries
  #rm_firstXyears the number of years which should be removed from the timeseries for rolling mean generation.
  #This is useful when one wants to omit the first 30 "spinup " years, so that they don't screw up the running mean calculations.
  # some years provide more than 30 years, so this must be a model-specific number that is passed into the function
  if(!is.null(rm_firstX_years)){
    maxlength <- dim(var_to_plot@data)[1]
    var_to_plot@data <- var_to_plot@data[rm_firstX_years:maxlength,]
  }
  if(total_only ==TRUE){
    var_roll_means <- var_to_plot@data
    var_roll_means_out <- as.data.frame(zoo::rollmean(var_roll_means$Total, k = k, fill = NA,align=align)) 
    #var_roll_means_tmp <- var_roll_means$Total[!is.na(var_roll_means$Total) & !is.null(var_roll_means$Total)]
    #var_to_plot@data$Total <- TTR::EMA(var_roll_means_tmp, n = k)
    var_to_plot@data$Total <- var_roll_means_out #`zoo::rollmean(var_roll_means$Total, k = k, fill = NA, align = align)`
    # p <- print(plotTemporal(var_to_plot, layers = c("Total"),sizes = c(1),
    #                         size.by = "Layer")) # sizes -> lines made thicker
  }else{
    var_roll_means <- var_to_plot@data
    end = dim(var_roll_means)[2]-1 # omit total. recalculate that one in update_total
    
    var_roll_means_out <- as.data.frame(zoo::rollmean(var_roll_means[,4:end], k = k, fill = NA,align=align)) 
    #reintroduce NAs as 0
    #var_roll_means_out[is.na(var_roll_means_out)] <- 0
    
    var_to_plot@data[,4:end] <- var_roll_means_out
    #determine layers (ls) for plotting- exclude Total:
    #https://stackoverflow.com/questions/40885360/remove-entries-from-string-vector-containing-specific-characters-in-r
    omit_strings = c("Year","Lat","Lon","Total")
    ns = names(var_to_plot@data)
    ls = ns[!ns %in% grep(paste0(omit_strings, collapse = "|"), ns, value = T)]
    var_to_plot@data <- update_total(var_to_plot@data)
    
    #p <- print(plotTemporal(var_to_plot, layers = ls,sizes = rep(1,length(ls)),
    #                        size.by = "Layer")) # sizes -> lines made thicker
  }
  var_to_plot@year.aggregate.method <- paste0("Rolling means: zoo::rollmean(var_roll_means$Total, k = k, fill = NA,align=align)) , which k=",k,",align=",align)
  var_to_plot@id <- paste0(var_to_plot@id,".rollmean.k=",k)
  return(var_to_plot)
  
}

# to calculate a sensible mean landscape value for BiomeEP
# deprecated
create_BiomeEP_means <- function(out){
  
  #because BiomeEP has regular disturbance intervals, I sample randomly across this to create the landscape mean,
  # as the random sample aspect.
  set.seed(143)
  
  lower_tmp  <- sample.int(390,60) # random lower bound of sampling range between 60 and 450 (=upper bound for BiomeEP timeseries)
  lower      = lower_tmp + 60
  upper_tmp  <- sample.int(60, 60) # random width of the sampling range
  upper      = upper_tmp + lower
  
  if(sum(upper >= 450)>1){ # resample if any upper value above 500
    idx <- which(upper >= 450)
    l <- length(lower[idx])
    maxrange = max(lower[idx])
    upper[idx] <- lower[idx] + sample.int(450 - maxrange+l, l) 
  } 
  if(sum(upper >= 450)>1){ # resample if any upper value above 500
    idx <- which(upper >= 450)
    l <- length(lower[idx])
    maxrange = max(lower[idx])
    upper[idx] <- lower[idx] + sample.int(450 - maxrange, l) 
  } 
  if(sum(upper >= 450)>1){ # resample if any upper value above 500
    idx <- which(upper >= 450)
    l <- length(lower[idx])
    maxrange = max(lower[idx])
    upper[idx] <- lower[idx] + sample.int(450 - maxrange, l) 
  } 
  #hopefully by now all values > 450 are taken care of
  
  #create means and record them 
  means <- ssize <- rep(NA,length(lower))
  for(ii  in 1:length(lower)){
    ssize[ii]      = length(seq(lower[ii],upper[ii]))
    means[ii]      = mean(out@data[lower[ii]:upper[ii],]$Total)
    
  }
  wm <- weighted.mean(means,w=ssize)
  return(wm)
  
}

# for self-thinning:
# model_name and site only needed for plotting.
identify_thinning_period <- function(agcwood,nstem,cmort=NULL,add_plot=FALSE,manadjust,cmort_plot=NULL,model_name, site){
  
  if(!is.na(manadjust$upper)){ #if thinning points were manually adjusted, just get them here now.
    # if they were not, then continue with the below options 
    # (whith are to either where a thinning mortality helps to determine a useful set of thinning line indeces,
    # or where the upper and lower boundary of a straight line determines something useful 
    # (probably not, but was helpful for determining the manually adjusted indeces))
    
    df <- data.frame(indiv_mass = log10(agcwood/nstem),nstem = log10(nstem),dnstem= NA,thin=NA,time=seq(1:length(nstem)))
    
    #here, index_* is the year after year 33 ( the time we start universally countin regrowth in this case)
    idx_lower <-  manadjust$lower
    idx_upper <-  manadjust$upper
    
    if(add_plot){
      points(df[idx_lower:idx_upper,c("nstem","indiv_mass") ],col="red",ylim=ylim_set,xlim=xlim_set,cex=point_cex)
      #mtext("method 3)",side=1,line=-1)
    }
    
  }else{
    if(!is.null(cmort)){# this means that the FATES LPJG, EDv3 or CABLEPOP model has been called. FATES in the tropics has a very strong slope that is not 
      # the result of self-thinning, but of other mortality mechanisms. 
      # The above filtering method can therefore not be reliably applied to distinguish between the self-thinning slope and the other mortality slope.
      # But we have the cstarv mortality mechanism, so with that column, we can create another dimension by which we select for the self-thinning curve.
      df <- data.frame(indiv_mass = log10(agcwood/nstem),nstem = log10(nstem), cmort = cmort,dnstem= NA,thin=NA,time=seq(1:length(nstem)) )
      # here, "time" is also the year after year 33. ( i.a. time=1, is simulation year 31)
      
      #do some cleaning, as  0's in some models have screwed the logging up, that is ok, 0s just mean no individuals, 
      #that happens at the beginning of the regrowth period.
      if(sum(is.nan(df$indiv_mass))>0){
        df[which(is.nan(df$indiv_mass)),]$nstem <- NA
        df[which(is.nan(df$indiv_mass)),]$indiv_mass <- NA
      }
      
      # must add another variable for evaluation for FATES:
      cstarv_med = quantile(cmort, c(.95))#c(.90))
      
      #first remove all non-thinning "slopes"
      #retain start df for plotting:
      df_start <- df
      df <- df[which(df$cmort >=cstarv_med),]
      
      
      if(!is.null(cmort_plot)){
        plot(df_start$cmort, ylab= "mortality rate (%/yr)", xlab="year", main=paste(model_name,site))
        #points(df[which(df$cmort >=cstarv_med),]$time,df[which(df$cmort >=cstarv_med),]$cmort,col="red")
        points(df$time,df$cmort,col="red")
        legend(legend=c("cmort","mort_95thperc"),fill=c(1,"red"),"topright",cex=0.6)
      }
      
      #to make compatible with the other selection options at return:
      #"find" upper edge and lower edge of the self-thinning line:
      idx_upper <- dim(df)[1]
      idx_lower <- 1
      # preps for plotting only:
      upper <- df[idx_upper,]
      lower <- df[idx_lower,]
      
      if(add_plot){
        
        #plot(df_start$nstem,df_start$indiv_mass)
        points(upper[,c("nstem","indiv_mass")],col="green",ylim=ylim_set,xlim=xlim_set,cex=point_cex)
        points(lower[,c("nstem","indiv_mass")],col="green",ylim=ylim_set,xlim=xlim_set,cex=point_cex)
        #from this, "triangulate" the values in that range that belong into the self-thinning line
        points(df[which(df$nstem >= upper$nstem & df$indiv_mass >= lower$indiv_mass),c("nstem","indiv_mass") ],col="red",ylim=ylim_set,xlim=xlim_set,cex=point_cex)
        #testing
        #points(df[which(df$time < upper$time & df$time > lower$time ),c("nstem","indiv_mass")],ylim=ylim_set,xlim=xlim_set,cex=point_cex,col="purple")
      }
      
    }else{# all other models:
      df <- data.frame(indiv_mass = log10(agcwood/nstem),nstem = log10(nstem),dnstem= NA,thin=NA,time=seq(1:length(nstem)))
      
      #do some cleaning, as  0's in some models have screwed the logging up, that is ok, 
      #0s just mean no individuals, that happens at the beginning of the regrowth period.
      if(sum(is.nan(df$indiv_mass))>0){
        df[which(is.nan(df$indiv_mass)),]$nstem <- NA
        df[which(is.nan(df$indiv_mass)),]$indiv_mass <- NA
      }
      
      if(sum(is.infinite(df$indiv_mass))>0){
        df[which(is.infinite(df$indiv_mass)),]$indiv_mass <- NA
        df[which(is.infinite(df$nstem)),]$nstem <- NA
      }
      
      #find upper edge and lower edge of the self-thinning line:
      upper <- df[which(df$indiv_mass == max(df$indiv_mass,na.rm = TRUE)),]
      lower <- df[which(df$nstem == max(df$nstem,na.rm = TRUE)),]
      
      #find index for upper edge and lower edge of the self-thining line: 
      idx_upper <-  which(df$indiv_mass == max(df$indiv_mass,na.rm = TRUE))
      idx_lower <-  min(which(df$nstem == max(df$nstem,na.rm = TRUE))) # hack to account for some models having multiple values for this == test
      
      if(add_plot){
        
        # plot(df$nstem,df$indiv_mass)
        points(upper[,c("nstem","indiv_mass")],col="red",ylim=ylim_set,xlim=xlim_set,cex=point_cex)
        points(lower[,c("nstem","indiv_mass")],col="red",ylim=ylim_set,xlim=xlim_set,cex=point_cex)
        #from this, "triangulate" the values in that range that belong into the self-thinning line
        #points(df[which(df$nstem >= upper$nstem & df$indiv_mass >= lower$indiv_mass),c("nstem","indiv_mass") ],col="green",ylim=ylim_set,xlim=xlim_set,cex=point_cex)
        #testing - this does not work for EDv3 FIN and BIA
        points(df[c(idx_lower:idx_upper),c("nstem","indiv_mass")],ylim=ylim_set,xlim=xlim_set,cex=point_cex,col="red")
        #mtext("method 3)",side=1,line=-1)
      }
      
    }
  }
  
  thinning_points <- df[c(idx_lower:idx_upper),c("nstem","indiv_mass","time") ]
  return(thinning_points)
  
}

# rates a %rate object for stemmort and cmort, that can be used within the DGVMTools plotting capabilities:
# mort_in  mortality flux in carbon or number of stems
# total total amount of carbon or number of stems
# var_in neeed to update metadata, either cmort_rate or stemmort_rate
create_mort_rate <- function(mort_in,total,var_in){
  
  rate.final <- mort_in
  #to be save only calculate mortrate on total:
  rate.final@data <- mort_in@data[,1:3]
  rate.final@data$Total <- mort_in@data$Total/total@data$Total*100
  
  #adjust metadata
  rate.final@quant <- get_quantity(var_in)
  
  return(rate.final)
  
}

# calculate_ensemble_mean_median()
# used for deriving mean carbon mass turnover etc values at equilibrium, to show 
# in the manuscript that demographic models increase in mortality, but not as much as the dummy model

# ns= list of model names in the ensemble
# var variable for which mean and median should be calculated
# n = length of time period before simulation end, for which the simulation output should be obtained. multiple of 30 to account for forcing variability.
# returns a list of two data-frames, element 1 is ambient and emelent 2 is elevated D-BEN runs
# site = simulation site
calculate_ensemble_mean_median <- function(var,site, ns, n=60, run="P0"){
  
  #initialise collection dataframes:
  collect_a <- collect_e <- as.data.frame(matrix(NA,ncol=length(ns),nrow=n))
  names(collect_a) <- ns
  names(collect_e) <- ns
  
  #median_col is a column index needed to be removed when calculating the row (model ensemble) mean of a given timestep:
  median_col =length(ns)+1
  
  #retrieve model output for ambient and elevated simulations
  for(model_name in ns){
    model_out <-  get_model_output(site=site, run = run, var=var,model_name=model_name, co2_levels = "412ppm")
    end <- length(model_out@data$Total)
    model_out <- model_out@data$Total[(end-(n-1)):end] # flexible end-point because models finish at different times, depending on when they were in equilibrium
    collect_a[model_name] <- model_out
    
    model_out <-  get_model_output(site=site, run = run, var=var,model_name=model_name, co2_levels = "562ppm")
    end <- length(model_out@data$Total)
    model_out <- model_out@data$Total[(end-(n-1)):end]
    collect_e[model_name] <- model_out
  }
  
  #create median and mean  for "var" of all selected models (ns):
  collect_a$median <-  apply(collect_a, 1,  median, na.rm = TRUE)
  collect_e$median <-  apply(collect_e, 1,  median, na.rm = TRUE)
  if(median_col == 2){# means that median and mean are the same, as only one model is tested
    collect_a$mean   <-  collect_a$median
    collect_e$mean   <-  collect_e$median
  }else{   # more than one model. this switch was just added as convenience to enable the upstream script to work with one. model. 
    # the analysis of mean and median here are not really relevant for only one model.
    collect_a$mean   <-  apply(collect_a[,-c(median_col)], 1,  mean, na.rm = TRUE)
    collect_e$mean   <-  apply(collect_e[,-c(median_col)], 1,  mean, na.rm = TRUE)
  }
  return(list(collect_a,collect_e))
}

# plotting ggplots
# convenience function for grid.arrange
# passes back ggplot object for each model_output (input)
# helper function to plot time-series
create_gobj <-function(model_out, y.lim = NULL,layers=NULL){
  if(!is.null(layers)){ # enhance title, so that the reader understands what layers are currently plotted
    p_out <- plotTemporal(model_out, legend.position = "none", layers=layers,title = paste(model_out@source@name," \n", layers), y.lim )
  }else{# if everything is plotted:
    p_out <- plotTemporal(model_out, legend.position = "none", layers=layers,title = model_out@source@name, y.lim )
    
  }
  
  return(p_out)
}

calculate_vegetation_percentage_fraction <- function(data){
  
  #Test 1: must contain Total column:
  sum(names(data@data) %in% c( "Total") == 1)
  
  #Test 2: must contain other columns:
  PFTs <- names(data@data)[!names(data@data) %in% c("Year", "Lat", "Lon", "Total")]
  df1 <- data@data %>%
    mutate(across(c(PFTs,"Total"), ~ ./data@data$Total * 100))
  
  data@data <- df1
  
  data@quant@name <- paste0("Vegetation fraction")
  data@quant@units <- "% total"
  
  return(data)
}

test_Wbudget_tolerance_threshold <- function(WBgrowth_in, cmort_in, cwood_in, model_name, site,time_range = 31:410) {
  # budget test across the whole simulation period from disturbance to the prescribed 450 years. 
  # equilibrium  may not be reached by then, and some models have run their models for longer. 
  # I use the mean of the last 30 years of the "official" simulations in this script.
  # must start from at least 32, because of (atmospheric) planting offset in ORCHIDEE
  
  # Subset the data to the desired time range
  WBgrowth <- WBgrowth_in[time_range]
  cmort    <- cmort_in[time_range]
  
  if(model_name=="BiomeE-Standalone" ){
    WBgrowth <-   c(cwood_in[1],WBgrowth[1:(length(WBgrowth)-1)])
    cmort <- c(0.0,cmort[1:(length(cmort)-1)])
  }
  
  cwood    <- cwood_in[time_range]
  if(model_name=="ORCHIDEE"){
    #We take carbon from the atmosphere on the first time step. 
    #This carbon is still there the second, third, ...
    #time step so this should be implemented as a bias correction on
    #the entire time series of the pool cwood
    cwood <- cwood - (cwood[1]-(WBgrowth[1]-cmort[1]))
  }
  if(model_name=="JULES-RED"){
    #  cwood <- cwood - cwood[1]
  }
  
  lgth <- length(cwood) # for equilibrium period 
  
  # Calculate flux-derived cwood pool
  flux_derived_cwood <-   cumsum(WBgrowth-cmort)
  
  # Calculate the mean of the deviation between actual cwood and flux-derived cwood
  deviation <- mean(abs(cwood- flux_derived_cwood))
  #deviation <- mean(abs(cwood[time_range] - flux_derived_cwood))
  
  # Calculate the mean equilibrium of cwood (mean over the period identified as equilibrium period)
  lower <- eq_values[[site]][which(eq_values[[site]]$model==model_name),]$lower-32
  # Address issue with >equilibrium period being omitted prior to this function
  upper <- min(eq_values[[site]][which(eq_values[[site]]$model==model_name),]$upper-32,lgth)
  mean_equilibrium_cwood <- mean(cwood[(lgth-30):lgth]) 
  #mean_equilibrium_cwood <- mean(cwood[lgth-30:lgth]) # do not omit NAs. if NAs (error message), then something wrong.
  # Set the tolerance threshold (10% of mean equilibrium cwood)
  tolerance_threshold <- 0.10 * mean_equilibrium_cwood
  
  # Test whether the deviation exceeds the tolerance threshold
  exceeds_tolerance <- deviation > tolerance_threshold
  
  if (deviation > tolerance_threshold) {
    print(paste ("Cumulative Carbon deviation ",deviation ,"exceeds tolerance ",tolerance_threshold, "at ",site))
    plotstring =">"
    passfail="threshold FAILED:"
  } else {
    print("Cumulative Carbon deviation is within the tolerance.")
    plotstring ="<"
    passfail="threshold PASSED:"
  }
  
  mname = model_name
  if(model_name=="BiomeE-Standalone"){
    mname = "BiomeE"
  }
  plot(cwood,lwd=1.5, type="l",ylim=c(0,50), xlim=c(0,450),main =  paste(model_name,site),cex.main=2, yaxt="n", yaxt="n",xaxt="n")
  lines(flux_derived_cwood,lwd=1.5, col=get_model_colour(model_name))
  
  #plotting aesthetics:
  axis(1,  line = -0.6, lwd = 0, cex.axis = 1.2)
  axis(1, lwd = 1, cex.axis = 0.9,labels = FALSE,tck=0.02)
  axis(1, tck = 0.01 ,cex.axis = 0.8,labels = FALSE,tck=0.01) 
  
  if(site=="FIN"){
    legend("topright",fill=c(1, col=get_model_colour(model_name)),legend=c("cwood","WBgrowth-cmort"),cex=0.8)
    mtext( paste0(passfail," \n ",round(deviation,digits=2) ,"(Mean deviation) ", plotstring ,round(tolerance_threshold,digits=2),"( 10% eq C)"),side=3,cex=0.8,line=-4.7) 
    axis(2,  line = -0.6, lwd = 0, cex.axis = 0.9)
    axis(2, lwd = 1, cex.axis = 0.9,labels = FALSE,tck=-0.02)
    #axis(1, tck = 0.01 ,cex.axis = 0.8,labels = FALSE,tck=0.01) 
    
  }else{
    mtext( paste0(passfail," \n ",round(deviation,digits=2) ,"(Mean deviation) ", plotstring ,round(tolerance_threshold,digits=2),"( 10% eq C)"),side=3,cex=0.8,line=-4.7) 
    #axis(1,  line = -0.6, lwd = 0, cex.axis = 0.9)
    #axis(1, lwd = 1, cex.axis = 0.9,labels = FALSE,tck=-0.02)
    axis(2, lwd = 1, cex.axis = 0.9,labels = FALSE,tck=0.02)
    #axis(2, tck = 0.01 ,cex.axis = 0.8,labels = FALSE,tck=0.01) 
    
  }
  
}

test_Wbudget_tolerance_threshold_ORCHIDEE <- function(WBgrowth_in, cmort_in, cwood_in, model_name, site,time_range = 31:300) {
  
  # subset the data to the desired time range
  WBgrowth <- WBgrowth_in[time_range]
  cmort <- cmort_in[time_range]
  cwood <- cwood_in[time_range]
  
  if(model_name=="ORCHIDEE"){
    # we take carbon from the atmosphere on the first time step. This carbon is still there the second, third, ... time step so this should be implemented as a bias correction on the entire time series of the pool cwood (not on the flux WBgrowth).
    #cwood[1] <- cwood - (cwood[1]-(WBgrowth[1]+cmort[1]))
    
    cwood <- cwood - (cwood[1]-(WBgrowth[1]-cmort[1]))
  }
  
  # Calculate flux-derived cwood pool
  flux_derived_cwood <-   cumsum(WBgrowth-cmort)
  
  # Calculate the mean of the deviation between actual cwood and flux-derived cwood
  deviation <- mean(abs(cwood - flux_derived_cwood))
  
  # Calculate the mean equilibrium of cwood (mean over the period you identify as equilibrium period)
  mean_equilibrium_cwood <- mean(cwood[250:200])
  
  # Set the tolerance threshold (10% of mean equilibrium cwood)
  tolerance_threshold <- 0.10 * mean_equilibrium_cwood
  
  # Test whether the deviation exceeds the tolerance threshold
  exceeds_tolerance <- deviation > tolerance_threshold
  
  if (deviation > tolerance_threshold) {
    print(paste ("Cumulative Carbon deviation ",deviation ,"exceeds tolerance ",tolerance_threshold, "at ",site))
    plotstring =">"
    passfail="threshold failed:"
  } else {
    print("Cumulative Carbon deviation is within the tolerance.")
    plotstring ="<"
    passfail="threshold passed:"
  }
  
  plot(cwood, type="l",ylim=c(0,50), main =  "ORCHIDEE")
  lines(flux_derived_cwood, col=get_model_colour("ORCHIDEE"))
  mtext(side=2,"cwood, KgC m-2",line=2,cex=0.7)
  legend("topright",fill=c(1, col=get_model_colour("ORCHIDEE")),legend=c("cwood","WBgrowth-cwood"))
  
  plot.new()
  mtext( paste0(passfail," \n ",round(deviation,digits=2) ,"(Mean deviation) ", plotstring ,round(tolerance_threshold,digits=2),"( 10% eq C)"),side=3,cex=0.8,line=-1) 
}

#####
# useful for another analysis:
# create and plot model output summary stats for stand structure across a specified timeperiod:

add_stand_structure_benchmarks_summary <- function(model_out,year_lower,year_upper){
  
  #subset for certain years:
  model_out@data <-  model_out@data %>%filter(Year >year_lower & Year < year_upper)
  
  # Select the size class columns (excluding Year, Lat, Lon, and Total)
  size_class_cols <- colnames(model_out@data)[!(colnames(model_out@data) %in% c("Year", "Lat", "Lon", "Total"))]
  
  # Create a dataframe with median values
  median_df <- model_out@data %>%
    summarise(across(all_of(size_class_cols), median, na.rm = TRUE)) %>%
    mutate(Year = NA, Lat = NA, Lon = NA, Total = NA) %>%  # Add the non-size-class columns back as NA
    select(colnames(model_out@data))  # Keep the column order same as the input
  
  # Create a dataframe with 10th percentile values
  p10_df <- model_out@data %>%
    summarise(across(all_of(size_class_cols), ~ quantile(., 0.10, na.rm = TRUE))) %>%
    mutate(Year = NA, Lat = NA, Lon = NA, Total = NA) %>%
    select(colnames(model_out@data))  # Keep the column order same as the input
  
  # Create a dataframe with 90th percentile values
  p90_df <- model_out@data %>%
    summarise(across(all_of(size_class_cols), ~ quantile(., 0.90, na.rm = TRUE))) %>%
    mutate(Year = NA, Lat = NA, Lon = NA, Total = NA) %>%
    select(colnames(model_out@data))  # Keep the column order same as the input
  
  points(dbh_classes_plot, 
         as.numeric(median_df[,c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                  "<90" ,  "<100",  "<150" , "<200",  ">=200")]),  
         col= "red")
  arrows(dbh_classes_plot,
         as.numeric(p10_df[1,c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                                "<90" ,  "<100",  "<150" , "<200",  ">=200")]),
         dbh_classes_plot,
         as.numeric(p90_df[,c( "<1" ,   "<5" ,   "<10", "<15",  "<20" ,  "<30",   "<40"  , "<50" ,  "<60" ,  "<70" ,  "<80" , 
                               "<90" ,  "<100",  "<150" , "<200",  ">=200")]),
         length = 0.00, angle = 90, code = 3)
  
}

benchmark_WBgrowth_flux <- function(model_name = model_name,site = site, var= "WBgrowth",eq_values_in=eq_values){
  #eq_values_in list of dataframes with equilibrium /mature forest values, which are used for subsetting for this benchmark
  
  # model output:
  if(model_name=="CABLE-POP" | model_name=="LPJ-GUESS"){
    run="P0"
  }
  else{
    run="P0"
  }
  WBgrowth    <- get_model_output(var = var,run=run,site = site,model_name = model_name,co2_levels = "412ppm") 
  WBgrowth    <- omit_equilibrium_phase(WBgrowth,model_name = model_name,site=site,selection_by_year = TRUE)
  
  if(model_name=="SEIB-DGVM"){
    #allow for 10-year smoothing in SEIB-DGVM:
    #H.Sato 29.10.2024: "This issue likely arises from the different timings used to calculate growth flux and mortality flux.
    #Since this is an individual-based model, the death of a large individual has a pronounced effect.
    #By using a 10-year moving average, these negative values are expected to disappear."
    WBgrowth <- create_rolling_means(WBgrowth,k=10,total_only=TRUE, align="left")
  }
  
  # extract the equilibrium phase, and use those years for the boxplot
  lower <- eq_values[[site]][which(eq_values[[site]]$model==model_name),]$lower-30# omit early equilibrium phase ( first 30 years)
  upper <- eq_values[[site]][which(eq_values[[site]]$model==model_name),]$upper-30# omit early equilibrium phase
  
  # apply 30 year smoothing to reduce repeated years climate-effect on the variability:
  #WBgrowth_smoothed <- zoo::rollmean(WBgrowth@data$Total[lower:upper], k=30, align="left")
  
  #WBgrowth_out <- na.omit(WBgrowth_smoothed) # rm NAs that emerged from smoothing on the right side.
  
  # collect turnover values, used in boxplot
  # collect_turnover_BCI[[model_name]] <- cwood@data$Total[lower:upper]/cmort@data$Total[lower:upper]
  
  # collect equilibrium content WBgrowth values
  WBgrowth_out<- WBgrowth@data$Total[lower:upper]
  
  return(WBgrowth_out)
  
}

omit_sizeclasses_below10dbh <- function(input){
  output <- input 
  
  # model_name_in <- strsplit(input@source@name,split=" ")[[1]][1]
  # site_in <- strsplit(input@source@name,split=" ")[[1]][2]
  # run_in <- strsplit(input@source@name,split=" ")[[1]][4]
  
  # test for which sizeclasses exist for this output:
  idx <- unlist(which( names(input@data) %in% c( "<1"  ,"<5" ,"<10" )))
  
  # remove small sizeclasses that exist
  tmp <- input@data[,-idx, with = FALSE]
  # now that some sizeclasses are missing, the "total" column has to be updated
  out_tmp  <- update_total(tmp) 
  
  # put the updated data-table back into @data object of model output:
  output@data <- out_tmp
  
  return(output)
}
