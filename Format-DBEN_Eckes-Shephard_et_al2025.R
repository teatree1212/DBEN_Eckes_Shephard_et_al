# Format file to handle D-BEN data (model output and observations) using DGVMTools.
# author: Annemarie Eckes-Shephard
# date: 2022-2024
# https://github.com/MagicForrest/DGVMTools.


Layers_DBEN <- list(
  
  #PFT1 - Pinus sylvestris
  new("Layer",
      id =   "BINE",
      name = "Boreal/Temperate Shade-Intolerant Needleleaved Evergreen Tree",
      colour = "dodgerblue3",
      properties = list("PFT",
                        growth.form =  "Tree",
                        leaf.form =  "Needleleaved",
                        phenology =  "Evergreen",
                        #climate.zone =  "Boreal", 
                        shade.tolerance = "BNE",
                        land.cover ="Natural")
  ),
  #PFT2 - Picea abies
  new("Layer",
      id =  "NE",
      name = "Boreal/Temperate Needleleaved Evergreen Tree",
      colour ="darkblue",
      properties = list( "PFT",
                         growth.form =  "Tree",
                         leaf.form = "Needleleaved",
                         phenology = "Evergreen",
                         #climate.zone = "Boreal",
                         shade.tolerance = "None",
                         land.cover = "Natural")
  ),
  
  #PFT3 - betula pendula
  new("Layer",
      id =  "IBS",
      name = "Shade-intolerant Broadleaved Summergreen Tree",
      colour ="chartreuse",
      properties = list( "PFT",
                         growth.form =  "Tree",
                         leaf.form = "Broadleaved",
                         phenology = "Summergreen",
                         # climate.zone = "Temperate",
                         shade.tolerance = "None",
                         land.cover =  "Natural")
  ),
  
  
  #!/////////////////////////////////////////////////////////////////////////////////
  #  ! DBEN: some PFTs have to be turned off, we only allow certain ones to run at
  #! the selected sites.
  #! Selecting only a subset of PFTs for these simulations, in the following order:
  #  !Pinus sylvestris (PFT1): BINE, (FIN,BIA)
  #!Picea abies (PFT2): BNE, (FIN)
  #!betula pendula (PFT3): IBS  (FIN,BIA)
  #!Carpinus betulus or Tilia
  #!cordata (intermediate) shade tolerant broadleaf deciduous (PFT4): TeBS (BIA)
  ##!Tropical broadleaf evergreen shade intolerant (PFT5): TrIBE (BCI)
  #!Tropical broadleaf evergreen shade tolerant (PFT6): TrBE (BCI)
  #!Tropical broadleaf deciduous (PFT7): TrBR (BCI)
  #!Grasses, to be lumped together in post-processing (PFT8): C3 C4 (FIN,BCI,BIA)
  #!/////////////////////////////////////////////////////////////////////////////////
  
  
  
  #PFT4 - Carpinus betulus or Tilia cordata
  new("Layer",
      id =  "TeBS",
      name = "(intermediate) shade tolerant broadleaf deciduous",
      colour ="purple",
      properties = list( "PFT",
                         growth.form =  "Tree",
                         leaf.form = "Broadleaved",
                         phenology = "Summergreen",
                         climate.zone = "Temperate",
                         shade.tolerance = "None",
                         land.cover =  "Natural")
  ),
  
  #PFT5 Tropical broadleaf evergreen shade intolerant (PFT5) 
  new("Layer",
      id =  "TrIBE",
      name = "Tropical Shade-intolerant Broadleaved Evergreen Tree",
      colour ="orchid",
      properties = list( "PFT",
                         growth.form =  "Tree",
                         leaf.form = "Broadleaved",
                         phenology = "Evergreen",
                         climate.zone = "Tropical",
                         shade.tolerance = "TrBE",
                         land.cover =  "Natural")
  ),
  
  #PFT6 Tropical broadleaf evergreen shade tolerant  
  new("Layer",
      id =  "TrBE",
      name = "Tropical Broadleaved Evergreen Tree",
      colour ="orchid4",
      properties = list( "PFT",
                         growth.form =  "Tree",
                         leaf.form = "Broadleaved",
                         phenology = "Evergreen",
                         climate.zone = "Tropical",
                         shade.tolerance = "None",
                         land.cover =  "Natural")
  ),
  
  # PFT7 Tropical broadleaf deciduous (PFT7) 
  new("Layer",
      id =  "TrBR",
      name = "Tropical Broadleaved Raingreen Tree",
      colour ="palevioletred",
      properties = list( "PFT",
                         growth.form =  "Tree",
                         leaf.form = "Broadleaved",
                         phenology = "Raingreen",
                         climate.zone = "Tropical",
                         shade.tolerance = "None",
                         land.cover =  "Natural")
  ),
  
  # PFT8 Grasses (C3+C4, whatever is present)
  new("Layer",
      id =  "Grasses",
      name = "C3 and/or C4 grasses",
      colour ="dark grey",
      properties = list( "PFT",
                         growth.form =  "Grass",
                         #leaf.form = "Broadleaved",
                         #phenology = "Raingreen",
                         #climate.zone = "Tropical",
                         shade.tolerance = "None",
                         land.cover =  "Natural")
  ),
  new("Layer",
      id = "Total",
      name = "Total",
      colour = "black",
      properties = list(type = "Sum",
                        land.cover = "All")
  )
  
  
  
)

# add size class layers: added separately and in a loop to be able to change colours easily
sc <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", "<250","300+")
#[TODO] added this 300 sizeclass for FATES, but this is not actually correct.
colfunc <- colorRampPalette(c("green", "brown"))
colours <- colfunc(length(sc))
sc_layers <- list()
for(i in 1:length(sc)){
  sc_layers[[i]] <- new("Layer",
                        id = sc[i],
                        name = "sizeclass_by_diameter_at_breastheight (cm)",
                        colour = colours[i],
                        properties = list("Sizeclass",
                                          growth.form =  "Tree",
                                          leaf.form = "Any",
                                          phenology = "Any",
                                          climate.zone = "Any",
                                          shade.tolerance = "Any",
                                          land.cover =  "Natural")
  )
}


Layers_DBEN <- append(Layers_DBEN,sc_layers)

##add ORCHIDEE_specific sizeclass layers:
sc_layers_ORCHIDEE <- list()
sc_ORCHIDEE <- c("1","2","3")
cols <- c("light green","dark green","brown")
for(i in 1:3){
  sc_layers_ORCHIDEE[[i]] <- new("Layer",
                                 id = sc_ORCHIDEE[i],
                                 name = "sizeclass_by_diameter_at_breastheight (cm)",
                                 colour = cols[i],
                                 properties = list("Sizeclass_ORCHIDEE",
                                                   growth.form =  "Tree",
                                                   leaf.form = "Any",
                                                   phenology = "Any",
                                                   climate.zone = "Any",
                                                   shade.tolerance = "Any",
                                                   land.cover =  "Natural")
  )
}
Layers_DBEN <- append(Layers_DBEN,sc_layers_ORCHIDEE)

# add layers for FATES:
# add size class layers: added separately and in a loop to be able to change colours easily
sc <- c("<1", "<5", "<10", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", "<250","<300","300+")
#[TODO] added this 300 sizeclass for FATES, but this is not actually correct.
colfunc <- colorRampPalette(c("purple","green","red","orange","blue","brown"))
colours <- colfunc(length(sc))
sc_layers_FATES <- list()
for(i in 1:length(sc)){
  sc_layers_FATES[[i]] <- new("Layer",
                              id = sc[i],
                              name = "sizeclass_by_diameter_at_breastheight (cm)",
                              colour = colours[i],
                              properties = list("Sizeclass_FATES",
                                                growth.form =  "Tree",
                                                leaf.form = "Any",
                                                phenology = "Any",
                                                climate.zone = "Any",
                                                shade.tolerance = "Any",
                                                land.cover =  "Natural")
  )
}

Layers_DBEN <- append(Layers_DBEN,sc_layers_FATES)

# add layers for new runs for Paper1, updated protocol:
# add size class layers: added separately and in a loop to be able to change colours easily
sc <- c("<1", "<5", "<10", "<15","<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200")
colfunc <- colorRampPalette(c("purple","green","red","orange","blue" ,"brown"))
colours <- colfunc(length(sc))
sc_layers_Paper1 <- list()
for(i in 1:length(sc)){
  sc_layers_Paper1[[i]] <- new("Layer",
                               id = sc[i],
                               name = "sizeclass_by_diameter_at_breastheight (cm)",
                               colour = colours[i],
                               properties = list("sc_layers_Paper1",
                                                 growth.form =  "Tree",
                                                 leaf.form = "Any",
                                                 phenology = "Any",
                                                 climate.zone = "Any",
                                                 shade.tolerance = "Any",
                                                 land.cover =  "Natural")
  )
}

Layers_DBEN <- append(Layers_DBEN,sc_layers_Paper1)

for(Layer in Layers_DBEN){
  print(Layer)
}


# add layers for new runs for Paper1, updated protocol:
# add size class layers: added separately and in a loop to be able to change colours easily
sc <- c("sdlg<1", "sdlg<5", "saplg<10", "saplg<15")
#colfunc <- colorRampPalette(c("green","orange","darkgreen" "brown"))
colours <- c("green","orange","darkgreen", "brown")
sc_layers_Paper1_estab <- list()
for(i in 1:length(sc)){
  sc_layers_Paper1_estab[[i]] <- new("Layer",
                                     id = sc[i],
                                     name = "sizeclass_by_diameter_at_breastheight (cm)",
                                     colour = colours[i],
                                     properties = list("Sizeclass_estab",
                                                       growth.form =  "Tree",
                                                       leaf.form = "Any",
                                                       phenology = "Any",
                                                       climate.zone = "Any",
                                                       shade.tolerance = "Any",
                                                       land.cover =  "Natural")
  )
}

Layers_DBEN <- append(Layers_DBEN,sc_layers_Paper1_estab)

for(Layer in Layers_DBEN){
  print(Layer)
}



quantities_DBEN <- list(
  
  # Just a couple of example quantities
  
  new("Quantity",
      id = "cveg",
      name = "Vegetation Carbon Mass  by PFT",
      units = "kgC m-2",
      colours = function(n) rev(viridis::viridis(n)),
      format = c("DBEN"),
      standard_name = c("Carbon_mass_in_vegetation_by_PFT")
  ),
  
  new("Quantity",
      id = "albedo",
      name = "surface reflectiveness of PFTs",
      units = "-",
      colours = function(n) rev(viridis::viridis(n)),
      format = c("DBEN"),
      standard_name = c("surface_reflectiveness_of_PFTs")
  ),
  
  new("Quantity",
      id = "z0",
      name = "roughness length (momentum)",
      units = "m",
      colours = function(n) rev(viridis::viridis(n)),
      format = c("DBEN"),
      standard_name = c("roughness length")
  ),
  
  new("Quantity",
      id = "cwood",
      name = "Wood Carbon Mass by PFT",
      units = "kgC m-2",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Carbon_mass_in_wood_by_PFT")
  ),
  
  new("Quantity",
      id = "AGcwood",
      name = "Aboveground woody biomass",
      units = "kgC m-2",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Aboveground_woody_biomas")
  ),
  
  new("Quantity",
      id = "AGcwood_size",
      name = "Aboveground woody carbon by sizeclass",
      units = "kgC m-2",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Aboveground_woody_carbon_by_sizeclass")
  ),
  
  new("Quantity",
      id = "AGB",
      name = "Aboveground total biomass",
      units = "kgC m-2",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Aboveground_total_biomas")
  ),
  
  new("Quantity",
      id = "cwood_size",
      name = "Wood Carbon Mass by size class",
      units = "kgC m-2",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon_mass_in_wood_by_size_class")
  ),
  
  new("Quantity",
      id = "nstem_size",
      name = "Stem number by size class",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Stem_number_by_size_class")
  ),
  
  new("Quantity",
      id = "stemmort_rate",
      name = "Stem mortality rate by size class",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Stem_mortality_rate_by_size_class")
  ),
  
  new("Quantity",
      id = "growth_rate",
      name = "Woody biomass growth rate",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Woody_biomass_growth_rate")
  ),
  
  new("Quantity",
      id = "lai",
      name = "Leaf area index",
      units = "m2 m-2",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("leaf_area_index")
  ),
  
  new("Quantity",
      id = "CA",
      name = "Crown area",
      units = "m2 ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("crown_area")
  ),
  
  new("Quantity",
      id = "BA",
      name = "Basal area",
      units = "m2 ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("basal_area")
  ),
  
  new("Quantity",
      id = "height",
      name = "Mean PFT height in the 95% height percentile",
      units = "m",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("mean_PFT_height")
  ),
  
  
  new("Quantity",
      id = "WBgrowth",
      name = "Woody biomass growth",
      units = "kgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Woody_biomass_growth")
  ),
  
  
  
  new("Quantity",
      id = "BAgrowth",
      name = "Basal area growth",
      units = "m2 ha-1 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Basal_area_growth")
  ),
  
  
  new("Quantity",
      id = "cmort",
      name = "Wood Carbon Mass flux lost by size class",
      units = "kgC m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon_mass_flux_lost_from_life_wood_due_to_mortality_or_other_turnover_processes_mass_in_wood_by_size_class")
  ),
  
  # same as above, but some models (BiomeEP and BiomeES) call it "_size" in filename. To account for that:
  new("Quantity",
      id = "cmort_size",
      name = "Wood Carbon Mass flux lost by size class",
      units = "kgC m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon_mass_flux_lost_from_life_wood_due_to_mortality_or_other_turnover_processes_mass_in_wood_by_size_class")
  ),
  
  #disturbance mortality from BiomeE
  new("Quantity",
      id = "cdstb_size",# id = "cmort_dist",
      name = "Wood Carbon Mass flux lost by size class",
      units = "kgC m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon_mass_flux_lost_from_life_wood_due_to_mortality_or_other_turnover_processes_mass_in_wood_by_size_class")
  ),
  
  new("Quantity",
      id = "cmort_rate",
      name = "% Wood Carbon Mass lost due to mortality, by size class",
      units = "% yr -1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("%Carbon_mass_lost_from_life_wood_due_to_mortality_or_other_turnover_processes_mass_in_wood_by_size_class")
  ),
  
  new("Quantity",
      id = "stemmort",
      name = "Stem number flux lost by size class",
      units = "Count ha-1 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Stem_number_Flux_lost_from_vegetation_due_to_mortality_or_other_turnover_process")
  ),
  
  # same as above, but some models call it "_size" in filename. To account for that:
  new("Quantity",
      id = "stemmort_size",
      name = "Stem number flux lost by size class",
      units = "Count ha-1 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Stem_number_Flux_lost_from_vegetation_due_to_mortality_or_other_turnover_process")
  ),
  
  
  
  new("Quantity",
      id = "gpp",
      name = "gpp",
      units = "kgC m-2 yr-1", # NB in Amsterdam , this unit was s-1, but changed to .yr in next post-processing  iteration
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Gross_primary_production")
  ),
  
  new("Quantity",
      id = "npp",
      name = "npp",
      units = "kgC m-2 yr-1", # NB in Amsterdam , this unit was s-1, but changed to .yr in next post-processing iteration
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Net_primary_production")
  ),
  
  
  new("Quantity",
      id = "nbp", # NB in Amsterdam , this unit was s-1, but changed to .yr in next post-processing iteration
      name = "Net Biospheric Production",
      units = "kgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon_Mass_Flux_out_of_Atmosphere_due_to_Net_Biospheric_Production")
  ),
  
  new("Quantity",
      id = "cveg_dyn",
      name = "Total Carbon mass Flux by Pft",
      units = "kgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon_mass_flux_in_vegetation_by_PFT")
  ),
  
  
  new("Quantity",
      id = "BA_size",
      name = "BA by sizeclass",
      units = "m2 ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("basal Area by sizeclass")
  ),
  
  
  
  
  ########ORCHIDEE -specific variables:
  new("Quantity",
      id = "AGcwood_pft",
      name = "Aboveground woody carbon",
      units = "kgC m-2",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Aboveground_woody_biomas_by_pft")
  ),
  
  ############FATES mortality and other quantities
  #"freezing_stemmort","csarvation_stemmort","background_stemmort","hydraulic_stemmort","impact_stemmort","termination_stemmort"
  #"cmort_hydro","cmort_cstarv", "Bgrowth"
  
  new("Quantity",
      id = "Bgrowth",
      name = "Basal area growth",
      units = "m2 ha-1 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Basal_area_growth")
  ),
  
  new("Quantity",
      id = "stemmort_freezing",
      name = "Stems lost due to freezing mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("freezing mortality by pft/size in number of plants per m2 per year")
  ),
  
  new("Quantity",
      id = "stemmort_freezing_rate",
      name = "% stems lost due to freezing mortality",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("freezing mortality by pft/size in number of plants per m2 per year")
  ),
  
  new("Quantity",
      id = "stemmort_cstarv",
      name = "Stems lost due to carbon starvation mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("carbon starvation mortality by pft/size in number of plants per m2 per year")
  ),
  
  new("Quantity",
      id = "stemmort_cstarv_rate",
      name = "% Stems lost due to carbon starvation mortality",
      units =  "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("carbon starvation mortality by pft/size in number of plants per m2 per year")
  ),
  new("Quantity",
      id = "stemmort_background",
      name = "Stems lost due to background mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("background mortality by pft/size in number of plants per m2 per year")
  ),
  
  new("Quantity",
      id = "stemmort_background_rate",
      name = " %Stems lost due to background mortality",
      units =  "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("background mortality by pft/size in number of plants per m2 per year")
  ),
  new("Quantity",
      id = "stemmort_impact",
      name = "Stems lost duee to impact mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Nstems lost duee to impact mortality")
  ),
  
  new("Quantity",
      id = "stemmort_impact_rate",
      name = "#Stems lost duee to impact mortality",
      units =  "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Nstems lost duee to impact mortality")
  ),
  
  new("Quantity",
      id   = "stemmort_termination",
      name = "Stems lost due to termination mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Stems lost due to termination mortality")
  ),
  
  new("Quantity",
      id   = "stemmort_termination_rate",
      name = "% Stems lost due to termination mortality",
      units =  "% yr-1",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Stems lost due to termination mortality, termination mortality that arises when the number of individuals in a cohort becomes so low as to cause numeric instability. Needham et al 2020, GCB")
  ),
  new("Quantity",
      id   = "stemmort_hydraulic",
      name = "Stems lost due to hydraulic failure(?) mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Stems lost due to hycraulic mortality")
  ),
  
  new("Quantity",
      id   = "stemmort_hydro_rate",
      name = "% Stems lost due to hydraulic failure(?) mortality",
      units =  "% yr-1",
      colours = viridis::viridis,
      format = c("DBEN"),
      standard_name = c("Stems lost due to hycraulic mortality")
  ),
  
  new("Quantity",
      id = "cmort_cstarv",
      name = "Wood Carbon Mass flux lost through carbon starvation mortality",
      units = "KgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Wood Carbon Mass flux lost due to\n carbon starvation mortality")
  ),
  
  new("Quantity",
      id = "cmort_cstarv",
      name = "%Wood Carbon Mass flux lost through carbon starvation mortality, from FATES",
      units = "% ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("%Wood Carbon Mass flux lost due to\n carbon starvation mortality,FATES-specific")
  ),
  
  new("Quantity",
      id = "cmort_cstarv_rate",
      name = " % Wood Carbon Mass lost through carbon starvation mortality",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c(" % Wood Carbon Mass lost due to\n carbon starvation mortality")
  ),
  new("Quantity",
      id = "cmort_hydro",
      name = "Wood Carbon Mass flux lost through hydarulic failur mort.",
      units = "KgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Wood Carbon Mass flux lost due to\n hydraulic failure mortality")
  ),
  
  new("Quantity",
      id = "cmort_hydro_rate",
      name = "% Wood Carbon Mass flux lost through hydarulic failur mort. FATES",
      units = "% m2?",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("%Wood Carbon Mass flux lost due to\n hydraulic failure mortality")
  ),
  
  #
  new("Quantity",
      id = "cmort_noturnover_pft", #needed to check carbon budget closure and create rates for the H0 model/LSM
      name = "Carbon Mass flux of all tissues lost through mortality. FATES",
      units = "kgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass flux lost through total tissue turnover")
  ),
  
  new("Quantity",
      id = "Bgrowth",
      name = "Biomass growth",
      units = "kgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Biomass_growth_FATES")
  ),  
  
  new("Quantity",
      id = "cveg_total",
      name = "cveg_total",# unclear how it differs from cveg. necessary for mortality rates
      units = "kgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("cveg_total_FATES")
  ), 
  
  ###############  BiomeEP mortality quantities
  new("Quantity",
      id = "cmort_pft",
      name = "Carbon Mass Flux lost from live wood due to mortality",
      units = "kg C m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to mortality")
  ),
  ###############  BiomeE-Standalone
  new("Quantity",
      id = "Hstar",
      name = "Hstar",
      units = "kg C m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("critical height, at which canopy is closed")
  ),
  
  
  ############### EDv3 mortality quantities
  #[TODO] Lei EDv3 mortality rate, check units
  new("Quantity",
      id = "cmort_disb",
      name = "carbon loss from trees due to disturbance events like tree fall and fire",
      units = "kg C m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("%Wood Carbon Mass flux lost due to\n hydraulic failure mortality")
  ),
  new("Quantity",
      id = "cmort_mort",
      name = "carbon loss from trees due to tree mortality which is related to tree density and carbon starvation",
      units = "kg C m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("tree density and carbon starvation mortality")
  ),
  
  new("Quantity",
      id = "stemmort_mort",
      name = "stem number decrease due to tree mortality which is related to tree density and carbon starvation",
      units = "count ha-1 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("tree density and carbon starvation mortality")
  ),
  
  new("Quantity",
      id = "stemmort_disb",
      name = "stem number decrease due to disturbance events like tree fall and fire",
      units = "count ha-1 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("%Wood Carbon Mass flux lost due to\n hydraulic failure mortality")
  ),
  
  
  ############### CABLE-POP mortality quantities
  new("Quantity",
      id = "cmort_res",
      name = "Wood Carbon Mass flux lost due to Resource mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Resource mortality")
  ),
  
  new("Quantity",
      id = "cmort_res_rate",
      name ="% Wood Carbon Mass lost, Resource mortality",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% wood carbon mass lost, Resource mortality")
  ),
  
  new("Quantity",
      id = "cmort_crowd",
      name = "Wood Carbon Mass flux lost due to crowding mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Crowding mortality")
  ),  
  
  new("Quantity",
      id = "cmort_crowd_rate",
      name = " % Wood Carbon Mass  lost due to crowding mortality",
      units = "% yr -1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% wood carbon mass lost Crowding mortality")
  ), 
  
  new("Quantity",
      id = "cmort_dist",
      name = "Wood Carbon Mass flux lost due to disturbance mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Disturbance mortality")
  ),
  
  new("Quantity",
      id = "cmort_dist_rate",
      name = "% Wood Carbon Mass lost, disturbance mortality",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Disturbance mortality")
  ),
  
  ####### LPJ-GUESS  mortality and other quantities
  #var_temp <- ncvar_def("cmort_age", "kg C m-2 s-1", dim=list(time,sizeclass_dim), 
  #                      longname="Carbon_Mass_Flux_lost_from_live_wood_due_to_age-based_mortality", mv,verbose=TRUE) 
  
  new("Quantity",
      id = "Cdebt",
      name = "Cdebt",
      units = "kgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Cdebt_borrowed_against_Csap_pool")
  ),
  
  new("Quantity",
      id = "cmort_age",
      name = "Wood Carbon Mass flux lost through old age mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to age-based mortality, by size class")
  ),
  
  new("Quantity",
      id = "cmort_age_rate",
      name = "% Wood Carbon Mass lost, old age mortality",
      units = "% s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to age-based mortality, by size class")
  ),
  
  new("Quantity",
      id = "cmort_dist",
      name = "Wood Carbon Mass flux lost through disturbance mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to disturbance, by size class")
  ),
  
  new("Quantity",
      id = "cmort_thin",
      name = "Wood Carbon Mass flux lost through thinning mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to thinning, by size class")
  ),
  
  new("Quantity",# [TODO] check whether this is harvest thinning or actually self-thinning mortality
      id = "cmort_thin",
      name = "Wood Carbon Mass flux lost through thinning mortality",
      units = "% s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c(" % Carbon Mass lost from live wood due to thinning, by size class")
  ),
  new("Quantity",
      id = "cmort_dist_rate",
      name = "% Wood Carbon Mass lost through disturbance mortality",
      units = "% s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% Carbon Mass lost from live wood due to disturbance, by size class")
  ),
  
  new("Quantity",
      id = "cmort_fire",
      name = "Wood Carbon Mass flux lost through fire mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to fire mortality, by size class")
  ),
  
  new("Quantity",
      id = "cmort_fire_rate",
      name = "% Wood Carbon Mass lost, fire mortality",
      units = "% s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("%Carbon Mass lost from live wood due to fire mortality, by size class")
  ),
  new("Quantity",
      id = "cmort_greff",
      name = "Wood Carbon Mass flux lost, low growth efficiency mortality",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to low growth efficiency mortality, by size class")
  ),
  
  new("Quantity",
      id = "cmort_greff_rate",
      name = "% Wood Carbon Mass lost, growth efficiency mortality",
      units = "% s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c(" % Carbon Mass lost from live wood due to low growth efficiency mortality, by size class")
  ),
  
  new("Quantity",
      id = "cmort_other",
      name = "Wood Carbon Mass flux lost, other mortality mechs.",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to other mortality mechs., by size class")
  ),
  
  new("Quantity",
      id = "cmort_res",
      name = "Wood Carbon Mass flux lost, other mortality mechs.",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Carbon Mass Flux lost from live wood due to other mortality mechs., by size class")
  ),
  
  new("Quantity",
      id = "cmort_other_rate",
      name = " % Wood Carbon Mass lost, other mortality mechs.",
      units = "kg C m-2 s-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% Carbon Mass lost from live wood due to other mortality mechs., by size class")
  ),
  
  
  new("Quantity",
      id = "stemmort_fire",
      name = "Number of stems lost through fire mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Number of stems lost from live wood due to fire mortality, by size class")
  ),
  
  new("Quantity",
      id = "stemmort_fire_rate",
      name = "fire mortality rate",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% stems lost from live wood due to fire mortality, by size class")
  ),
  
  new("Quantity",
      id = "stemmort_greff",
      name = "Number of stems lost through low growth efficiency mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Number of stems lost from live wood due to low growth efficiency mortality, by size class")
  ),
  
  new("Quantity",
      id = "stemmort_greff_rate",
      name = "low growth efficiency mortality rate",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% of stems lost from live wood due to low growth efficiency mortality, by size class")
  ),
  
  
  new("Quantity",
      id = "stemmort_other",
      name = "Number of stems lost through other mortality mechs.",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Number of stems from live wood due to other mortality mechs., by size class")
  ),
  
  
  new("Quantity",
      id = "stemmort_other_rate",
      name = "mortality rate, other mechs.",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% of stems from live wood due to other mortality mechs., by size class")
  ),
  new("Quantity",
      id = "stemmort_age",
      name = "Number of stems lost through age related mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Number of stems from live wood due to age related mortality., by size class")
  ),
  
  new("Quantity",
      id = "stemmort_age_rate",
      name = "old age related mortality rate",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% of stems from live wood due to age related mortality., by size class")
  ),
  
  new("Quantity",
      id = "stemmort_dist",
      name = "Number of stems lost through disturbance mortality",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Number of stems from live wood due to disturbance mortality., by size class")
  ),
  
  new("Quantity",
      id = "stemmort_dist_rate",
      name = " disturbance mortality rate",
      units = "% yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("% of stems from live wood due to disturbance mortality, by size class")
  ),
  
  ##########BiomeE-Standalone quantities for carbon budget closure:
  
  new("Quantity",
      id = "WDgrowTot",
      name= "Total wood growth",
      units = "KgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("woody growth, including reproduction")
  ), # this is a temporary quantity derived from BiomeE output: WDgrow+ WDrepr  (new seedlings from reproduction)+ WDmgrt (migrated plants, once in year 31), just for traceability, that will be updated as WBgrowth, as that is what "Total wood growth" means for DBEN. 
  
  new("Quantity",
      id = "WDmortTot",
      name= "Total wood mortality",
      units = "KgC m-2 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Total wood mortality")
  ), # a temporary quantity derived from BiomeE output: WDmort+WDkill (killed trees because of low density)+WDdstb (removed by disturbance, year 31)
  
  ##########ORCHIDEE quantities for distinguishing between representative cohort and overall output:
  # ORCHIDEE submitted 
  new("Quantity",
      id = "nstem",
      name= "number_of_stems",
      units = "count ha-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("woody growth, including reproduction")
  ), # 
  
  new("Quantity",
      id = "stemmort",
      name= "mortality stems",
      units = "count ha-1 yr-1",
      colours = viridis::viridis,
      format =c("DBEN"),
      standard_name = c("Total wood mortality")
  ) # 
  
)

# Now take a look at them
for(quant in quantities_DBEN ) {
  print(quant)
}

#not sure I need this:
availablequantities_DBEN <- function(x, names = TRUE, additional.args){ 
  
  # typical stuff 
  # * get a list of files in the directory
  # * scan the files for different variables
  # * build a list of qunatities present in the model output / dataset
  
  # dummy code
  Quantities.present <- quantities_DBEN
  return(Quantities.present)
  
} 


#temporary hack- but maybe not even necessary? did not create STA.info object and removed it from the new getField function.
#sta.info is now return.sta.info within the function
#have to deal with this later though

#getField switches between PFTs and sizeclasses automatically
getField_DBEN <- function(source, quant, file.name,verbose = FALSE, model_name = "model_name_placeholder", ...){  
  # for some models, I need post-processing in here, so I need model_name here as switches 
  # varname is added as hack, as otherwise I get an issue with getField function and quant@id towards the end.
  #quant <- availablequantities_DBEN(which)
  # open the file and read the data with ncvar_get() or read.table() or whatever
  
  print(model_name)
  ##Handle .nc output:
  
  if((model_name != "Obs") & (model_name != "seib")){ # don't evaluate for either OBS or SEIB-DGVM
    
    possible.time.dim.names <- c("time", "times","Time", "Times", "year", "Years", "Year", "month", "Month", "t", "T","Base/Time")
    possible.sizeclass.dim.names <- c("Class", "sizeclass","DIAM_BINS","fates_levscls","fates_levscpf")
    # note: "fates_levscpf" includes bofth PFT and sizeclasses
    # [TODO], at some point later, decide how to create a switch for PFT -sizeclass output?
    possible.pft.dim.names <- c("pft","PFT","fates_levpft")
    
    this.nc <- nc_open(file.name)
    nc_var  <-  ncvar_get(this.nc ,quant@id)
    #if(model_name == "BiomeE"){
    #  nc_var  <-  t(ncvar_get(this.nc ,quant@id))
    #}
    
    # pick up sizeclass
    #loop through all dimensions:
    sizeclass.string <- NULL
    for (this.dimension in this.nc$dim){
      
      # adress the fact that CABLE-POP sizeclass dimension cannot be picked up explicitly
      # when reading in netcdfs:
      if(this.dimension$name %in% possible.sizeclass.dim.names) {
        sizeclass.string <- this.dimension$name
      }
    }
    
    # pick up PFT
    #loop through all dimensions:
    pft.string <- NULL
    for (this.dimension in this.nc$dim){
      # adress the fact that CABLE-POP sizeclass dimension cannot be picked up explicitly
      # when reading in netcdfs:
      if(this.dimension$name %in% possible.pft.dim.names) {
        pft.string <- this.dimension$name
      }
    }
    
    # pick up Time
    #loop through all dimensions:
    for (this.dimension in this.nc$dim){
      
      if(this.dimension$name %in% possible.time.dim.names) {
        #print("here, in Time dimension retrieval")
        time.string <- this.dimension$name
        if(model_name == "BiomeE"){
          nc_time  <- ncdf4::ncvar_get(this.nc, "Years")
        }else if(model_name == "FATES" ){ # get year
          nc_time  <- ncdf4::ncvar_get(this.nc, time.string)/365+1
        }else if(model_name =="EDv3"){ # raw data starts with year 0, so doing +1
          nc_time <- ncdf4::ncvar_get(this.nc, time.string)+1 
        }else if(model_name == "CABLE-POP" & quant@id != "albedo"
                 &  quant@id != "cmort_crowd" &   quant@id!= "cmort_dist" 
                 & quant@id !="cmort_res"   & quant@id != "AGcwood"  & quant@id != "nbp"){
          #Address the fact that some variables in CABLE-POP have no time variable:
          nc_time  <- seq(1,dim(nc_var)[2])
        }else{ # LPJ-GUESS, CABLE-POP(some vars), BiomeEP
          nc_time  <- ncdf4::ncvar_get(this.nc, time.string)
        }
        if(verbose) message(paste0("** Confirmed time dimension: '", time.string, "', with ", length(all.time.intervals), " values."))
      }
      # if(model_name=="JULES-RED"){
      #   nc_time <-floor(nc_time)
      # }
      
    }
    
    #AHES legacy from initial/exploratory simulations for Amsterdam; maybe still useful when other models are added, so
    # just keeping it commented out for now.
    #Now modify time variable content to correspond to years in p0 or p1 simulations:
    #if(length(nc_time)<=114){ #p1 simulations
    #  nc_time <- seq(1901,length(nc_time)+1900)
    #}
    #if(length(nc_time) == 116){ #p1 simulations
    #  nc_time <- seq(1899,2014)
    #}else{#p0 simulations.
    #  #here, remain flexible to account for the fact that CABLEPOP does a 700 year simulation
    #  nc_time <- seq(1:length(nc_time))
    #}
    
    
    if(!is.null(pft.string)){ # if pft dimension exists, proceed with the below
      
      #hack_for_JULES-RED -because we only need stemmort as sizeclass for now:
      if(model_name=="JULES-RED" & quant@id =="stemmort"){
        #do nothing here and then continue to sizeclasses
      }else {
        
        df <- data.frame(Year=nc_time)
        # reintroduce manually here, just so that there is no trouble in any DGVMTools-related functions later on:
        if(source@id == "FIN" || source@id == "FI" ){
          df$Lat    <- 62.25
          df$Lon    <- 23.25
        }
        if(source@id == "BIA"){
          df$Lat    <- 52.75
          df$Lon    <- 23.75
        }
        if(source@id == "BCI"){
          df$Lat    <- 9.25
          df$Lon    <- -79.75
        }
        
        
        # another legacy code from initial "Amsterdam" simulations 
        #  if(!is.null(model_name)){
        #    if(model_name == "FATES"){ # move PFTs directly into correct column, to skip a post-processing script
        #      #prepare data frame into which to put all pfts:
        #      df <- cbind(df,as.data.frame(matrix(0,ncol = 6, nrow = length(nc_time) )))
        #      if(source@id == "BCI"){
        #        #bug in raingreen, so is excluded from simulations at this time (email from J Needham, 24.Aug.2022)
        #        idx <- 7:8 
        #      }
        #      if(source@id == "BIA"|source@id =="Bia"){  
        #        idx <- 4:5
        #      }
        #      if(source@id == "FI"|source@id == "Fi2"){
        #        idx <- 4:6
        #      }
        #    }
        #  } 
        # the above legacy code updated for paper simulations:
        if(!is.null(model_name)){
          if(model_name == "FATES"){ # move PFTs directly into correct column, to skip a post-processing script
            #prepare data frame into which to put all pfts:
            df <- cbind(df,as.data.frame(matrix(0,ncol = 8, nrow = length(nc_time) )))
            if(source@id == "BCI"){
              # c("BINE","NE","IBS","TeBS","TrIBE","TrBE","TrBR","Grasses")
              idx <- c(5:8)+3 # +3 to account for the Year, Lat and Lon column, that is already part of the df at that point
            }
            if(source@id == "BIA"|source@id =="Bia"){  
              idx <- c(2:4,8)+3
            }
            if(source@id == "FI"|source@id == "Fi2"|source@id == "FIN"){
              idx <- c(1:3,8)+3
            }
          }
          # if(model_name == "BiomeEP" & site != "FIN"){ # move PFTs directly into correct column, to skip a post-processing script
          #prepare data frame into which to put all pfts:
          #  df <- cbind(df,as.data.frame(matrix(0,ncol = 8, nrow = length(nc_time) )))
          #   if(source@id == "BCI"){
          #     # c("BINE","NE","IBS","TeBS","TrIBE","TrBE","TrBR","Grasses")
          #     idx <- c(5:8)+3 # +3 to account for the Year, Lat and Lon column, that is already part of the df at that point
          #   }
          #   if(source@id == "BIA"|source@id =="Bia"){  
          #     idx <- c(2:4,8)+3
          #   }
          # Finland is in the correct format
          # }
          if(model_name =="EDv3"){
            # mapping:
            #      DBEN            EDv3
            #      PFT1 (B)INE -  1; shade intolerant needleleaf;
            #      PFT2  NE       2: shade tolerant needleleaf;
            #      PFT3  IBS      3: shade intolerant broadleaf
            #      PFT4  TeBS     4: intermediate shade tolerant broadleaf;
            #      PFT5  TrIBE    3: shade intolerant broadleaf
            #      PFT6  TrBE     5: shade tolerant broadleaf;
            #      PFT7  TrBR      ?   
            #      PFT8  Grasses  6:  grasses
            
            #prepare data frame into which to put all pfts:
            df <- cbind(df,as.data.frame(matrix(0,ncol = 8, nrow = length(nc_time) )))
            
            if(source@id == "FIN"){
              idx <- c(1,2,3,8)+3 # +3 to account for the Year, Lat and Lon column, that is already part of the df at that point
            }
            if(source@id == "BIA"|source@id =="Bia"){  
              idx <- c(2:4,8)+3
            }
            if(source@id == "BCI"){
              idx <- c(5,6,8)+3 # no explicit raingreen in EDv3, all broadleaves are raingreen if conditions are bad.
              idx <- c(4:6,8)+3 #testing
            }
            
          }
        }
        
        
        
        # pick up pft names and fill df with values from netcdf files
        # loop through all dimensions:
        for (this.dimension in this.nc$dim){
          
          if(this.dimension$name %in% possible.pft.dim.names) {
            pft.string <- this.dimension$name
            
            if(model_name =="EDv3"){
              #now, fill with nc-data,but depends on the site:
              if(source@id == "FIN"){
                df[,idx] <- nc_var[,c(1:3,6)]
              }
              if(source@id == "BIA"|source@id =="Bia"){ 
                df[,idx] <- nc_var[,c(2:4,6)]
              }
              if(source@id == "BCI"){
                #df[,idx] <- nc_var[,c(3,5:6)]
                df[,idx] <- nc_var[,c(3:6)]
              }
              
            }else if(model_name == "BiomeE" & !is.null(sizeclass.string) ){
              
              
              nc_pft <- ncdf4::ncvar_get(this.nc,quant@id)
              collect_pft <- data.frame(matrix(ncol=8,nrow = this.nc$dim$time$len))
              
              for(i in 1:this.nc$dim$time$len){
                collect_pft[i,1:8] <- rowSums(nc_pft[, , i])
              }
              
              df <- cbind(df,collect_pft) # nstem by pft is now technically possible here.
              
            }else{
              #adress the fact that FATES pft dimension cannot be picked up through standard procedure
              # when reading in netcdfs:
              ## not necessary? 28.04.2023 account for the fact that FATES has stemmort and others also outputted by pft.
              # if(pft.string ==  "fates_levscls" | pft.string ==  "fates_levpft"){ # FATES
              #in Amsteram  data "fates_levscls"  was  "fates_levscpf", and the format was different.
              #if(pft.string ==  "fates_levpft" |pft.string ==  "fates_levscls" ){ # FATES
              if(pft.string ==  "fates_levpft" ){ # FATES
                if(quant@id == "stemmort" | quant@id == "stemmort_cstarv" | quant@id == "stemmort_other"
                   | quant@id == "stemmort_hydraulic"  | quant@id == "stemmort_freezing"  
                   | quant@id == "stemmort_impact"  | quant@id == "stemmort_background"
                   | quant@id == "stemmort_termination" ){
                  nc_pft <- ncdf4::ncvar_get(this.nc,quant@id)
                  collect_pft <- data.frame(matrix(0.0,ncol=8,nrow = this.nc$dim$time$len))
                  
                  for(i in 1:this.nc$dim$time$len){
                    collect_pft[i,(idx-3)] <- colSums(nc_pft[, , i])
                  }
                  # a bit hacky - as above I have already created an empty pft-df, 
                  # but here for stemmort, in the loop, it makes sense
                  df <- cbind(df[1:3],collect_pft) # nstem by pft is no technically possible here.
                  
                }else{# for other vars, which don't have nstem and pft in one file, but only pft:
                  #now, fill with nc-data
                  df[,idx] <- t(nc_var)
                }
                
              }else if(pft.string =="PFT" & model_name =="CABLE-POP"){
                df <- cbind(df,as.data.frame(t(nc_var)))
              }else
              { # LPJGUESS,JULES-RED, BiomeE for PFT-only output
                nc_pft <- ncdf4::ncvar_get(this.nc,pft.string)
                if(model_name == "BiomeE"){
                  df <- cbind(df,as.data.frame(t(nc_var))) #[TODO] should this be nc_pft?!?!?!?! not nc_var?
                }
                else{# LPJGUESS,JULES-RED,BiomeEP
                  #pft + sizeclass dimensions here create problem for JULES-RED, jumped out of the loop, but here JULES-RED probably. needs
                  # a different treatment
                  df <- cbind(df,as.data.frame(nc_var)) #[TODO] should this be nc_pft?!?!?!?! not nc_var?
                  ##NOTE: stemmort pft dimension has a NAN for what I think are grasses
                }
                
                
              } 
              if(verbose) message(paste0("** Confirmed pft dimension: '", pft.string, "', with ", length(all.pft.intervals), " values."))
              
            }
          }
          
        }
        
        
        #add pft names to df
        names(df)[4:11] <- c("BINE","NE","IBS","TeBS","TrIBE","TrBE","TrBR","Grasses")
        
        #JULES-RED has Grasses in the WBgrowth variable, set to 0 before calculating total:
        if(model_name == "JULES-RED" & quant@id =="WBgrowth"){
          df$Grasses = 0
        }
        if(model_name == "JULES-RED" & quant@id =="cwood"){
          df$Grasses = 0
        }
        # EDv3 has values in WBgrowth for grasses, too, set to 0:
        if(model_name == "EDv3" & quant@id =="WBgrowth"){
          df$Grasses = 0
        }
        
        # add a new column:
        #df$Total <- rowSums(df[4:11],na.rm = TRUE)
        #df$Total <- rowSums(df[4:11],na.rm = TRUE)
        df$Total <- rowSums(round(df[4:11],digits = 10),na.rm = TRUE)
        
        
        
        if(!is.null(model_name)){  # Daniel's observations
          if(model_name == "Obs"){
            df <- read.csv(file.name, header=TRUE)
          }
        }
        
      }#hack_for_JULES-RED end
    }
    
    
    
    if(!is.null(sizeclass.string)){ # if sizeclass dimension exists, proceed with the below
      # had to be separate if statement to evaluate BiomeE-Standalone's handling of PFTs and sizeclasses in the same file.
      # [TODO] devise more elegant method, in which we also get BiomeE-PFts outputted (not needed at the moment)
      # pick up sizeclass
      # loop through all dimensions:
      for (this.dimension in this.nc$dim){
        
        if(this.dimension$name %in% possible.sizeclass.dim.names) {
          sizeclass.string <- this.dimension$name
          #print(sizeclass.string)
          # address the fact that CABLE-POP sizeclass dimension cannot be picked up explicitly
          # when reading in netcdfs:
          if(sizeclass.string ==  "DIAM_BINS"){ # cable-pop
            nc_size <- c(1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250)
            
          }else if(sizeclass.string ==  "fates_levscls" |sizeclass.string ==  "fates_levpft"  ){ # FATES
            ### https://github.com/JessicaNeedham/FATES-DBEN/blob/main/Variable_extraction_scripts/BCI_p0.ipynb
            ## removed first size bin (0-1) see box 17. 
            ## " Bin one is 1-5 and the last bin should be >200. " J needham email from 16.04.2023
            nc_size  <- c(1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250) # so only 15 sizeclasses are actually read in, but I add the 16th further below. 
            # [TODO] tidy this... waiting to hear back grom JN re removal of first sizeclass
            # 07.06.2023  suddenly only 15 sizeclasses again. so removed 250 as last sc.
            # 24.09.2024 16 sc again
            
            
          }else if(sizeclass.string == "Class"){ #BiomeE-Standalone 
            
            nc_size  <- c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2)*100 # cm 
          }else
          { # ORCHIDEE, JULES-RED, LPJ-GUESS,..
            nc_size  <- ncdf4::ncvar_get(this.nc, sizeclass.string)
          } 
          if(verbose) message(paste0("** Confirmed sizeclass dimension: '", sizeclass.string, "', with ", length(nc_size), " values."))
        }
        
      }
      
      df <- data.frame(Year=nc_time)
      # reintroduce manually here, just so that there is no trouble in any DGVMTools-related functions later on:
      if(source@id == "FIN" | source@id == "FI" ){
        df$Lat    <- 62.25
        df$Lon    <- 23.25
      }
      if(source@id == "BIA"){
        df$Lat    <- 52.75
        df$Lon    <- 23.75
      }
      if(source@id == "BCI"){
        df$Lat    <- 9.25
        df$Lon    <- -79.75
      }
      
      # adress some models' sizeclass organisation:
      # 28.04.2023 FATES seems to have changed its sizeclass organisation for stemmort and no longer falls into
      # this test. Therefore adding this FATES exception quickly to get on...:
      if(dim(nc_var)[2]>16 | (model_name == "FATES" & quant@id =="stemmort" |  quant@id == "stemmort_cstarv" | quant@id == "stemmort_other"
                              | quant@id == "stemmort_hydraulic"  | quant@id == "stemmort_freezing"  
                              | quant@id == "stemmort_impact"  | quant@id == "stemmort_background"
                              | quant@id == "stemmort_termination")){
        #transpose sizeclass-input:
        # address CABLE-POP-sizeclass organisation
        if(model_name == "CABLE-POP" ){
          df <- cbind(df,as.data.frame(t(nc_var)))
        }
        if( model_name == "FATES"){
          # there are two 'sizeclass' representations: fates_levscpf and fates_levscls, which must be handled differently here:
          # fates_levscls can be directly used, and just needs to be transformed, like in CABLE-POP, but, where upper sizeclasses are missing, we need to add it
          if(sizeclass.string == "fates_levscls"){
            #28.04.2023 FATES seems to have changed its sizeclass organisation for stemmort, adressing this here:
            # 07.06.2023  suddenly only 15 sizeclasses again., so changed to ncol 15
            if(quant@id == "stemmort"|  quant@id == "stemmort_cstarv" | quant@id == "stemmort_other"
               | quant@id == "stemmort_hydraulic"  | quant@id == "stemmort_freezing"  
               | quant@id == "stemmort_impact"  | quant@id == "stemmort_background"
               | quant@id == "stemmort_termination"){
              nc_size <- ncdf4::ncvar_get(this.nc,quant@id)
              collect_sizeclass <- data.frame(matrix(ncol=16,nrow = this.nc$dim$time$len)) # adjust when 16 sc
              
              # retrieve size-class dimension, by summing across all PFTs per year:
              for(i in 1:this.nc$dim$time$len){
                collect_sizeclass[i,1:16] <- rowSums(nc_size[,,i ]) # adjust when 16 sc
              }
              df <- cbind(df,collect_sizeclass)
              # in JN's postprocessing she removed sizeclass <1, but not for stemmort.
              # in that sizeclass there are waay too many trees to be realistic. #
              # I imagine they probably die instantly and may not even be counted as biomass
              # [TODO] make sure this is sensible 
              # set <1 to 0 values. What else could I do?
              #df <- cbind(df[1:3], new_column = 0, df[4:ncol(df)]) # add column in 4th place [TODO] check with Jessie
            }else{
              # df <- cbind(df,data.frame(V16=rep(0,dim(df)[1]) )) # add 16th size class: <1cm, fill with 0s.[TODO] ask JN why removed? solved 24092024
              df <- cbind(df,as.data.frame(t(nc_var)))
            }
            
          }
          
        }# 28.04.2023 JULES-RED has new pft dimension, so here access stemmort-sizeclass dimension
        if(model_name == "JULES-RED" & quant@id =="stemmort"){ 
          
          nc_size <- ncdf4::ncvar_get(this.nc,quant@id)
          collect_sizeclass <- data.frame(matrix(ncol=16,nrow = this.nc$dim$time$len))
          
          # retrieve size-class dimension, by summing across all PFTs per year:
          for(i in 1:this.nc$dim$time$len){
            collect_sizeclass[i,1:16] <- colSums(nc_size[1:7,i,]) # 1:7 to exclude grasses
          }
          df <- cbind(df,collect_sizeclass)
          
          
        }
      }else{
        #account for the fact that BiomeE has nstem_size and others also outputted by pft.
        if(model_name == "BiomeE" & !is.null(pft.string)){
          
          nc_size <- ncdf4::ncvar_get(this.nc,quant@id)
          
          collect_sizeclass <- data.frame(matrix(ncol=16,nrow = this.nc$dim$time$len))
          
          # retrieve size-class dimension:
          for(i in 1:this.nc$dim$time$len){
            collect_sizeclass[i,1:16] <- colSums(nc_size[, , i])
          }
          # for(i in 1:this.nc$dim$time$len){
          #    collect_pft[i,1:8] <- rowSums(nc_pft[, , i])
          # }
          
          
          df <- cbind(df,collect_sizeclass)
          
        }else{# all other models
          df <- cbind(df,as.data.frame(nc_var))
        }
      }
      
      ##[TODO] 28.04.2023 the below this can be done more elegantly by using df for these evaluations I think..
      ## (i.e. as done for first evaualtion with dim(df)) all models should have the same number and types of sizeclasses by this iteration, so some of these exceptions are 
      # probably also no longer necessary.
      
      # Distinguish between ORCHIDEE sizeclasses (3)
      # and the ones from the protocol:
      # and adress the fact that nc_size for BiomeE is multidimensional, so the first test has to be:
      if(model_name != "BiomeE"){
        if(length(nc_size) == 16 | (model_name == "FATES" & quant@id == "stemmort" |  quant@id == "nstem_size" |quant@id == "cwood_size" | quant@id == "stemmort" | quant@id == "stemmort_cstarv" | quant@id == "stemmort_other"
                                    | quant@id == "stemmort_hydraulic"  | quant@id == "stemmort_freezing"  
                                    | quant@id == "stemmort_impact"  | quant@id == "stemmort_background"
                                    | quant@id == "stemmort_termination")){ # CABLE_POP,JULES-RED, LPJ-GUESS, BiomeEP, FATES(paper1, also adjusted for stemmort)
          if(model_name == "FATES" & quant@id == "stemmort" | quant@id == "stemmort_cstarv" | quant@id == "stemmort_other"
             | quant@id == "stemmort_hydraulic"  | quant@id == "stemmort_freezing"  
             | quant@id == "stemmort_impact"  | quant@id == "stemmort_background"
             | quant@id == "stemmort_termination"){
            #adjust to 15 dimensions for FATES again 6.06.2023
            names(df)[4:dim(df)[2]] <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200")
          }else{
            
            names(df)[4:dim(df)[2]] <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200")
          }
          
        }else if(length(nc_size) > 16 ){# JULES-RED: 
          # names(df)[4:dim(df)[2]] <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200")
        }
        else{ #ORCHIDEE - should not be called. deprecated.
          names(df)[4:(4+length(nc_size)-1)] <- c("1", "2", "3")
        }
      }else{# if this is BiomeE:
        names(df)[4:(4+length(nc_size[1,,8])-1)] <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200")
      }
      
      
      #address the fact that there are no <1 sizeclasses in JULES-RED:
      #in stemmmort variable, there is 0 available.
      if(model_name =="JULES-RED"| model_name =="FATES"){ # omit small trees <1 in FATES; they are so many that they interfere with the "real" ecological signals.
        df[,4] <- NA # will be removed later, but here set to NA for rowSums to work downstream ( otherwise NaNf value)
        
      }
      
      # add a new column:
      df$Total <- rowSums(round(df[,4:(dim(df)[2])],digits = 10),na.rm = TRUE)
      
      
      #address the fact that there are no <1 sizeclasses in JULES-RED:
      if(model_name =="JULES-RED" ){
        df[,4] <- NULL # NA for no data. not "no trees"[TODO] allow for NA data in trimFieldsForPlotting, right now I am re-introducing this column later to avoid this error message.
      }# AHES this now contains "something", but Arthur sais it's grasses. trees only occur >3 dbh.
    }
    
    # neither sizeclass nor pft dimension exist, so proceed with the below ( valid for nbp- and for AGcwood /AGB variable): 
    if(is.null(sizeclass.string) & is.null(pft.string)){ 
      #can't do an else statement, because BiomeE can handle both sizeclass and pft strings in one run.
      
      df <- data.frame(Year=nc_time)
      # reintroduce manually here, just so that there is no trouble in any DGVMTools-related functions later on:
      if(source@id == "FI" | source@id == "FIN"){
        df$Lat    <- 62.25
        df$Lon    <- 23.25
      }
      if(source@id == "BIA"){
        df$Lat    <- 52.75
        df$Lon    <- 23.75
      }
      if(source@id == "BCI"){
        df$Lat    <- 9.25
        df$Lon    <- -79.75
      }
      
      
      #transpose sizeclass-input:
      # address CABLE-POP-sizeclass organisation:
      # if(model_name == "JULES-RED"){
      #  df <- cbind(df,as.data.frame(t(nc_var[,,2])))
      # }else{
      #if(model_name == "BiomeE" ){ # just anticipating this here. it may not be correct.
      #    df <- cbind(df,as.data.frame(t(nc_var)))
      #  }else{
      df <- cbind(df,as.data.frame(nc_var))
      # }
      
      names(df)[4]<- "Total"
      
      if(!is.null(model_name)){ # Daniel's observations
        if(model_name == "Obs"){
          df <- read.csv(file.name, header=TRUE)
        }
      }
    }
    
    
    #adress the fact that EDv3 has NaN values. Will be set to 0 e.g. in nstem_size there are NAN values. raised issue on this on github
    if(model_name == "EDv3"){
      df <- replace(df, is.na(df), 0)
    }
    
    
    # code needs to get the data as a data.table so reformat it
    # data.tables should be melted down to one column for Lat, Lon, Year, Month and Day (as appropriate).
    # days run from 1-365, and if data is daily no Month column should be used
    
    # as data table
    dt <- data.table(df)
    
    
    # also an STAInfo object
    # define this properly based on the actual data in the data.table
    return.sta.info <- new("STAInfo",
                           first.year = nc_time[1], # let's say
                           last.year = nc_time[length(nc_time)], # let's say
                           year.aggregate.method = "none", # see NOTE 1 below
                           spatial.extent = extent(dt), # raster::extent function has been define for data.tables, might be useful
                           spatial.extent.id = "site", # see NOTE 2 below
                           spatial.aggregate.method = "none", # see NOTE 1 below
                           subannual.resolution = "daily", # let's say
                           subannual.aggregate.method = "none", # see NOTE 1 below
                           subannual.original = "annually" # let's say
    )
    
    nc_close(this.nc) 
  }
  
  ##post-processing necessary for LPJ-GUESS, because it simulates more than 30 years before disturbance:
  if (model_name=="LPJGUESS"){
    dt <- dt[-c(1:3),]
    
  }
  
  #handle SEIB-DGVM:
  if(model_name=="seib"){
    
    df <- data.frame(Year=1:450)
    # reintroduce manually here, just so that there is no trouble in any DGVMTools-related functions later on:
    if(source@id == "FIN" || source@id == "FI" ){
      df$Lat    <- 62.25
      df$Lon    <- 23.25
    }
    if(source@id == "BIA"){
      df$Lat    <- 52.75
      df$Lon    <- 23.75
    }
    if(source@id == "BCI"){
      df$Lat    <- 9.25
      df$Lon    <- -79.75
    }
    
    df <- cbind(df,read.csv(file.name,header=FALSE))
    #add latlons and sim year back in:
    if(dim(df)[2] < 19){ # PFT variable
      if(dim(df)[2] == 4){ #neither sizeclass nor pft dimension exist, so proceed with the below ( valid for nbp- and for AGcwood variable)
        names(df)[4] <- "Total"
      }else{
        if(quant@id =="cwood" |quant@id =="AGcwood" | quant@id =="CA" | quant@id =="BA" | quant@id =="height" | quant@id =="WBgrowth"
           | quant@id =="BAgrowth"){ # manually add column with 0s for Grasses
          df$Grasses <- 0
        }
        names(df)[4:11] <- c("BINE","NE","IBS","TeBS","TrIBE","TrBE","TrBR","Grasses")
        df$Total <- rowSums(df[,c("BINE","NE","IBS","TeBS","TrIBE","TrBE","TrBR","Grasses")])
      }
      
      
    }
    if(dim(df)[2] == 19){ # sizeclasses variable
      names(df)[4:19] <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200")
      df$Total <- rowSums(df[,c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200")])
      
      # adress warning. these  sizelacces normally have 0 values in them, which get interpreted as integers, but best  have them as numeric, like the rest.
      df$`<1` <- as.numeric(df$`<1`)
      df$`>=200` <- as.numeric( df$`>=200`)
      #str(df)
      #which(str(df)==int)
      df <- data.frame(sapply(df,as.numeric))
      names(df)[4:20] <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=200","Total")# quick hack to have names back
      
    }
    
    
    #create this object so sta.info is filled in correctly:
    #NB- [TODO] the below should probably be consolidated. currently this is created three times in this code. But it works for now, moving on..
    dt <- data.table(df)
    
    return.sta.info <- new("STAInfo",
                           first.year = 1, # let's say
                           last.year = 450, # let's say
                           year.aggregate.method = "none", # see NOTE 1 below
                           spatial.extent = extent(dt), # raster::extent function has been define for data.tables, might be useful
                           spatial.extent.id = "Full", # see NOTE 2 below
                           spatial.aggregate.method = "none", # see NOTE 1 below
                           subannual.resolution = "annual", # let's say
                           subannual.aggregate.method = "none", # see NOTE 1 below
                           subannual.original = "annual" # let's say
    )
  }
  
  # Daniel's Observations
  if(model_name == "Obs"){
    
    df <- read.csv(paste0(source@dir,file.name), header=TRUE)
    if(names(df)[4] == "X.1.5."){# this is to deal with the updated files with uncertainties, but to make it backwards-compatible with the
      # previous files used for the initial benchmarking results workshop in Amsterdam.
      # for now i have removeds the "<1" column,as not in the data
      names(df)[4:dim(df)[2]] <- c( "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", ">=250")  
    }else{
      end <- dim(df)[2]-1 # excl total column
      names(df)[4:end] <- c("<1", "<5", "<10", "<15", "<20", "<30", "<40", "<50", "<60", "<70", "<80", "<90","<100","<150", "<200", "<250")
      
    }
    
    
    # code needs to get the data as a data.table so reformat it
    # data.tables should be melted down to one column for Lat, Lon, Year, Month and Day (as appropriate).
    # days run from 1-365, and if data is daily no Month column should be used
    
    
    dt <- data.table(df)
    
    return.sta.info <- new("STAInfo",
                           first.year = df$Year[1], # let's say
                           last.year = df$Year[length(df$Year)], # let's say
                           year.aggregate.method = "none", # see NOTE 1 below
                           spatial.extent = extent(dt), # raster::extent function has been define for data.tables, might be useful
                           spatial.extent.id = "Full", # see NOTE 2 below
                           spatial.aggregate.method = "none", # see NOTE 1 below
                           subannual.resolution = "daily", # let's say
                           subannual.aggregate.method = "none", # see NOTE 1 below
                           subannual.original = "daily" # let's say
    )
    
  }
  
  # NOTE 1: normally there is no reason to do any aggregation in this function since it will be done later (and doesn't save any disk reading)
  # NOTE 2: if no spatial cropping at this stage set this to be the character string "Full", if spatial cropping was done, set it to the spatail.extent.id argument of the target sta.info object
  
  # make the Field ID based on the STAInfo and 
  field.id <- makeFieldID(source = source, quant.string = quant@id, sta.info = return.sta.info)
  
  
  
  
  # build a new Field
  return.Field <- new("Field",
                      id = field.id,
                      quant = quant,
                      data = dt,
                      source = source,
                      return.sta.info
  )
  
  
  
  
  return(return.Field)
  
} 


DBEN <- new("Format", 
            id = "DBEN",
            predefined.layers = Layers_DBEN, 
            quantities = quantities_DBEN, 
            availableQuantities = availablequantities_DBEN, 
            getField = getField_DBEN)


# function that browses through the list of quantities and returns the metadata for 
# the quantity that corresponds to "var"
get_quantity <- function(var){
  
  i=1
  
  while(quantities_DBEN[[i]]@id != var){
    i=i+1
  }
  
  out <- quantities_DBEN[[i]]
  
  
  return(out)
}

#prints all quantity names possible
print_quantities <- function(info=FALSE){
  ll <- length(quantities_DBEN)
  var = quantities_DBEN[[ll]]@id
  while(quantities_DBEN[[i]]@id != var){
    i=i+1
    
    if(info==FALSE){
      # only print the name of the quantity
      print(quantities_DBEN[[i]]@id)
    }else{ # also additional information
      print(quantities_DBEN[[i]])
    }
    
    
    
  }
  
}




