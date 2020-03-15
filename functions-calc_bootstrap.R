# functions-calc_bootstrap

# The bootstrapping functions to derive GPP and Rs

# Do the static (no random draws) steps that calculate GPP
# We define a function to re-use later for variance decomposition
calc_bootstrap_gpp <- function(x) {
  x %>%
    mutate(
      # Scenario 2: Rc separate into groups: agriculture, deciduous forest, grassland, etc
      # all other ecosystems with obs<40, take samples from all RC records
      # area of rest ecosystem = 1-sum(agriculture, forest, grassland, etc)
      Rc2 = (Rc_agr * ag_area + Rc_df * df_area + Rc_ef * ef_area + Rc_mf * mf_area +
               Rc_gra * gra_area + Rc_shr * shr_area + Rc_rest *
               (1 - ag_area - df_area - ef_area - mf_area - gra_area - shr_area)),
      
      FsFr = (1 - Froot) / Froot, # calculate Froot to Fshoot ratio of scenario 1
      # Bottom up estimate of GPP, scenario 1
      Rroot = Rs_raw * Rc, # root respiration
      Rshoot = Rroot * FsFr, # shoot respiration
      
      # Scenario 2: separate into different ecosystems and weighted by area
      # Froot %>% count(Ecosystem)
      Froot2 = Froot_df * df_area + Froot_ef * ef_area + Froot_mf * mf_area +
        Froot_rest * (1 - df_area - ef_area - mf_area),
      FsFr2 = (1 - Froot2) / Froot2, # Froot to Fshoot ratio of scenario 2
      # Bottom up estimate of GPP, scenario 2
      Rroot2 = Rs_raw * Rc2,
      Rshoot2 = Rroot2 * FsFr2,
      
      GPP = NPP + Rroot + Rshoot, # scenario 1: GPP
      GPP2 = NPP + Rroot2 + Rshoot2 # scenario 2: GPP2
    )
}

# Do the static (no random draws) steps that calculate Rs
# We define a function to re-use later for variance decomposition
calc_bootstrap_rs <- function(x) {
  x %>%
    mutate(
      # Scenario 1
      RaGpp = RaGpp_rest,
      # Scenario 2
      RaGpp2 = RaGpp_df * df_area + RaGpp_ef * ef_area + RaGpp_mf * mf_area +
        RaGpp_gra * gra_area + RaGpp_rest * (1 - df_area - ef_area - mf_area - gra_area),
      
      # get Froot (weight by area of ecosystem type)
      Froot2 = Froot_df * df_area + Froot_ef * ef_area + Froot_mf * mf_area +
        Froot_rest * (1 - df_area - ef_area - mf_area),
      
      Ra1 = GPP_raw * RaGpp, # first way to calculate Ra, Ra = GPP * RaGPP
      Ra2 = GPP_raw - NPP, # second way to calculate Ra, Ra = GPP - NPP
      Ra3 = if_else(Ra2 < 0, Ra1, Ra2),
      Ra_avg1 = (Ra1 + Ra3) / 2,
      
      Rroot = Ra_avg1 * Froot2,
      Rshoot = Ra_avg1 * (1 - Froot2),
      Rs_topdown = NPP - HerbComsum - Fire - sink - DOC - BVOCs + Rroot,
      
      # Scenario 2
      Ra4 = GPP_raw * RaGpp2, # first way to calculate Ra, Ra = GPP * RaGPP
      Ra_avg2 = (Ra4 + Ra3) / 2,
      Rroot2 = Ra_avg2 * Froot2,
      Rshoot2 = Ra_avg2 * (1 - Froot2),
      Rs_topdown2 = NPP - HerbComsum - Fire - sink - DOC - BVOCs + Rroot2
    )
}
