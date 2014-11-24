Calculate_energy_bal_error <- function(ameriflux_df, ameriflux_meta, r_surf, temp_surf){
  
  N_data = nrow(ameriflux_df)
  
  #------- SH error threshold -------#
  sh_error_threshold = 0.5 # (W/m^2)
  max_number_iterations = 11 
  
  #------- Constants -------#
  BOLTZMAN = 1.380658e-23
  AVOGADRO = .602214199e24
  MD = 28.9644e-3
  MV = 18.0153e-3
  r_v = (AVOGADRO)*(BOLTZMAN)/(MV) # gas constant f or water vapor (J/(kg-K))
  
  r_d = (AVOGADRO)*(BOLTZMAN)/(MD)  # gas constant for dry air (J/(kg-K))
  cp = 7/2*(r_d)    # specific heat of air (J/(kg-K))
  v = 1.4531e-05    # kinematic viscosity (m^2/s]
  k = 0.41          # von Karman constant ()
  g = 9.81          # gravitational acceleration (m/s^2)
  eps = r_d/r_v     # ratio of universal gas constant dry air to water vapor ()
  lv = 2.5008e6     # latent heat of vaporization (J/kg)
  emis = 0.98       # emissivity of the ground-surface ()
  sb = 5.6704*10^(-8)   # stephan-boltzman constant (W/(m^2 K^4))
  
  # Extract met and flux data from ameriflux_df
  temp = ameriflux_df$ta # air temp (K)
  q = ameriflux_df$q # specific humidity
  p = ameriflux_df$press # pressure (Pa)
  u_str = ameriflux_df$u_str # ustr (m/s)
  
  net_SW = ameriflux_df$rg - ameriflux_df$rgout # netSW = SWin - SWout
  net_SW[net_SW<0]=0 # make sure all net_SW > 0
  r_ld = ameriflux_df$rgl # longwave down
  r_ld[r_ld<0]=0  # make sure all longwave down > 0
  r_ld = r_ld - ((1-emis)*r_ld) # subract out relected longwave
  ground = ameriflux_df$gf # ground heat flux
  
  # Extract constand site parameters
  z_veg = ameriflux_meta$z_veg # vegetation height
  z_m = ameriflux_meta$z_m  # measurement height
  KB = ameriflux_meta$KB  # KB^-1 
  
  ## ------------- Test Case #2 -------------  
  #   temp = ameriflux_df[,1]
  #   q =  ameriflux_df[,2]
  #   net_SW =  ameriflux_df[,3]
  #   r_ld = ameriflux_df[,4]
  #   p = ameriflux_df[,7]
  #   ground = ameriflux_df[,9]
  #   u_str = ameriflux_df[,10]
  # -----------------------------------------
  
  # Caluclate density
  den_m = p/(r_d*temp)
  
  # Calculate roughness heights
  z_o = 0.1*z_veg; # momentum roughness height(m)
  d = 0.7*z_veg; # displacement height (m)
  z_ov = (z_o/exp(KB))  # vapor roughness height(m)
  z_oh = (z_o/exp(KB))  # heat roughness height(m)
  
  # Near the surface, the potential temperature is ~= to the temperature.
  ts_pot = temp_surf
  
  # Guess L = infinity (so that the atmoshpere is neutral and stability
  # functions go to zero). Note: L is not a function of z.
  L = 10^10
  i = 0
  
  h_diff = 100
  h = matrix(10^10,N_data,max_number_iterations+1)
  et = matrix(10^10,N_data,max_number_iterations+1)
  
  while((h_diff > sh_error_threshold) & (i < max_number_iterations)){
    i = i + 1
    
    # Calculate three dimensionless hieghts
    xi_vh = ((z_m-d)/L); # dummy-variable for vapor and temp
    xi_v = ((z_ov)/L); # dummy-variable for vapor
    xi_h = ((z_oh)/L); # dummy-variable for temp
    
    # Determine stability based on dimensionless hieghts (above) and equation(s) S4
    # Stabilitly...
    stbl_vh = rep(NA,length.out = N_data)
    stbl_v = rep(NA,length.out = N_data)
    stbl_h = rep(NA,length.out = N_data)
    
    phi_vh = rep(NA,length.out = N_data)
    phi_v = rep(NA,length.out = N_data)
    phi_h =  rep(NA,length.out = N_data)
    
    # VH
    II=which(xi_vh<=0)
    phi_vh[II]=(1-(16*xi_vh[II]))^0.5
    stbl_vh[II] = 2*log((1+phi_vh[II])/2)
    stbl_vh[xi_vh>0&xi_vh<=1] = -5*xi_vh[xi_vh>0&xi_vh<=1]
    stbl_vh[xi_vh>1] =  -5-(5*log(xi_vh[xi_vh>1]))
    
    # V
    II=which(xi_v<=0)
    phi_v[II]=(1-(16*xi_v[II]))^0.5
    stbl_v[II] = 2*log((1+phi_v[II])/2)
    stbl_v[xi_v>0&xi_v<=1] = -5*xi_v[xi_v>0&xi_v<=1]
    stbl_v[xi_v>1] =  -5-(5*log(xi_v[xi_v>1]))
    
    # H
    II=which(xi_h<=0)
    phi_h[II]=(1-(16*xi_h[II]))^0.5
    stbl_h[II] = 2*log((1+phi_h[II])/2)
    stbl_h[xi_h>0&xi_h<=1] = -5*xi_h[xi_h>0&xi_h<=1]
    stbl_h[xi_h>1] =  -5-(5*log(xi_h[xi_h>1]))
    
    # Second part of equation
    chiv = (log((z_m-d)/z_ov)- stbl_vh + stbl_v)
    chih = (log((z_m-d)/z_oh)- stbl_vh + stbl_h)
    
    # Calculate sensible heat
    h[,(i+1)] = ((den_m*cp*k*u_str)*(ts_pot-temp))/chih
    
    # Calculate saturateion vapor pressure using the integrated
    # clausius clayperon relation
    es_sat = 611.2*exp((17.67*(ts_pot-273.15))/(ts_pot-29.65)) #S10
    
    #TO KEEP FROM BOILING %%%%%
    es_sat = min(es_sat,0.95*(p))
    
    # Calculate saturated specific humidity at the surface
    qs_sat = (eps*es_sat)/(p-((1-eps)*es_sat)) #S9
    
    # Calculate et
    et[,(i+1)] = (qs_sat-q)/((r_surf/den_m)+(chiv/(k*u_str*den_m)))
    
    # Use the sensible heat flux to update L (S2) --> iterate until L equillibrates.
    L = -(den_m*cp*(u_str*u_str*u_str)*temp)/(k*g*h[,(i+1)]);
    
    h_diff = max(max(abs(h[,(i+1)] - h[,i])))
    
    if((h_diff <= sh_error_threshold) || (i > (max_number_iterations-1))){
      
      # Calculate three dimensionless hieghts
      xi_vh = ((z_m-d)/L) # dummy-variable for vapor and temp
      xi_v = ((z_ov)/L) # dummy-variable for vapor
      xi_h = ((z_oh)/L) # dummy-variable for temp
      
      stbl_vh = rep(NA,length.out = N_data)
      stbl_v = rep(NA,length.out = N_data)
      stbl_h = rep(NA,length.out = N_data)
      
      phi_vh = rep(NA,length.out = N_data)
      phi_v = rep(NA,length.out = N_data)
      phi_h =  rep(NA,length.out = N_data)
      
      # VH
      II=which(xi_vh<=0)
      phi_vh[II]=(1-(16*xi_vh[II]))^0.5
      stbl_vh[II] = 2*log((1+phi_vh[II])/2)
      stbl_vh[xi_vh>0&xi_vh<=1] = -5*xi_vh[xi_vh>0&xi_vh<=1]
      stbl_vh[xi_vh>1] =  -5-(5*log(xi_vh[xi_vh>1]))
      
      # V
      II=which(xi_v<=0)
      phi_v[II]=(1-(16*xi_v[II]))^0.5
      stbl_v[II] = 2*log((1+phi_v[II])/2)
      stbl_v[xi_v>0&xi_v<=1] = -5*xi_v[xi_v>0&xi_v<=1]
      stbl_v[xi_v>1] =  -5-(5*log(xi_v[xi_v>1]))
      
      # H
      II=which(xi_h<=0)
      phi_h[II]=(1-(16*xi_h[II]))^0.5
      stbl_h[II] = 2*log((1+phi_h[II])/2)
      stbl_h[xi_h>0&xi_h<=1] = -5*xi_h[xi_h>0&xi_h<=1]
      stbl_h[xi_h>1] =  -5-(5*log(xi_h[xi_h>1]))
      
      # Second part of equation
      chiv = (log((z_m-d)/z_ov)- stbl_vh + stbl_v)
      chih = (log((z_m-d)/z_oh)- stbl_vh + stbl_h)
      
      # Calculate sensible heat
      h_final = ((den_m*cp*k*u_str)*(ts_pot-temp))/chih;
      
      # Calculate saturateion vapor pressure using the integrated
      # clausius clayperon relation
      es_sat = 611.2*exp((17.67*(ts_pot-273.15))/(ts_pot-29.65)) #S10
      
      # TO KEEP FROM BOILING #
      es_sat = apply(cbind((es_sat),(0.95*p)),1,min)
      
      # Calculate saturated specific humidity at the surface
      qs_sat = (eps*es_sat)/(p-((1-eps)*es_sat)) #S9
      
      # Calculate et
      et_final = (qs_sat-q)/((r_surf/den_m)+(chiv/(k*u_str*den_m)));
      
    } # end if loop
    
  } # end while loop
  
  # Calculate latent heat
  le = lv*et_final;
  
  # Calculate sensible heat
  r_lu = emis*sb*(ts_pot*ts_pot*ts_pot*ts_pot);
  
  # Calculate energy balance error
  balance_error = net_SW + r_ld - r_lu - le - ground - h_final
  
  balance_error_df = data.frame(balance_error = balance_error, le = le, h = h_final)
  
  return(balance_error_df)
  
} # end funtion



