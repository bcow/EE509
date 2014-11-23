Bisection_ts <- function(ameriflux_df, ameriflux_meta, r_surf){
  
  ############## ------ just for testing  ------ ############################
  
  # ------ Test Case #1  ------ #
  #   load('Level_2_Ameriflux_Formatted/Duke_Forest_Hardwoods.Rda')
  #   ameriflux_df = ameriflux_df[which(ameriflux_df$year == 2007 
  #                                     & ameriflux_df$doy == 100),]
  #   load('Level_2_Ameriflux_Formatted/ameriflux_meta.Rda')
  #   ameriflux_meta = ameriflux_meta[1,]
  #   
  #   N_data = nrow(ameriflux_df)
  #   r_surf = rep(403.4288,N_data)
  # --------------------------- #
  
  # ------ Test Case #2  ------ #
  #   matlab_data = readMat("data_test_dk_dy100_yr2007.mat")
  #   ameriflux_df = matlab_data$data
  # --------------------------- #
  ###########################################################################
  
  source("Calculate_energy_bal_error.R")
  
  # Number of data points 
  N_data = nrow(ameriflux_df)
  
  # energy balance error threshold 
  enbal_error_threshold = 1 # (W/m^2)
  max_number_iterations = 500 
  
  # Final temperature array
  TS = rep(0,length.out = N_data) 
  
  # Initialize surface temp for bisection
  ts1 = rep(100,N_data)
  ts2 = rep(450,N_data)
  
  arg1 = Calculate_energy_bal_error(ameriflux_df, ameriflux_meta, r_surf, ts1)
  arg2 = Calculate_energy_bal_error(ameriflux_df, ameriflux_meta, r_surf, ts2)
  
  bal_err1 = arg1$balance_error
  bal_err2 = arg2$balance_error
  
  Z = 48
  iter = 0
  while(iter < max_number_iterations & Z != 0){
    
    # Find the midpoint of the bounding temperatures
    ts_mid = (ts1 + ts2)/2
    arg_mid = Calculate_energy_bal_error(ameriflux_df, ameriflux_meta, r_surf, ts_mid)
    bal_errmid = arg_mid$balance_error
    
    # Check and see if mid-point is within error.
    II = which(abs(bal_errmid) < enbal_error_threshold)
    TS[II] = ts_mid[II]
    
    II = which(bal_errmid < 0)
    ts2[II]  = ts_mid[II]
    bal_err2[II] = bal_errmid[II] 
    
    II = which(bal_errmid > 0)
    ts1[II] = ts_mid[II] 
    bal_err1[II] = bal_errmid[II] 
    
    Z = length(which(TS==0))
    iter=iter+1
    
  }

if(iter == max_number_iterations){
  warning("Surface temperature bisection did not converge")
}else{
  ts_mid = (ts1 + ts2)/2
  arg_mid = Calculate_energy_bal_error(ameriflux_df, ameriflux_meta, r_surf, ts_mid)
  
  le_final = arg_mid$le
  return(le_final) # return the LE for the given Rsurf =)
}



# ------------- Test Case #1 & #2 should be (#2 exact, #1 approx): ------------- 
# le_mat = c(6.5201,6.8378,6.2153,5.7688,-1.7328,-2.1290,-1.4413,0.2755,-0.2770,-0.3795,-1.2142, 0.7976,
#            -1.4350, 1.6084, 8.3337, 21.2209, 42.5672, 45.1662, 56.6124, 71.5459, 82.3013, 95.3135,
#            109.6752,119.1091,109.0845,145.1613,149.5670,140.2238,132.9496,130.5271,118.4632,
#            91.1841,89.8645,80.5134,53.1395,35.4334, 19.6095,9.2215,7.0197,8.8415,7.3777,11.0465,
#            14.4573,13.6881,17.6692,17.7705,16.8264,15.1921)
# -------------------------------------------------------------------------------
  
} # end function