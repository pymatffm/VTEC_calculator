import numpy as np
import constants as cs
import vtec_plotting


class Calculate_VTEC:
        def __init__(self, RINEXfile, elevation_array):                   
            self.calculation_STEC(RINEXfile, elevation_array)

       def calculation_STEC(self, RINEXfile, elevation_array):
            for sv in range(self.itemSize):
                    sat  = self.obs_data_chunks_dataframe[sv, :]
                    satNumberIndex = self.obs_data_chunks_dataframe.items
                    satNumber = satNumberIndex[sv]
                    
                    gnssType = satNumber[0]
                    satelliteNumber = int(satNumber[1:].lstrip('0'))

                    if gnssType == 'G':
                        self.gnssType = 'GPS'
                        signal_one = 'L1'                                                 
                        signal_two = 'L2'  
                        pseudorange_1 = 'C1'
                        pseudorange_2 = 'P2'

                        CP1 = sat[signal_one]
                        CP2 = sat[signal_two]
                        phi1 = CP1 * cs.GPS_LAMBDA_1         #Change units of L1 to meters
                        phi2 = CP2 * cs.GPS_LAMBDA_2         #Change units of L2 to meters
                        PR1 = sat[pseudorange_1]
                        PR2 = sat[pseudorange_2]
    
                        phase_diff_pseudoranges = PR1 - PR2
                        phase_diff = phi2 - phi1
                        # Set carrier phase to zero at the start due to integer ambiguity. This is relative phase.
                        if phase_diff.isnull().all():
                            print "Sat {0} contains only NaNs for the phase_diff array. Consider using just L1-L2 values.".format(satNumber)
                            scale_to_zero = np.NaN
                        else:
                            scale_to_zero = next(x for x in phase_diff if not isnan(x))           
                            phase_diff = phase_diff - scale_to_zero 

                        iono_teqc = cs.COEFF_GPS * phase_diff 
                        iono_teqc_original = iono_teqc
                        iono_teqc_pr = cs.COEFF_GPS * phase_diff_pseudoranges
                        
                        #TODO: Retrieve the DCM (m) value for this particular PRN                            
                        #dcb_value = self.final_dict.get(satNumber)
                        #iono_teqc_pr = iono_teqc_pr - dcb_value
                        
                        sv_column = satelliteNumber -1
                        sv_ele_array = pd.Series(elevation_array[sv_column])
                        
                        sv_ele_array = sv_ele_array.reset_index(drop=True)

                        high_ele_index = np.argmax(sv_ele_array) 
                                   
                        high_ele_value_pr = iono_teqc_pr.iloc[high_ele_index]            
                        high_ele_value_cp = iono_teqc.iloc[high_ele_index]
                        shift_value = abs(high_ele_value_pr - high_ele_value_cp)
                        iono_teqc = iono_teqc + shift_value
                    
                        TECU_units = iono_teqc / cs.TECU_L1L2

                            
                    elif gnssType == 'R':
                        self.gnssType = 'GLONASS'
                        signal_one = 'L1'                                                 
                        signal_two = 'L2'  
                        pseudorange_1 = 'C1'
                        pseudorange_2 = 'P2'

                        GLONASS_FREQUENCIES, COEFF_GLONASS, TECU_L1L2_glo  = cs.glonass_frequencies(satNumber)
                        GLONASS_L1 = GLONASS_FREQUENCIES[0]
                        GLONASS_L2 = GLONASS_FREQUENCIES[1]
                        RUS_LAMBDA_1 = cs.LIGHT_V / GLONASS_L1
                        RUS_LAMBDA_2 = cs.LIGHT_V / GLONASS_L2
                              
                        CP1 = sat[signal_one]
                        CP2 = sat[signal_two]
                        phi1 = CP1 * RUS_LAMBDA_1          
                        phi2 = CP2 * RUS_LAMBDA_2          
                        PR1 = sat[pseudorange_1]
                        PR2 = sat[pseudorange_2]

                        phase_diff_pseudoranges = PR1 - PR2
                        phase_diff = phi2 - phi1
                        if phase_diff.isnull().all():
                            print "Sat {0} contains only NaNs for the phase_diff array.".format(satNumber)
                            scale_to_zero = np.NaN
                        else:
                            scale_to_zero = next(x for x in phase_diff if not isnan(x))   
                        phase_diff = phase_diff - scale_to_zero 
    
                        iono_teqc = COEFF_GLONASS * phase_diff  
                        
                        # Retrieve the DCM (m) value for this particular PRN                            
                        dcb_value = self.final_dict.get(satNumber)
                        iono_teqc = iono_teqc - dcb_value
                        
                        sv_column = satelliteNumber -1  
                        sv_ele_array = pd.Series(elevation_array[sv_column])    
                        
                        sv_ele_array = sv_ele_array.reset_index(drop=True)
                                 
                        high_ele_index = np.argmax(sv_ele_array)     
                        high_ele_value_pr = phase_diff_pseudoranges.iloc[high_ele_index]
                        high_ele_value_cp = iono_teqc.iloc[high_ele_index]
                        shift_value = abs(high_ele_value_pr - high_ele_value_cp)
                        iono_teqc = iono_teqc + shift_value
                        
                        TECU_units = iono_teqc / TECU_L1L2_glo                        
                        
                    #TODO: Complete this block. Not urgent till there are Galileo svs in view for RINEXv2
                    elif gnssType == 'E':
                        self.gnssType = 'GALILEO'
                        self.galileo_signal_1 = 'L1'
                        self.galileo_signal_2 = 'L8'

                        phi1 = sat[self.galileo_signal_1] * cs.GALILEO_LAMBDA_E1 #Change both parameters if signal is changed
                        phi2 = sat[self.galileo_signal_2] * cs.GALILEO_LAMBDA_E5  #Change both parameters if signal is changed
                        phase_diff = phi1 - phi2
                        iono_teqc = cs.COEFF_GALILEO_GAMMA * phase_diff
                        TECU_units = iono_teqc / cs.TECU_E1E5

                    self.calculation_VTEC(TECU_units, satelliteNumber, satNumber, signal_one, signal_two) 
            
            print "Processing of both STEC and VTEC plots completed." 


        def calculation_VTEC(self, TECU_units, satelliteNumber, satNumber, signal_one, signal_two):

            mf_array = self.mf_dataframe 
            revised_mf_array = []
            revised_mf_array.append(pd.DataFrame(mf_array).dropna(axis=0, how='all').dropna(axis=1, how='all'))
            revised_mf_array = revised_mf_array[0] #This variable is a dataframe
            
            iono_delay_STEC = TECU_units.values
            timeline = TECU_units.index
            
            column = satelliteNumber -1
            specific_sv_mf = revised_mf_array[column]  
            VTEC = iono_delay_STEC * specific_sv_mf
            
            vtec_plotting.Calculate_VTEC(timeline, VTEC, satNumber, signal_one, signal_two)
            
            
  

