import numpy as np
import matplotlib.pylab as plt
import matplotlib.colors as mcolors
import BMX
import scipy.signal as sps
import numpy as np
import sys
import os


def high_Pass_Filter(data, filter_Freq):
    #########################################################
    """
        Beginnings of a filter function.
    """
    #########################################################
    #filter_Freq = 1e6/(40)  # cutoff frequency
    fs = 125e6 # This is 1/(time step) 8ns
    N = 5  # Order
    Wn = 2.0*filter_Freq/fs
    B, A = sps.butter(N, Wn, btype = 'highpass')
    output = sps.filtfilt(B, A, data)
    
    return output

def process_Data(data, start_Dex, apply_Filt = False, filt = 100e3):
    ##############################
    ##############################
    """
        A function that does the 
        common theme of filtering
        and mean subtracting.
    """
    ##############################
    ##############################
    data = data[start_Dex:]
    if apply_Filt:
        data = high_Pass_Filter(data, filt)
    #data -= np.mean(data)
    
    return data

def create_Directory(path):
    
    try:
        if os.path.isdir(path):
            pass
        else:
            os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)
        
    return None

def load_Fluctuation_Data(pos, shot, data):
    #####################################################
    """ Loads the fluctuation data """
    #####################################################
    if (data == '11182019'):
        print('Not setup yet')
        sys.exit()
        """data_Path =  '../11182019/'
        filename1 = data_Path + 'port' + str(pos) + '/20191118-(Pos' + str(pos)+ 'S'+ str(shot) + ').npy'"""
    elif (data == '10142019'):
        print('Not setup yet')
        sys.exit()
        """data_Path =  '../10142019/'
        filename1 = data_Path + 'port' + str(pos) + '/20191014-(Pos' + str(pos)+ 'S'+ str(shot) + ').npy'"""
    if (data == '12092019'):
        data_Path =  'fluctuation_Data/'
        filename1 = data_Path + 'port' + str(pos) + '/20191209-(Pos' + str(pos)+ 'S'+ str(shot) + ').npy'
    elif (data == '03102020'):
        data_Path =  'Fluctuation_Data/'
        filename1 = data_Path + 'port' + str(pos) + '/20200310-(Pos' + str(pos)+ 'S'+ str(shot) + ')-fluctuation.npy'
    data1 = np.load(filename1, allow_pickle = True)
    time = data1.item().get('time')
    Br1 = data1.item().get('br')
    #Br1 = process_Data(Br1, time)
    Bt1 = data1.item().get('bt')
    #Bt1 = process_Data(Bt1, time)
    Bz1 = data1.item().get('bz')
    #Bz1 = process_Data(Bz1, time)
    
    return Br1, Bt1, Bz1, time

def correlation(sing_1, sing_2, time, normalized = True):
    ####################################################################
    """ Different Normalization based on the number of point considered """
    ####################################################################
    corr = np.correlate(sing_1, sing_2, mode = 'same')
    dt = time[1]-time[0]
    tau = dt*(np.arange(corr.size) - corr.size/2)

    if normalized:
        normalization = np.zeros(tau.shape)    
        time_Window = time[-1] - time[0]
        if (time_Window - np.abs(tau).any())/dt >= 16000:
            normalization = dt/(time_Window - np.abs(tau) - dt)
        else:
            normalization = 1/16000
        corr = normalization*corr
        
    return tau ,corr

def pos_Norm(pos1, pos2, filt_Data = False, which_Filt = 1e5):
    ######################################################################################
    ######################################################################################
    """
        Calculating the peak normalization. That is, the correlations are normalized
        by the maximums of the two-point correlations based on their global maximum.
        
        Additions: 03/06/2020
        I am adding the standard deviation computation
    """
    ######################################################################################
    ######################################################################################
    ######################################################################################
    """
        The lines below are to find the correct
        path to the data.
    """
    ######################################################################################
    if filt_Data:
        file_Path = 'corr_Data/Filter_' + str(which_Filt) + '/'
    else:
        file_Path = 'corr_Data/Unfilt/'
    ######################################################################################
    if (pos1 == pos2):
        filename = file_Path + 'autocorrelation_P'+ str(pos1) + '_Dict.npy'
    else:
        filename = file_Path + 'correlation_P' + str(pos1) + '-' + str(pos2) + '_Dict.npy'
    my_Dict = np.load(filename, allow_pickle = True)
    
    corr_R = my_Dict.item().get('corr_R')    # < corr_R >
    corr_T = my_Dict.item().get('corr_T')
    corr_Z = my_Dict.item().get('corr_Z')
    ######################################################################################
    """ The standard deviation calculation """
    ######################################################################################
    
    corr_2R = my_Dict.item().get('corr_2R')   # < corr_R**2 >
    corr_2T = my_Dict.item().get('corr_2T')
    corr_2Z = my_Dict.item().get('corr_2Z')
    
    norm = np.max(np.abs(corr_R)) + np.max(np.abs(corr_T)) + np.max(np.abs(corr_Z))
    #print(corr_R.shape)
    #print(corr_T.shape)
    #print(corr_Z.shape)
    #corr_Total = np.abs(corr_R) + np.abs(corr_T) + np.abs(corr_Z)  # corr_Total is direction averaged, in a sense
    #norm = np.max(corr_Total)
    #index_Max = np.where(corr_Total == norm)
    
    #std_R = np.sqrt(corr_2R[index_Max[0]] - corr_R[index_Max[0]]**2)
    #std_T = np.sqrt(corr_2T[index_Max[0]] - corr_T[index_Max[0]]**2)
    #std_Z = np.sqrt(corr_2Z[index_Max[0]] - corr_Z[index_Max[0]]**2)
    
    norm_STD = norm*0.1   #Assuming 10% error
    
    return norm, norm_STD

def pos_Component_Norm(pos1, pos2, component, filt_Data = False, which_Filt = 1e5):
    ######################################################################################
    ######################################################################################
    """
        Calculating the peak normalization. That is, the correlations are normalized
        by the maximums of the two-point correlations based on their global maximum.
        
        Additions: 03/06/2020
        I am adding the standard deviation computation
    """
    ######################################################################################
    ######################################################################################
    ######################################################################################
    """
        The lines below are to find the correct
        path to the data.
    """
    ######################################################################################
    if filt_Data:
        file_Path = 'corr_Data/Filter_' + str(which_Filt) + '/'
    else:
        file_Path = 'corr_Data/Unfilt/'
    ######################################################################################
    if (pos1 == pos2):
        filename = file_Path + 'autocorrelation_P'+ str(pos1) + '_Dict.npy'
    else:
        filename = file_Path + 'correlation_P' + str(pos1) + '-' + str(pos2) + '_Dict.npy'
    my_Dict = np.load(filename, allow_pickle = True)
    
    corr_R = my_Dict.item().get('corr_R')    # < corr_R >
    corr_T = my_Dict.item().get('corr_T')
    corr_Z = my_Dict.item().get('corr_Z')
    ######################################################################################
    """ The standard deviation calculation """
    ######################################################################################
    
    corr_2R = my_Dict.item().get('corr_2R')   # < corr_R**2 >
    corr_2T = my_Dict.item().get('corr_2T')
    corr_2Z = my_Dict.item().get('corr_2Z')
    if (component.lower() == 'r'):
        norm = np.max(np.abs(corr_R))
    elif (component.lower() == 't'):
        norm = np.max(np.abs(corr_T))
    elif (component.lower() == 'z'):
        norm = np.max(np.abs(corr_Z))
    #norm = np.max(np.abs(corr_R)) + np.max(np.abs(corr_T)) + np.max(np.abs(corr_Z))
    #print(corr_R.shape)
    #print(corr_T.shape)
    #print(corr_Z.shape)
    #corr_Total = np.abs(corr_R) + np.abs(corr_T) + np.abs(corr_Z)  # corr_Total is direction averaged, in a sense
    #norm = np.max(corr_Total)
    #index_Max = np.where(corr_Total == norm)
    
    #std_R = np.sqrt(corr_2R[index_Max[0]] - corr_R[index_Max[0]]**2)
    #std_T = np.sqrt(corr_2T[index_Max[0]] - corr_T[index_Max[0]]**2)
    #std_Z = np.sqrt(corr_2Z[index_Max[0]] - corr_Z[index_Max[0]]**2)
    
    norm_STD = norm*0.1   #Assuming 10% error
    
    return norm, norm_STD

def two_Point_Correlations(pos1, pos2, apply_Filt = False, filt = 100e3):
    ########################################################################################
    ########################################################################################
    """
        I want to create a function that generates the ensemble correlation, unnormalized.
    """
    ########################################################################################
    ########################################################################################
    print(pos1)
    print(pos2)
    #data_1 = '11182019'
    data_1 = '03102020'
    ########################################################################################
    """
        The lines below are setting up the shots corresponding to each dataset of 1.5ms 
        stuffing delay.
    """
    ########################################################################################
    
    avoid_Shots = []
    all_Shots = list(range(1, 26))
    considered_Shots = list(set(all_Shots).symmetric_difference(set(avoid_Shots)))


    for shot in considered_Shots:
    
        Br_1, Bt_1, Bz_1, time = load_Fluctuation_Data(pos1, shot, data_1)
        Br_2, Bt_2, Bz_2, time = load_Fluctuation_Data(pos2, shot, data_1)
        
        if (shot == considered_Shots[0]):
            start_Dex = BMX.finding_Index_Time(time*1e-6, 50)
        time = time[start_Dex:]
        
        
        Br_1 = process_Data(Br_1, start_Dex, apply_Filt, filt)
        Bt_1 = process_Data(Bt_1, start_Dex, apply_Filt, filt)
        Bz_1 = process_Data(Bz_1, start_Dex, apply_Filt, filt)
        
        Br_2 = process_Data(Br_2, start_Dex, apply_Filt, filt)
        Bt_2 = process_Data(Bt_2, start_Dex, apply_Filt, filt)
        Bz_2 = process_Data(Bz_2, start_Dex, apply_Filt, filt)
        
        #### Fixing probe components ####
        if (pos1 == 7 or pos1 == 5):
            dum = Bt_1
            Bt_1 = Bz_1
            Bz_1 = dum
        if (pos2 == 7 or pos2 == 5):
            dum = Bt_2
            Bt_2 = Bz_2
            Bz_2 = dum
            
        if (shot == considered_Shots[0]):
            print("Went through")
            tau, corr_R = correlation(Br_1, Br_2, time)
            tau, corr_T = correlation(Bt_1, Bt_2, time)
            tau, corr_Z = correlation(Bz_1, Bz_2, time)
            
            corr_2R = corr_R**2
            corr_2T = corr_T**2
            corr_2Z = corr_Z**2
            
        else:
            tau, corr_R_Hold = correlation(Br_1, Br_2, time)
            tau, corr_T_Hold = correlation(Bt_1, Bt_2, time)
            tau, corr_Z_Hold = correlation(Bz_1, Bz_2, time)
            
            corr_R += corr_R_Hold
            corr_T += corr_T_Hold
            corr_Z += corr_Z_Hold
            
            corr_2R += corr_R**2
            corr_2T += corr_T**2
            corr_2Z += corr_Z**2
    
    corr_R = corr_R/(len(considered_Shots) - 1)
    corr_T = corr_T/(len(considered_Shots) - 1)
    corr_Z = corr_Z/(len(considered_Shots) - 1)
    
    corr_2R = corr_2R/(len(considered_Shots) - 1)
    corr_2T = corr_2T/(len(considered_Shots) - 1)
    corr_2Z = corr_2Z/(len(considered_Shots) - 1)
    
    print(corr_R.shape)
    print(corr_T.shape)
    print(corr_Z.shape)
    
    if (corr_R.shape != corr_T.shape):
        print('R and T')
        print(pos1)
        print(pos2)
    elif (corr_R.shape != corr_Z.shape):
        print('R and Z')
        print(pos1)
        print(pos2)
    elif (corr_T.shape != corr_Z.shape):
        print('T and Z')
        print(pos1)
        print(pos2)
        
    my_Dict = {
        'corr_R': corr_R, 'corr_T': corr_T, 'corr_Z': corr_Z, 'tau': tau, 'corr_2R': corr_2R,
        'corr_2T': corr_2T, 'corr_2Z': corr_2Z
    }
    
    if apply_Filt:
        path = 'Filter_%s'%filt
        create_Directory('corr_Data/' + path)
        
        if (pos1 == pos2):
            filename = 'corr_Data/' + path + '/autocorrelation_P' + str(pos1) + '_Dict'
        else:    
            filename = 'corr_Data/' + path + '/correlation_P' + str(pos1) + '-' + str(pos2) + '_Dict'
    else:
        create_Directory('corr_Data/UnFilt')
        if (pos1 == pos2):
            filename = 'corr_Data/UnFilt/autocorrelation_P' + str(pos1) + '_Dict'
        else:    
            filename = 'corr_Data/UnFilt/correlation_P' + str(pos1) + '-' + str(pos2) + '_Dict'
    np.save(filename, my_Dict)
    
    return None

def save_Normalized_Correlation(filesave, filt_Data = True, which_Filt = 1e5):
    position = [1, 3, 5, 7, 9, 11, 13, 15]
    #####################################################################################
    """
        Calculating Norm
    """
    #####################################################################################

    #####################################################################################
    """
        Calculating Normalized Correlations; Standard Normalization
    """
    #####################################################################################
    corr_Structure = np.zeros(8)
    std_Structure = [[], [], [], [], [], [], [], []]
    #r_Structure = np.ones((4, 4))
    #t_Structure = np.ones((4, 4))
    #z_Structure = np.ones((4, 4))
    div_0 = 0
    div_1 = 0
    div_2 = 0
    div_3 = 0
    div_4 = 0
    div_5 = 0
    div_6 = 0
    div_7 = 0
    
    if filt_Data:
        file_Path = 'corr_Data/Filter_' + str(which_Filt) + '/'
    else:
        file_Path = 'corr_Data/Unfilt/'
    pos1 = 1
    for pos2 in position:
            if (pos1 <= pos2):
                #################################################
                """ Loading data """
                #################################################
                if (pos1 == pos2):
                    filename = file_Path + 'autocorrelation_P' + str(pos1) + '_Dict.npy'
                else:
                    filename = file_Path + 'correlation_P' + str(pos1) + '-' + str(pos2) + '_Dict.npy'
                my_Dict = np.load(filename, allow_pickle = True)
                """
                my_Dict = {'corr_R': corr_R, 'corr_T': corr_T, 'corr_Z': corr_Z, 'tau': tau
                    , 'corr_2R': corr_2R, 'corr_2T': corr_2T, 'corr_2Z': corr_2Z'}
                """
                #norm = standard_Norm()
                norm, norm_STD = pos_Norm(pos1, pos2, filt_Data, which_Filt)  # This should also output the standard deviation
                #print(pos1)
                #print(pos2)
                #print(norm)
                #print(norm_STD)
                
                corr_R = my_Dict.item().get('corr_R')
                print(corr_R.shape)
                corr_T = my_Dict.item().get('corr_T')
                print(corr_T.shape)
                corr_Z = my_Dict.item().get('corr_Z')
                print(corr_Z.shape)
                
                corr_2R = my_Dict.item().get('corr_2R')
                corr_2T = my_Dict.item().get('corr_2T')
                corr_2Z = my_Dict.item().get('corr_2Z')
                
                
                #######################################
                """
                    Finding the zero time delay spatial
                    correlation.
                """
                #######################################
                tau = my_Dict.item().get('tau')
                zero_Index = BMX.finding_Index_Time(tau*1e6, 0.0)
                corr_Index = int((pos2 - pos1)/2)
                #corr_Structure[corr_Index] = np.abs(np.max(corr_T))/norm
                #corr_Structure[corr_Index] = (np.max(np.abs(corr_R)) + np.max(np.abs(corr_Z)) + np.max(np.abs(corr_T)))/norm
                corr_Total = (np.abs(corr_R) + np.abs(corr_Z) + np.abs(corr_T))
                print(corr_Total[zero_Index]/norm)
                corr_Total = corr_Total[zero_Index]/norm
                ###########################################################################
                """ The uncertainty calculation """
                ###########################################################################
                
                std_R = np.sqrt(corr_2R[zero_Index] - corr_R[zero_Index]**2)
                std_T = np.sqrt(corr_2T[zero_Index] - corr_T[zero_Index]**2)
                std_Z = np.sqrt(corr_2Z[zero_Index] - corr_Z[zero_Index]**2)
                std_Total = np.sqrt(std_R**2 + std_T**2 + std_Z**2)
                print(corr_Total/norm)
                corr_Structure[corr_Index] += corr_Total
                
                if (corr_Index == 0):
                    std_Structure[0].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))   # Fractional uncertainties 
                    div_0 += 1
                elif (corr_Index == 1):
                    std_Structure[1].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_1 += 1
                elif (corr_Index == 2):
                    std_Structure[2].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_2 += 1
                elif (corr_Index == 3):
                    std_Structure[3].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_3 += 1
                elif (corr_Index == 4):
                    std_Structure[4].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_4 += 1
                elif (corr_Index == 5):
                    std_Structure[5].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_5 += 1
                elif (corr_Index == 6):
                    std_Structure[6].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_6 += 1
                elif (corr_Index == 7):
                    std_Structure[7].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_7 += 1
                    
                std_Values = [0, 0, 0, 0, 0, 0, 0, 0]
                for i in range(0, 8):
                    for values in std_Structure[i]:
                        std_Values[i] += values**2
                for i in range(0,8):
                    std_Values[i] = np.sqrt(std_Values[i])
                
                
                
                print(pos1)
                print(pos2)
                #print(corr_Structure[corr_Index])
                #print(std_Values[corr_Index])
                print([div_0, div_1, div_2, div_3, div_4, div_5, div_6, div_7])
                ############################--Previous Stuff--########################################
                #index1 = position.index(pos1)
                #index2 = position.index(pos2)
                #r_Structure[index1,index2] = np.max(np.abs(corr_R))
                #t_Structure[index1,index2] = np.max(np.abs(corr_T))
                #z_Structure[index1,index2] = np.max(np.abs(corr_Z))
                #####################################################################################
    corr_Structure[0] = corr_Structure[0]/div_0
    corr_Structure[1] = corr_Structure[1]/div_1
    corr_Structure[2] = corr_Structure[2]/div_2
    corr_Structure[3] = corr_Structure[3]/div_3
    corr_Structure[4] = corr_Structure[4]/div_4
    corr_Structure[5] = corr_Structure[5]/div_5
    corr_Structure[6] = corr_Structure[6]/div_6
    corr_Structure[7] = corr_Structure[7]/div_7
    
    std_Values[0] = std_Values[0]/div_0
    std_Values[1] = std_Values[1]/div_1
    std_Values[2] = std_Values[2]/div_2
    std_Values[3] = std_Values[3]/div_3
    std_Values[4] = std_Values[4]/div_4
    std_Values[5] = std_Values[5]/div_5
    std_Values[6] = std_Values[6]/div_6
    std_Values[7] = std_Values[7]/div_7
    
    
    print(corr_Structure)
    print(std_Values)
    separations = 2.6*np.asarray([0, 1, 2, 3, 4, 5, 6, 7])
    new_Dict = {
        'spatial_Corr': corr_Structure, 'separation': separations, 'uncertainties': np.asarray(std_Values)
    } 
    np.save(filesave, new_Dict)
    #np.savetxt("r_Correlations.csv", r_Structure, delimiter=",", header = '[1, 3, 5, 7] X [1, 3, 5, 7] Position Format')
    #np.savetxt("t_Correlations.csv", t_Structure, delimiter=",", header = '[1, 3, 5, 7] X [1, 3, 5, 7] Position Format')
    #np.savetxt("z_Correlations.csv", z_Structure, delimiter=",", header = '[1, 3, 5, 7] X [1, 3, 5, 7] Position Format')

    return None

def save_Normalized_Component_Correlation(filesave, component, filt_Data = True, which_Filt = 1e5):
    position = [1, 3, 5, 7, 9, 11, 13, 15]
    #####################################################################################
    """
        Calculating Norm
    """
    #####################################################################################

    #####################################################################################
    """
        Calculating Normalized Correlations; Standard Normalization
    """
    #####################################################################################
    corr_Structure = np.zeros(8)
    std_Structure = [[], [], [], [], [], [], [], []]
    #r_Structure = np.ones((4, 4))
    #t_Structure = np.ones((4, 4))
    #z_Structure = np.ones((4, 4))
    div_0 = 0
    div_1 = 0
    div_2 = 0
    div_3 = 0
    div_4 = 0
    div_5 = 0
    div_6 = 0
    div_7 = 0
    
    if filt_Data:
        file_Path = 'corr_Data/Filter_' + str(which_Filt) + '/'
    else:
        file_Path = 'corr_Data/Unfilt/'
    pos1 = 1
    for pos2 in position:
            if (pos1 <= pos2):
                #################################################
                """ Loading data """
                #################################################
                if (pos1 == pos2):
                    filename = file_Path + 'autocorrelation_P' + str(pos1) + '_Dict.npy'
                else:
                    filename = file_Path + 'correlation_P' + str(pos1) + '-' + str(pos2) + '_Dict.npy'
                my_Dict = np.load(filename, allow_pickle = True)
                """
                my_Dict = {'corr_R': corr_R, 'corr_T': corr_T, 'corr_Z': corr_Z, 'tau': tau
                    , 'corr_2R': corr_2R, 'corr_2T': corr_2T, 'corr_2Z': corr_2Z'}
                """
                #norm = standard_Norm()
                norm, norm_STD = pos_Component_Norm(pos1, pos2, component, filt_Data, which_Filt)  # This should also output the standard deviation
                #print(pos1)
                #print(pos2)
                #print(norm)
                #print(norm_STD)
                
                corr_R = my_Dict.item().get('corr_R')
                print(corr_R.shape)
                corr_T = my_Dict.item().get('corr_T')
                print(corr_T.shape)
                corr_Z = my_Dict.item().get('corr_Z')
                print(corr_Z.shape)
                
                corr_2R = my_Dict.item().get('corr_2R')
                corr_2T = my_Dict.item().get('corr_2T')
                corr_2Z = my_Dict.item().get('corr_2Z')
                
                
                #######################################
                """
                    Finding the zero time delay spatial
                    correlation.
                """
                #######################################
                tau = my_Dict.item().get('tau')
                zero_Index = BMX.finding_Index_Time(tau*1e6, 0.0)
                corr_Index = int((pos2 - pos1)/2)
                #corr_Structure[corr_Index] = np.abs(np.max(corr_T))/norm
                #corr_Structure[corr_Index] = (np.max(np.abs(corr_R)) + np.max(np.abs(corr_Z)) + np.max(np.abs(corr_T)))/norm
                ################# Minor Change to incorporate the component correlation ##################
                if (component.lower() == 'r'):
                    corr_Total = (np.abs(corr_R))
                elif (component.lower() == 't'):
                    corr_Total = (np.abs(corr_T))
                elif (component.lower() == 'z'):
                    corr_Total = (np.abs(corr_Z))
                print(corr_Total[zero_Index]/norm)
                corr_Total = corr_Total[zero_Index]/norm
                ###########################################################################
                """ The uncertainty calculation """
                ###########################################################################
                
                std_R = np.sqrt(corr_2R[zero_Index] - corr_R[zero_Index]**2)
                std_T = np.sqrt(corr_2T[zero_Index] - corr_T[zero_Index]**2)
                std_Z = np.sqrt(corr_2Z[zero_Index] - corr_Z[zero_Index]**2)
                if (component.lower() == 'r'):
                    std_Total = std_R
                elif (component.lower() == 't'):
                    std_Total = std_T
                elif (component.lower() == 'z'):
                    std_Total = std_Z
                    
                #std_Total = np.sqrt(std_R**2 + std_T**2 + std_Z**2)
                print(corr_Total/norm)
                corr_Structure[corr_Index] += corr_Total
                
                if (corr_Index == 0):
                    std_Structure[0].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))   # Fractional uncertainties 
                    div_0 += 1
                elif (corr_Index == 1):
                    std_Structure[1].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_1 += 1
                elif (corr_Index == 2):
                    std_Structure[2].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_2 += 1
                elif (corr_Index == 3):
                    std_Structure[3].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_3 += 1
                elif (corr_Index == 4):
                    std_Structure[4].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_4 += 1
                elif (corr_Index == 5):
                    std_Structure[5].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_5 += 1
                elif (corr_Index == 6):
                    std_Structure[6].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_6 += 1
                elif (corr_Index == 7):
                    std_Structure[7].append(corr_Total/norm*np.sqrt((norm_STD/norm**2 + std_Total/corr_Total**2)))
                    div_7 += 1
                    
                std_Values = [0, 0, 0, 0, 0, 0, 0, 0]
                for i in range(0, 8):
                    for values in std_Structure[i]:
                        std_Values[i] += values**2
                for i in range(0,8):
                    std_Values[i] = np.sqrt(std_Values[i])
                
                
                
                print(pos1)
                print(pos2)
                #print(corr_Structure[corr_Index])
                #print(std_Values[corr_Index])
                print([div_0, div_1, div_2, div_3, div_4, div_5, div_6, div_7])
                ############################--Previous Stuff--########################################
                #index1 = position.index(pos1)
                #index2 = position.index(pos2)
                #r_Structure[index1,index2] = np.max(np.abs(corr_R))
                #t_Structure[index1,index2] = np.max(np.abs(corr_T))
                #z_Structure[index1,index2] = np.max(np.abs(corr_Z))
                #####################################################################################
    corr_Structure[0] = corr_Structure[0]/div_0
    corr_Structure[1] = corr_Structure[1]/div_1
    corr_Structure[2] = corr_Structure[2]/div_2
    corr_Structure[3] = corr_Structure[3]/div_3
    corr_Structure[4] = corr_Structure[4]/div_4
    corr_Structure[5] = corr_Structure[5]/div_5
    corr_Structure[6] = corr_Structure[6]/div_6
    corr_Structure[7] = corr_Structure[7]/div_7
    
    std_Values[0] = std_Values[0]/div_0
    std_Values[1] = std_Values[1]/div_1
    std_Values[2] = std_Values[2]/div_2
    std_Values[3] = std_Values[3]/div_3
    std_Values[4] = std_Values[4]/div_4
    std_Values[5] = std_Values[5]/div_5
    std_Values[6] = std_Values[6]/div_6
    std_Values[7] = std_Values[7]/div_7
    
    
    print(corr_Structure)
    print(std_Values)
    separations = 2.6*np.asarray([0, 1, 2, 3, 4, 5, 6, 7])
    new_Dict = {
        'spatial_Corr': corr_Structure, 'separation': separations, 'uncertainties': np.asarray(std_Values)
    } 
    np.save(filesave, new_Dict)
    #np.savetxt("r_Correlations.csv", r_Structure, delimiter=",", header = '[1, 3, 5, 7] X [1, 3, 5, 7] Position Format')
    #np.savetxt("t_Correlations.csv", t_Structure, delimiter=",", header = '[1, 3, 5, 7] X [1, 3, 5, 7] Position Format')
    #np.savetxt("z_Correlations.csv", z_Structure, delimiter=",", header = '[1, 3, 5, 7] X [1, 3, 5, 7] Position Format')

    return None

filt = 2e5
pos1 = 1
for pos2 in range(1, 16, 2):
    if (pos1 <= pos2):
        two_Point_Correlations(pos1, pos2, apply_Filt = True, filt = filt)
save_Normalized_Component_Correlation('tp_rCorrelation_100kHz_P1', component = 'r', filt_Data = True, which_Filt= filt)
for pos2 in range(1, 16, 2):
    if (pos1 <= pos2):
        two_Point_Correlations(pos1, pos2, apply_Filt = False, filt = filt)
#save_Normalized_Correlation('tp_Correlation_200kHz_P1', filt_Data = True, which_Filt= filt)
save_Normalized_Component_Correlation('tp_rCorrelation_P1', component = 'r', filt_Data = False, which_Filt= filt)

filt = 1e5
"""for pos2 in range(1, 16, 2):
    if (pos1 <= pos2):
        two_Point_Correlations(pos1, pos2, apply_Filt = True, filt = filt)"""
            
#save_Normalized_Correlation('tp_Correlation_200kHz_P1', filt_Data = True, which_Filt= filt)
save_Normalized_Component_Correlation('tp_rCorrelation_200kHz_P1', component = 'r', filt_Data = True, which_Filt= filt)
filt = 4e5

"""for pos2 in range(1, 16, 2):
    if (pos1 <= pos2):
        two_Point_Correlations(pos1, pos2, apply_Filt = True, filt = filt)"""
            
#save_Normalized_Correlation('tp_Correlation_400kHz_P1', filt_Data = True, which_Filt= filt)
save_Normalized_Component_Correlation('tp_rCorrelation_400kHz_P1', component = 'r', filt_Data = True, which_Filt= filt)
filt = 3.5e5

"""for pos2 in range(1, 16, 2):
    if (pos1 <= pos2):
        two_Point_Correlations(pos1, pos2, apply_Filt = True, filt = filt)"""
            
#save_Normalized_Correlation('tp_Correlation_350kHz_P1', filt_Data = True, which_Filt= filt)
save_Normalized_Component_Correlation('tp_rCorrelation_350kHz_P1', component = 'r', filt_Data = True, which_Filt= filt)
filt = 5e5

"""for pos2 in range(1, 16, 2):
    if (pos1 <= pos2):
        two_Point_Correlations(pos1, pos2, apply_Filt = True, filt = filt)"""
            
#save_Normalized_Correlation('tp_Correlation_500kHz_P1', filt_Data = True, which_Filt= filt)
save_Normalized_Component_Correlation('tp_rCorrelation_500kHz_P1', component = 'r', filt_Data = True, which_Filt= filt)
filt = 1e6

"""for pos2 in range(1, 16, 2):
    if (pos1 <= pos2):
        two_Point_Correlations(pos1, pos2, apply_Filt = True, filt = filt)"""
            
#save_Normalized_Correlation('tp_Correlation_1MHz_P1', filt_Data = True, which_Filt= filt)            
save_Normalized_Component_Correlation('tp_rCorrelation_1MHz_P1', component = 'r', filt_Data = True, which_Filt= filt)
#save_Normalized_Component_Correlation('tp_zCorrelation_100KHz', component = 'z', filt_Data = True, which_Filt=1e5)
