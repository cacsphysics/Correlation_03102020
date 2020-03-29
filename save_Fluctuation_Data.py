import numpy as np
import Carlos_Functions as BMX
import sys  
import os
###############################################################################################
"""
    This script will load the data and save fluctuation data. 
"""
###############################################################################################
def load_Data(pos, shot):
    
    ###################################################
    """
        Piece of code that reads in the data.
    """
    ###################################################
    filename = '03102020/port'+ str(pos) + '/20200310-(Pos' + str(pos)+ 'S'+ str(shot) + ').npy'
    
    data1 = np.load(filename, allow_pickle = True)
    time = data1.item().get('Time_B')
    #### Field Data ###
    Br1 = data1.item().get('Br')*1e4
    Bt1 = data1.item().get('Bt')*1e4
    Bz1 = data1.item().get('Bz')*1e4
    #Br1 = process_Data(Br1, time, pos = pos)
    #Bz1 = process_Data(Bz1, time, pos = pos)
    #Bmag1 = data1.item().get('Bmag')
    print('Br, Bt, Bz, time')
    return Br1, Bt1, Bz1, time

def ensemble_Average(pos):
    #####################################################################
    """
        This functions calculates and saves the ensemble average 
    """
    #####################################################################
    avoid_Shots = []
    all_Shots = list(range(1,26))
    #print(all_Shots[-1])
    #shots_1ms = list(range(7,27))
    #shots_10ms = list(range(27,48))
    #print(all_Shots)
    considered_Shots = list(set(all_Shots).symmetric_difference(set(avoid_Shots)))
    for shot in considered_Shots:
        Br_Hold, Bt_Hold, Bz_Hold, time = load_Data(pos, shot)
        #print(Br_Hold)
        if shot == considered_Shots[0]:
            Br = Br_Hold/(len(considered_Shots) - 1)
            Bt = Bt_Hold/(len(considered_Shots) - 1)
            Bz = Bz_Hold/(len(considered_Shots) - 1)
        else:
            Br += Br_Hold/(len(considered_Shots) - 1)
            Bt += Bt_Hold/(len(considered_Shots) - 1)
            Bz += Bz_Hold/(len(considered_Shots) - 1)
        
    ensemble_Z = Bz
    ensemble_T = Bt
    ensemble_R = Br
    print('[ensemble_R, ensemble_T, ensemble_Z]')
    np.save('ensemble_Pos%i'%int(pos), [ensemble_R, ensemble_T, ensemble_Z, time])
    return None

def save_Fluctuation_Data(pos, shot):
    ###########################################################################################
    """
        This function will save the data subtracted by the ensemble average.
    """
    ###########################################################################################
    ###########################################################################################
    """ Checking if ensemble save file exists """
    ###########################################################################################
    if (os.path.exists('ensemble_Pos%i.npy'%int(pos))):
        ens_Data = np.load('ensemble_Pos%i.npy'%int(pos), allow_pickle = True)
    else:
        ensemble_Average(pos)
        ens_Data = np.load('ensemble_Pos%i.npy'%int(pos), allow_pickle = True)
    ense_R = ens_Data[0]
    ense_T = ens_Data[1]
    ense_Z = ens_Data[2]
    Br, Bt, Bz, time = load_Data(pos, shot)
    
    Br_fluct = Br - ense_R
    Bt_fluct = Bt - ense_T
    Bz_fluct = Bz - ense_Z
    ###########################################################################################
    """ Checking if directory exists """
    ###########################################################################################
    if (os.path.isdir('Fluctuation_Data')):
        if (os.path.isdir('Fluctuation_Data/port' + str(pos))):
            filename = 'Fluctuation_Data/port' + str(pos) + '/20200310-(Pos' + str(pos) + 'S' + str(shot) + ')-fluctuation'
        else:
            os.mkdir('Fluctuation_Data/port' + str(pos))
            filename = 'Fluctuation_Data/port' + str(pos) + '/20200310-(Pos' + str(pos) + 'S' + str(shot) + ')-fluctuation'
    else:
        os.mkdir('Fluctuation_Data')
        if (os.path.isdir('Fluctuation_Data/port' + str(pos))):
            filename = 'Fluctuation_Data/port' + str(pos) + '/20200310-(Pos' + str(pos) + 'S' + str(shot) + ')-fluctuation'
        else:
            os.mkdir('Fluctuation_Data/port' + str(pos))
            filename = 'Fluctuation_Data/port' + str(pos) + '/20200310-(Pos' + str(pos) + 'S' + str(shot) + ')-fluctuation'
    my_Dict = {'br': Br_fluct, 'bt': Bt_fluct, 'bz': Bz_fluct,'time': time}
    
    np.save(filename, my_Dict)
    
    return None


###############################################################################################
""" Functions Above this Line """
###############################################################################################


for shot in range(1, 26):
    for pos in range(1, 16, 2):
        save_Fluctuation_Data(pos = pos, shot = shot)
