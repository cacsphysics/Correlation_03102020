import numpy as np


def finding_Index_Time(time_Data, time_Interest):
    ##############################################################################################
    """ This function returns the index corresponding to a given time. """
    ##############################################################################################
    """data_Path = '04232019/04232019pico1/'
    filename = data_Path + '20190423-(1).npy'
    print('Loading ' + filename)
    data_Dict = np.load(filename)"""
    """
    my_Dict = {'Time_B': timeB_Sec, 'Bz': Bz, 'Bt': Bt, 
                            'Bz_prime': Bz_prime, 'Bt_prime': Bt_prime}
    """
    time = time_Data
    output_Time = 0
    for i in range(0,time.size - 1):
        if (time[i]*1e6 <= time_Interest):
            output_Time = i

    #print('Outputing the index of interest')
    return output_Time