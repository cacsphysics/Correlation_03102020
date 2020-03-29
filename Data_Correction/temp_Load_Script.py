import numpy as np
import Carlos_Functions as BMX
import scipy.integrate as sp
import scipy.signal as sps

"""
This script is a temporary holding place for the code. Essentially the previous code was not general
enough. 
"""
#################################################################
""" Below are some inductive-probe parameters """
#################################################################
probe_Dia = 0.003175
#ztprobe_Dia = 0.00158755
hole_Sep = 0.001016

rprobe_Area = np.pi*(probe_Dia/2)**2
ztprobe_Area = np.pi*probe_Dia*hole_Sep/4
#################################################################
def which_Indicator(pico_Num, shot_Num, starting_Pos):
    ####################################################
    """
        Quick function to generate the indicator
    """
    ####################################################
    return '(Pos' + str(starting_Pos + 2*(pico_Num - 1)) + 'S' + str(shot_Num) + ')'

def correct_Integration(input_Data, time, data_Component):
    
    if (data_Component == 'z' or data_Component == 't'):
        output_Data = sp.cumtrapz(input_Data/ztprobe_Area, time*1e-6)
    elif (data_Component == 'r'):
        output_Data = sp.cumtrapz(input_Data/rprobe_Area, time*1e-6)
    else:
        output_Data = input_Data
    return output_Data

def averaging_Data(input_Data):
    
    avg_Data = np.zeros(input_Data.size - 1)
    for i in range(0, input_Data.size):
        avg_Data[i] = (input_Data[i + 1] + input_Data[i])/2
        return avg_Data
    
def correct_Time(time):
    
    time_B = np.zeros(time.size - 1)
    for i in range(0, time_B.size):
        time_B[i] = (time[i+1] + time[i])/2
    return time_B

def mean_Subtraction(data, time):
    
    ending_Index = BMX.finding_Index_Time(time, 0)
    data -= np.mean(data[:ending_Index])
    return data

def correct_Area(input_Data, data_Component):
    
    if (data_Component == 'z' or data_Component == 't'):
        output_Data = input_Data/ztprobe_Area
    elif (data_Component == 'r'):
        output_Data = input_Data/rprobe_Area
    else:
        output_Data = input_Data
    return output_Data

shot_Max = 25
pico_Max = 8

data_Path = '03102020/'
for shot in range(1, shot_Max + 1):
    for pico in range(1, pico_Max + 1):
        filename = data_Path + '03102020pico' + str(pico) + '/20200310-0001 (' + str(shot) + ').txt'
        data = BMX.BMX_Pico_Read(filename = filename)
        ###################################################
        """ The picoscope-probe structure """
        ###################################################
        """         A           B           C           D
        pico1       1R          1T          1Z          3R
        pico2       3T          3Z          5R          5T
        pico3       5Z          7R          7T          7Z
        pico4       --          N           DC          HV
        N --- Density
        DC --- Discharge Current
        HV --- High Voltage
        """
        ####################################################
        
        data = BMX.BMX_Pico_Read(filename)
            
        time = data[0]
        
        end_Index = BMX.finding_Index_Time(time*1e-6, 13)  # Avoiding the EMP from initial discharge
        
        Br_dot = mean_Subtraction(data[1], time)
        Bt_dot = mean_Subtraction(data[2], time)
        Bz_dot = mean_Subtraction(data[3], time)
        my_Port = {'Time': time, 'Br_dot': Br_dot, 'Bt_dot': Bt_dot, 'Bz_dot': Bz_dot}
        Br = correct_Integration(Br_dot[end_Index:], time[end_Index:], 'r')
        Bt = correct_Integration(Bt_dot[end_Index:], time[end_Index:], 't')
        Bz = correct_Integration(Bz_dot[end_Index:], time[end_Index:], 'z')
        
        Br_dot = correct_Area(Br_dot, 'r')
        Bt_dot = correct_Area(Bt_dot, 't')
        Bz_dot = correct_Area(Bz_dot, 'z')
        
        time_B = correct_Time(time[end_Index:])
        
        my_Port['Br'] = Br
        my_Port['Bt'] = Bt
        my_Port['Bz'] = Bz
        my_Port['Time_B'] = time_B
        
        if (pico == 1):
            DC = averaging_Data(data[4]) # Discharge Current
            my_Volt = {'Current': DC, 'Time_B': time_B}    
        elif (pico == 8):
            HV = averaging_Data(data[4]) # High Voltage
            my_Volt['Volt'] = HV
        indicator1 = 1 + 2*(pico - 1)
        savefile = data_Path + 'port' + str(indicator1) + '/20200310-' + which_Indicator(pico_Num = pico, shot_Num = shot, starting_Pos = 1)
        np.save(savefile, my_Port)
    savefileHV = data_Path + 'diagnostic/20200310-(DenCurrVoltS' + str(shot) + ')'
    np.save(savefileHV, my_Volt)

