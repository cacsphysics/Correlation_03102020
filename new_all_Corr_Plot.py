import numpy as np
import matplotlib.pylab as plt


def taylor_Scale_Fits_Plot():
    ######################################################################
    ######################################################################
    """
        Fitting Gaussian and obtaining outerscale length

        R[z] = exp[-z**2/(2*L**2)]
    """
    ######################################################################
    ######################################################################

    def parabolic_Fit(z, L):
        ###################################################
        """
            The functional form of the natural log of the
            gaussian, the form is not general but specific
            to this situation.
        """
        ###################################################
        nt_Corr = 1 - z**2/(2*L**2)

        return nt_Corr

    def linear_Fit(x, m, b):

        y = m*x + b

        return y

    def guassian_Fit(z, L):
        ###################################################
        """
            The functional form the gaussian, the form is
            not general but specific to this situation.
        """
        ###################################################
        corr = np.exp(-z**2/(2*L**2))

        return corr

    def richardson_extrapolation_v1(correlations, separations):
        taylor_List = []
        error_List = []
        max_List = []

        for i in range(1, separations.size):
            print(separations[:i+1])
            popt, pcov = curve_fit(parabolic_Fit, separations[:i+1], correlations[:i+1])

            taylor_List.append(popt[0])
            error_List.append(pcov[0])
            max_List.append(separations[i])

        taylor_Array = np.asarray(taylor_List)
        max_Sep_Array = np.asarray(max_List)
        error_Array = np.transpose(np.sqrt(np.asarray(error_List)))

        print("Returning tarylor_Array, error_Array, max_Sep_Array")

        return taylor_Array, error_Array, max_Sep_Array
    
    
    x = np.arange(0, 20, 0.1)
    plt.figure()
    plt.title('Two-Point Correlation: \n Position 1 Starting point')
    
    filename0 = '../Large_Scale/Delay_15ms/tp_Correlation_200kHz_P1.npy'
    filename1 = '../Large_Scale/Delay_15ms/tp_Correlation_350kHz_P1.npy'
    filename2 = '../Large_Scale/Delay_15ms/tp_Correlation_400kHz_P1.npy'
    filename3 = '../Large_Scale/Delay_15ms/tp_Correlation_500kHz_P1.npy'
        
    my_Dict = np.load(filename0, allow_pickle = True)
    """new_Dict = {'spatial_Corr': corr_Structure, 'separation': separations}"""
    correlations = my_Dict.item().get('spatial_Corr')
    separations = my_Dict.item().get('separation')
    plt.plot(separations, correlations, '--', label = '200kHz', color = 'blue')
    plt.plot(separations, correlations, 'o', color = 'blue')
    
    my_Dict = np.load(filename1, allow_pickle = True)
    """new_Dict = {'spatial_Corr': corr_Structure, 'separation': separations}"""
    correlations = my_Dict.item().get('spatial_Corr')
    separations = my_Dict.item().get('separation')
    plt.plot(separations, correlations, '--', label = '350kHz', color = 'green')
    plt.plot(separations, correlations, 'o', color = 'green')
    
    my_Dict = np.load(filename2, allow_pickle = True)
    """new_Dict = {'spatial_Corr': corr_Structure, 'separation': separations}"""
    correlations = my_Dict.item().get('spatial_Corr')
    separations = my_Dict.item().get('separation')
    plt.plot(separations, correlations, '--', label = '400kHz', color = 'violet')
    plt.plot(separations, correlations, 'o', color = 'violet')
    
    my_Dict = np.load(filename3, allow_pickle = True)
    """new_Dict = {'spatial_Corr': corr_Structure, 'separation': separations}"""
    correlations = my_Dict.item().get('spatial_Corr')
    separations = my_Dict.item().get('separation')
    plt.plot(separations, correlations, '--', label = '500kHz', color = 'grey')
    plt.plot(separations, correlations, 'o', color = 'grey')
    #taylor_Array, error_Array, max_Sep_Array = richardson_extrapolation_v1(correlations, separations)
    #popt, pcov = curve_fit(linear_Fit, xdata = max_Sep_Array, ydata = taylor_Array, sigma = error_Array[0])
    
    """y1 = parabolic_Fit(x, taylor_Array[0])
    y2 = parabolic_Fit(x, taylor_Array[1])
    y3 = parabolic_Fit(x, taylor_Array[2])
    y4 = parabolic_Fit(x, taylor_Array[3])
    y5 = parabolic_Fit(x, taylor_Array[4])
    y6 = parabolic_Fit(x, taylor_Array[5])
    y7 = parabolic_Fit(x, taylor_Array[6])"""


    #plt.plot(x, y1, label = '%0.2fcm'%taylor_Array[0], alpha = 0.75)
    #plt.plot(x, y2, label = '%0.2fcm'%taylor_Array[1], alpha = 0.75)
    #plt.plot(x, y3, label = '%0.2fcm'%taylor_Array[2], alpha = 0.75)
    #plt.plot(x, y4, label = '%0.2fcm'%taylor_Array[3], alpha = 0.75)
    #plt.plot(x, y5, label = '%0.2fcm'%taylor_Array[4], alpha = 0.75)
    #plt.plot(x, y6, label = '%0.2fcm'%taylor_Array[5], alpha = 0.75)
    #plt.plot(x, y7, label = '%0.2fcm'%taylor_Array[6], alpha = 0.75)
    
    x = list(np.arange(1.3, 15*1.3, 2.6))
    plt.ylim([0, 1])
    plt.xticks(x)
    plt.xlim([0, 15*1.3])
    plt.legend()
    plt.ylabel('Correlation')
    plt.xlabel('Separation (cm)')
    #plt.text(x = 9.1, y = 0.2, s = 'scales.py/taylor_Scale_Fits_Plot', color = 'grey', alpha = 0.5)
    #plt.text(x = 3.9, y = 1.1, s = '03102020 dataset: 1.5ms stuffing delay')
    #plt.legend()
    plt.show()
    plt.close()
    
    return None
taylor_Scale_Fits_Plot()