from astroquery.vizier import Vizier
import numpy as np
import scipy.optimize as op



catalog_list = Vizier.find_catalogs('New white dwarf SDSS DR12')
print({k:v.description for k,v in catalog_list.items()})

Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs(list(catalog_list))


# Our model has a temperature effect, then absorption features on top of the model


# Constants
h_planck = 4.136e-15 # in eV*s
c_light = 2.998e10   # in cgs
k_boltz = 8.617e-5   # in eV/K


# Balmer line centers (in cgs):
lambda_balmer = np.array([6563, 4861, 4341, 4102, 3970, 3889])
window_dict = {}



# Lorentz distribution:
def lorentz(C_coeff, x_center, gamma, x):
    return -C_coeff / (np.pi * gamma * (1.0 + ((x-x_center)/gamma)**2))

#We must define our probability function in order to use emcee
def like(p, x, y):
    C_coeff, x_center, gamma = p
    model = lorentz(C_coeff, x_center, gamma, x)
    sigma = (model - y)**2
    return np.sum(sigma)

def my_line(points, x_vec):
    x_coords, y_coords = zip(*points)
    A = np.vstack([x_coords,np.ones(len(x_coords))]).T
    m, c = np.linalg.lstsq(A, y_coords)[0]
    return m*x_vec + c

def FWHM(height,norm):
    gam = 1/(np.pi*np.abs(height))
    return gam*norm



nll = lambda *args: -like(*args)

def spectrum_fit(lambda_balmer,all_lambda,all_flux,num):
    line = lambda_balmer[num]
    window = np.where(np.logical_and(all_lambda>= (line - 30), all_lambda<= (line + 30)))[0]

    #Get an educated guess for C_true
    thirds = len(all_flux[window])/3
    left_max = np.max(all_flux[window][:thirds])
    right_max = np.max(all_flux[window][-thirds:])
    left_point = all_lambda[window][:thirds][np.argmax(all_flux[window][:thirds])]
    right_point = all_lambda[window][2*thirds + np.argmax(all_flux[window][-thirds:])]

    points = [(left_point,left_max),(right_point,right_max)]
    vert_fix = my_line(points,all_lambda[window])
    int_y = all_flux[window]-vert_fix
    result = np.trapz(int_y)

    left_mean = np.mean(all_flux[window][:thirds])
    right_mean = np.mean(all_flux[window][-thirds:])
    left_lambda = np.mean(all_lambda[window][:thirds])
    right_lambda = np.mean(all_lambda[window][-thirds:])

    points = [(left_lambda,left_mean),(right_lambda,right_mean)]
    vert_shift = my_line(points,all_lambda[window])

    #Define initial guesses
    C_true = result*-1
    x_true = (all_lambda[window][0] + all_lambda[window[-1]])/2
    gamma_true = FWHM(np.min(all_flux[window]-vert_shift),C_true)
    print(gamma_true)

    result = op.minimize(nll, [C_true, x_true, gamma_true], args=(all_lambda[window],all_flux[window]))
    C_coeff_ls, x_center_ls, gamma_ls = result["x"]

    return C_coeff_ls, x_center_ls, gamma_ls, window, vert_fix, vert_shift




#Main loop
with open('SDSS_lorentz_params.txt', 'w') as lorentz_params:
    lorentz_params.write('Plate C1 X1 G1 C2 X2 G2 C3 X3 G3 C4 X4 G4 C5 X5 G5 C6 X6 G6 \n')
# lorentz_params = open('SDSS_lorentz_params.txt', 'w')
# lorentz_params.write('Plate C1 X1 G1 C2 X2 G2 C3 X3 G3 C4 X4 G4 C5 X5 G5 C6 X6 G6 \n')

C_coeff_ls = np.zeros(6)
x_center_ls = np.zeros(6)
gamma_ls = np.zeros(6)
total = 1

for i in range(len(catalogs[0]['SpType'])):

    if 'DA' in catalogs[0]['SpType'][i].decode('utf8'):

        plate = catalogs[0]['PMF'][i].decode('utf8')

        #Gets data from file
        directory = "../data/"
        filename = directory + "spectra-"+str(plate)+".npy"

        try:
            data = np.load(filename)
        except:
            print("Missing", filename)
            continue

        all_lambda = data['wavelength']
        all_flux = data['flux']

        with open('SDSS_lorentz_params.txt', 'a') as lorentz_params:
            lorentz_params.write(str(total) + ' ' + str(plate))

        total = total + 1
        if total%100 == 0:  print(total)

        for ind,line in enumerate(lambda_balmer):
            #Make the windows for lorentzian fitting

            window = np.where(np.logical_and(all_lambda>= (line - 30), all_lambda<= (line + 30)))[0]

            #Get an educated guess for C_true
            if not all_flux[window].size:
                C_coeff_ls[ind] = 0.0
                x_center_ls[ind] = 0.0
                gamma_ls[ind] = 0.0

                with open('SDSS_lorentz_params.txt', 'a') as lorentz_params:
                    lorentz_params.write(' ' + str(C_coeff_ls[ind]) + ' ' + str(x_center_ls[ind]) + ' ' + str(-gamma_ls[ind]))
                continue

            thirds = int(len(all_flux[window])/3)
            left_max = np.max(all_flux[window][:thirds])
            right_max = np.max(all_flux[window][-thirds:])
            left_point = all_lambda[window][:thirds][np.argmax(all_flux[window][:thirds])]
            right_point = all_lambda[window][2*thirds + np.argmax(all_flux[window][-thirds:])]

            points = [(left_point,left_max),(right_point,right_max)]
            vert_fix = my_line(points,all_lambda[window])
            int_y = all_flux[window]-vert_fix
            result = np.trapz(int_y)

            left_mean = np.mean(all_flux[window][:thirds])
            right_mean = np.mean(all_flux[window][-thirds:])
            left_lambda = np.mean(all_lambda[window][:thirds])
            right_lambda = np.mean(all_lambda[window][-thirds:])

            points = [(left_lambda,left_mean),(right_lambda,right_mean)]
            vert_shift = my_line(points,all_lambda[window])

            #Define initial guesses
            C_true = result*-1
            x_true = (all_lambda[window][0] + all_lambda[window[-1]])/2
            gamma_true = FWHM(np.min(all_flux[window]-vert_shift),C_true)
#             print gamma_true

            result = op.minimize(nll, [C_true, x_true, gamma_true], args=(all_lambda[window],all_flux[window]))
            C_coeff_ls[ind], x_center_ls[ind], gamma_ls[ind] = result["x"]

            with open('SDSS_lorentz_params.txt', 'a') as lorentz_params:
                lorentz_params.write(' ' + str(C_coeff_ls[ind]) + ' ' + str(x_center_ls[ind]) + ' ' + str(gamma_ls[ind]))


        with open('SDSS_lorentz_params.txt', 'a') as lorentz_params:
            lorentz_params.write('\n')
