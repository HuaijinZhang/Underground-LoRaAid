
#############################################################
# Set precision
# at least 150 bits for SF7, up to 4000 bits for SF12
#############################################################
import json
import logging
import gmpy2
import math
from gmpy2 import mpfr
from scipy.special import comb


gmpy2.get_context().precision = 5000

def main():
    # Set/initialize parameters
    sf = 7 # Spreading factor
    n_fft = 2**sf  # Corresponding FFT size
    distance_start = 1  # low-bound of SNR range
    distance_end = 140  # upper-bound of SNR range
    p_error = []    # list containing the error probability for Rayleigh channel
    alpha = 1.9899  # parameter for soil propagation
    beta = 75.2017  # parameter for soil propagation
    Pt = 20         # transmit power
    BW = 125000     # bandwidth
    DiversityOrder = 3   # number of nodes
    for d in range(distance_start, distance_end+1):
        distance = d*0.1
        Lug = 6.4 + 20 * math.log10(distance) + 20 * math.log10(beta) + 8.69 * alpha * distance
        snr_dB = Pt-Lug+174-10*math.log10(BW)-6
        snr1 = mpfr(n_fft*(10 ** (snr_dB / 10.0)))
        snr2 = snr1
        snr3 = snr1
        error = mpfr(0.0)  # Initialise error
        error_avg = mpfr(0.0)  # Initialise error
        for k in range(1, n_fft-1):
            nchoosek = mpfr(comb(n_fft-1, k, exact=True))
            c2 = (1+k+k*snr1)**(-1) * (1+k+k*snr2)**(-1) * (1+k+k*snr3)**(-1)
            error = error - mpfr(nchoosek * (-1)**k * (k+1)**(DiversityOrder-1) * c2)
        error_avg = mpfr(2**(sf-1)/(2**sf-1)*error)
        error_avg = mpfr(error_avg,32)  # Limit precision for printing/saving
        p_error.append(float(error_avg))
    file = open("BER_Rayleigh_sf"+str(sf)+".txt", "w")
    file.write(json.dumps(p_error))
    file.close()
    print(p_error)

if __name__ == "__main__":
    main()
