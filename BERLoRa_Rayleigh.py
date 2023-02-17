
#############################################################
# Set precision
# at least 150 bits for SF7, up to 4000 bits for SF12
#############################################################
import matplotlib.pyplot as plt
import numpy as np
import json
import logging
import gmpy2
import math
from gmpy2 import mpfr
from scipy.special import comb


gmpy2.get_context().precision = 5000

def main():
    # Set/initialize parameters
    sf = [7, 10, 12]      # Spreading factor
    n_fft_sf7 = 2**sf[0]  # Corresponding FFT size
    n_fft_sf10 = 2 ** sf[1]  # Corresponding FFT size
    n_fft_sf12 = 2 ** sf[2]  # Corresponding FFT size
    distance_start = 1  # low-bound of SNR range
    distance_end = 100  # upper-bound of SNR range
    p_error_sf7 = []    # list containing the error probability for Rayleigh channel
    p_error_sf10 = []
    p_error_sf12 = []
    Pt = 20         # transmit power
    BW = 125000     # bandwidth
    DiversityOrder = 3   # number of nodes
    mv = 0.15       # water content
    S = 0.9         # Proportion of sand
    C = 0.1         # Proportion of clay

    # Calculate alpha and beta
    def UndergroundPathLoss(mv, S, C):
        f = 900 * pow(10, 6)
        ew_inf = 4.9
        ew0 = 80.1
        e_air = 8.854 * pow(10, -12)
        rho_s = 2.65
        rho_b = 1.5
        pi = 3.1415926
        delta_eff = 0.0467 + 0.2204 * rho_b - 0.4111 * S + 0.6614 * C
        e_water_real = ew_inf + (ew0 - ew_inf) / (1 + pow((0.58 * pow(10, -10) * f), 2))
        e_water_imag = 0.58 * pow(10, -10) * f * (ew0 - ew_inf) / (
                    1 + pow((0.58 * pow(10, -10) * f), 2)) + delta_eff / (2 * pi * f * e_air) * (rho_s - rho_b) / (
                                   rho_s * mv)
        e0 = 1.15 * pow((1 + rho_b / rho_s * (pow((pow((1.01 + 0.44 * 2.65), 2) - 0.062), 0.65) - 1) + pow(mv, (
                    1.2748 - 0.519 * S - 0.152 * C)) * pow(e_water_real, 0.65) - mv), (1 / 0.65)) - 0.68
        e1 = pow(pow(mv, (1.33797 - 0.603 * S - 0.166 * C)) * pow(e_water_imag, 0.65), (1 / 0.65))
        e0 = e0 * e_air
        e1 = e1 * e_air

        mu = 1.0006 * 4 * pi * pow(10, -7)
        alpha = 2 * pi * f * math.sqrt((mu * e0 / 2) * (math.sqrt(1 + pow((e1 / e0), 2)) - 1))
        beta = 2 * pi * f * math.sqrt((mu * e0 / 2) * (math.sqrt(1 + pow((e1 / e0), 2)) + 1))
        lamda = 2 * pi / beta
        return alpha, beta


    for d in range(distance_start, distance_end+1):
        distance = d*0.1
        print(d/distance_end*100, "%")
        alpha, beta = UndergroundPathLoss(mv, S, C) # parameter for soil propagation
        Lug = 6.4 + 20 * math.log10(distance) + 20 * math.log10(beta) + 8.69 * alpha * distance
        snr_dB = Pt-Lug+174-10*math.log10(BW)-6
        # Get the SNR for each SF
        snr1_sf7 = mpfr(n_fft_sf7*(10 ** (snr_dB / 10.0)))
        snr2_sf7 = snr1_sf7
        snr3_sf7 = snr1_sf7
        snr1_sf10 = mpfr(n_fft_sf10 * (10 ** (snr_dB / 10.0)))
        snr2_sf10 = snr1_sf10
        snr3_sf10 = snr1_sf10
        snr1_sf12 = mpfr(n_fft_sf12 * (10 ** (snr_dB / 10.0)))
        snr2_sf12 = snr1_sf12
        snr3_sf12 = snr1_sf12
        # Initialise error
        error_sf7 = mpfr(0.0)
        error_avg_sf7 = mpfr(0.0)
        error_sf10 = mpfr(0.0)
        error_avg_sf10 = mpfr(0.0)
        error_sf12 = mpfr(0.0)
        error_avg_sf12 = mpfr(0.0)

        # SF=7
        for k in range(1, n_fft_sf7-1):
            nchoosek = mpfr(comb(n_fft_sf7 - 1, k, exact=True))
            c2 = (1 + k + k * snr1_sf7) ** (-1) * (1 + k + k * snr2_sf7) ** (-1) * (1 + k + k * snr3_sf7) ** (-1)
            error_sf7 = error_sf7 - mpfr(nchoosek * (-1) ** k * (k + 1) ** (DiversityOrder - 1) * c2)
        # SF=10
        for k in range(1, n_fft_sf10 - 1):
            nchoosek = mpfr(comb(n_fft_sf10 - 1, k, exact=True))
            c2 = (1 + k + k * snr1_sf10) ** (-1) * (1 + k + k * snr2_sf10) ** (-1) * (1 + k + k * snr3_sf10) ** (-1)
            error_sf10 = error_sf10 - mpfr(nchoosek * (-1) ** k * (k + 1) ** (DiversityOrder - 1) * c2)
        # SF=12
        for k in range(1, n_fft_sf12 - 1):
            nchoosek = mpfr(comb(n_fft_sf12 - 1, k, exact=True))
            c2 = (1 + k + k * snr1_sf12) ** (-1) * (1 + k + k * snr2_sf12) ** (-1) * (1 + k + k * snr3_sf12) ** (-1)
            error_sf12 = error_sf12 - mpfr(nchoosek * (-1) ** k * (k + 1) ** (DiversityOrder - 1) * c2)

        # Save the error probability of SF=7
        ErrorAvg_sf7 = mpfr(2**(sf[0]-1)/(2**sf[0]-1)*error_sf7)
        ErrorAvg_sf7 = mpfr(ErrorAvg_sf7, 32)  # Limit precision for printing/saving
        p_error_sf7.append(float(ErrorAvg_sf7))
        # Save the error probability of SF=10
        ErrorAvg_sf10 = mpfr(2 ** (sf[1] - 1) / (2 ** sf[1] - 1) * error_sf10)
        ErrorAvg_sf10 = mpfr(ErrorAvg_sf10, 32)  # Limit precision for printing/saving
        p_error_sf10.append(float(ErrorAvg_sf10))
        # Save the error probability of SF=12
        ErrorAvg_sf12 = mpfr(2 ** (sf[2] - 1) / (2 ** sf[2] - 1) * error_sf12)
        ErrorAvg_sf12 = mpfr(ErrorAvg_sf12, 32)  # Limit precision for printing/saving
        p_error_sf12.append(float(ErrorAvg_sf12))

    plot_error_sf7 = p_error_sf7
    plot_error_sf10 = p_error_sf10
    plot_error_sf12 = p_error_sf12
    # Save to file
    file = open("BER_Rayleigh_sf"+str(sf[0])+".txt", "w")
    file.write(json.dumps(p_error_sf7))
    file.close()
    file = open("BER_Rayleigh_sf" + str(sf[1]) + ".txt", "w")
    file.write(json.dumps(p_error_sf10))
    file.close()
    file = open("BER_Rayleigh_sf" + str(sf[2]) + ".txt", "w")
    file.write(json.dumps(p_error_sf12))
    file.close()
    # Plot figure
    x = np.linspace(0.1, 10, 100)
    fig, ax = plt.subplots()
    ax.semilogy(x, plot_error_sf7)
    ax.semilogy(x, plot_error_sf10)
    ax.semilogy(x, plot_error_sf12)
    plt.show()

if __name__ == "__main__":
    main()





