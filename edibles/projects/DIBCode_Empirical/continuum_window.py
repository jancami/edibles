import numpy as np
import matplotlib.pyplot as plt
"""
Display the continuum window to ensure it is flat and 
Args:
    common_wave:
    coadded_flux:
    ContinuumMin:
Returns:
    Plot of the continuum window
"""
def continuum_window(common_wave, coadded_flux, ContinuumMin, ContinuumMax, wavelength_target):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(common_wave, coadded_flux, 'k-', lw=1, label='Coadded spectrum')
    ax.axvspan(ContinuumMin, ContinuumMax, alpha=0.3, color='green', label=f'SNR window [{ContinuumMin:.1f}–{ContinuumMax:.1f} Å]')
    ax.axvline(wavelength_target, color='red', linestyle='--', lw=1, label=f'Target {wavelength_target} Å')
    ax.set_xlabel('Wavelength [Å]')
    ax.set_ylabel('Normalized Flux')
    ax.set_title(f'SNR Window Check — {wavelength_target} Å')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
