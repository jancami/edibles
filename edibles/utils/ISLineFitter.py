import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class ISLineFitter():
    def __init__(self, wave, flux):
        # wave and flux from edibles spectrum or elsewhere
        self.wave = wave
        self.flux = flux

    def fit(self):
        # this will do the fitting
        pass


if __name__ == "__main__":
    print("Hello Word!")