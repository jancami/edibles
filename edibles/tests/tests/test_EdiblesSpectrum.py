from edibles.edibles.functions.edibles_spectrum import EdiblesSpectrum
import matplotlib.pyplot as plt


def testEdiblesSpectrum():
    sp = EdiblesSpectrum("/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits")
    print("Barycentric Velocity is", sp.v_bary)
    wave, flux = sp.getSpectrum()
    plt.plot(wave, flux)
    plt.show()


if __name__ == "__main__":
    testEdiblesSpectrum()
