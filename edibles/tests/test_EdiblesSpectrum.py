def testEdiblesSpectrum():

    import matplotlib.pyplot as plt
    from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum

    sp = EdiblesSpectrum("/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits")
    print("Barycentric Velocity is", sp.v_bary)
    subset = sp.getSpectrum()
    plt.plot(subset["wave"], subset["flux"])
    plt.show()


if __name__ == "__main__":
    testEdiblesSpectrum()
