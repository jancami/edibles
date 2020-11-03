EDIBLES Reference
=================

Subpackages
-----------

.. toctree::

   edibles.data
   edibles.utils


edibles.models module
---------------------

.. automodule:: edibles.models
   :members:
   :undoc-members:
   :show-inheritance:


edibles.sightline module
------------------------

.. automodule:: edibles.sightline
   :members:
   :undoc-members:
   :show-inheritance:


edibles.continuum module
------------------------

The Continuum class is a tool that was developed for creating continuum models, and saving and accessing continuum data in the form of csv files. 

.. Note::
   The only method available currently is a spline fit through the spectrum to a number of anchor points.

.. Warning::
   Continuum.add_to_csv() saves the data LOCALLY in the ediblesdr4 GitHub repository. For any saved data to be reflected online, you must push your changes manually.


Example:

.. code-block:: python

   # Create the spectrun using EdiblesSpectrum
   sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")
   sp.getSpectrum(xmin=3270, xmax=3305)

   # build a 4 anchor points spline
   cont = Continuum(sp, method="spline", n_anchors=4, plot=False, verbose=2)
   # Guess the model parameters
   params = cont.model.guess(sp.flux, x=sp.wave)
   # Fit the model
   result = cont.model.fit(data=sp.flux, params=params, x=sp.wave)
   # Get the output of the fit model
   out = result.eval(params=result.params, x=sp.wave)
   # Print the result parameters
   print(result.params)
   # Plot
   plt.plot(sp.wave, sp.flux)
   plt.plot(sp.wave, out)
   plt.show()

   # Add the model parameters to a csv file
   cont.add_to_csv(
     user="Mario", comments="Test of 4 anchor points spline"
   )

   # reinitialize the edibles spectrum class, to get the continuum_filename in it
   cont = Continuum(sp)
   cont.prebuilt_model(chosen_save_num=None, plot=True, verbose=1)







.. automodule:: edibles.continuum
   :members:
   :undoc-members:
   :show-inheritance:
