:py:mod:`observing_suite.imaging`
=================================

.. py:module:: observing_suite.imaging


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   observing_suite.imaging.ScalePlot



Functions
~~~~~~~~~

.. autoapisummary::

   observing_suite.imaging.implot



.. py:function:: implot(image, figsize=(15, 13), cmap='gray_r', scale=0.5, colorbar=False, header=None, wcs=None, **kwargs)

   Plot an astronomical image, setting default options and easy tweaking of parameters

   :param image: 2D array containing an astronomical image to be plotted. Cutouts should be input as cutout.data.
   :type image: array_like
   :param figsize: figure size to use. Default: (15,13)
   :type figsize: tuple, optional
   :param cmap: Colormap to use for the image. Default: 'gray_r'
   :type cmap: str, optional
   :param scale: By default, function will scale image to some number of standard deviations about the mean pixel value. Scale sets this number (or fraction). Default: 0.5.
   :type scale: float, optional
   :param colorbar: Whether to add a colorbar or not. Default: False
   :type colorbar: bool, optional
   :param header: If input, function will attempt to create a WCS object from header and plot in celestial coordinates. Default: None
   :type header: dict, optional
   :param wcs: If input, the function will plot using a projection set by the WCS. Default: None
   :type wcs: WCS object
   :param \*\*kwargs: Additional arguments are passed to matplotlib plotting commands. Currently supported: vmin, vmax.

   :returns: figure and axes objects containing currently plotted data.
   :rtype: fig, ax


.. py:class:: ScalePlot(image)

   .. py:method:: plot(self, meanScalar=10, stdScalar=10, **kwargs)

      Function for plotting astronomical images, adding the convenience of scaling sliders and a colormap inverter.

      :param image: the image to be plotted.
      :type image: array_like
      :param kwargs: any specific keywords for the `myutils.plotting.implot()` function
      :type kwargs: dict, optional

      :returns: **fig, ax**
      :rtype: the figure and axes objects



