:py:mod:`observing_suite`
=========================

.. py:module:: observing_suite


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   imaging/index.rst
   observing_log/index.rst
   observing_plan/index.rst
   target/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   observing_suite.Target
   observing_suite.ScalePlot
   observing_suite.ObservingPlan



Functions
~~~~~~~~~

.. autoapisummary::

   observing_suite.implot



.. py:class:: Target(name, parse_name=True, coordinates=None, coord_units=None)

   .. py:method:: add_configuration(self, config_name, obstype=None, coordinates=None, coord_units=None, **kwargs)

      Add an observing configuration for this target, specifying as many fields as desired.

      :param config_name: Name for this configuration. As names are eventually used in the exporting of targetlists, it is worth keeping the name short-ish, as many observatories have character limits on this column
      :type config_name: str
      :param obstype: For now, either 'imaging' or 'spectroscopy'. Some features later on depend on this.
      :type obstype: str, optional, default None
      :param coordinates: If the coordinates of this configuration differ from the object coordinates or from other configurations, supply coordinates (either SkyCoord or string). If string, coord_units must be provided.
      :type coordinates: str or SkyCoord, optional, default None
      :param coord_units: If coordinates are provided as a string, a unit (e.g., (u.hourangle, u.deg) or 'deg') as accepted by SkyCoord is required.
      :type coord_units: tuple or str, optional, default None
      :param \*\*kwargs: Any desired fields for this configuration one wants displayed later, e.g., slit pa, slit width, etc., can be added as keyword arguments with values, and will be stored.
      :type \*\*kwargs: optional

      :returns: * *None*
                * *Sets*
                * *----*
                * **self.configs** (*dict*) -- dictionary of all configuration specifications.


   .. py:method:: remove_configuration(self, config_name)

      Remove a configuration from the list

      :param config_name: the configuration name to remove
      :type config_name: str


   .. py:method:: edit_configuration(self, config_name, quantity, value)

      Edit a configuration by changing the value in one of the columns.

      :param config_name: the name of the configuration to edit
      :type config_name: str
      :param quantity: the name of the quantity (e.g., 'obstype', or a quantity added via keyword argument) to edit
      :type quantity: str
      :param value: updated value. As a note, we recommend only using this for simple string/display values. Editing, e.g., coordinates this way does not run the code to make a new SkyCoord. To change the coordinates associated with a configuration, we suggest re-adding it (with the same name) but new coords to overwrite it.
      :type value: Any


   .. py:method:: add_offset_star(self, coordinate, coord_units=None, configurations='all')

      Add an offset star to the configuration. Offset stars are used to execute blind offsets when a source is too faint to see in typical aquisition exposures.
      If an offset star is provided, the offsets between the star and the configurations coordinates (in arcsec east and north) is automatically calculated and added to the configuration.

      :param coordinate: coordinates of the offset star. Either SkyCoord object or string. If string provided, must also provide coord_units for creation of SkyCoord object.
      :type coordinate: str or SkyCoord
      :param coord_units: if coordinates provided as a string, units acceptable by SkyCoord (e.g., (u.hourangle, u.deg) or 'deg') must be provided here.
      :type coord_units: tuple or str, optional, default None
      :param configurations: Which configurations to apply this offset star to. Default is 'all', one can pass individual configuration names as strings, or a list of configuration names (as strings).
      :type configurations: str or list, optional, default 'all'

      :returns: * *None*
                * *Sets*
                * *----*
                * *Sets the 'offset star' key for the chosen configuration(s) as the star coordinates and the 'offsets' key to the offsets, visible via view_configurations().*


   .. py:method:: set_survey(self, survey_name)


   .. py:method:: retrieve_finder_chart(self, config_name, size, pixels=500, show_aperture=True, **implot_kwargs)

      Retrieve a DSS image (finder chart) around the target. If obsmode is spectroscopy, optionally show the location of the slit or circular fiber on the image.

      :param config_name: name of the configuration to retrieve finder for
      :type config_name: str
      :param size: dimensions of the finder box to use. Box is square.
      :type size: astropy Quantity
      :param pixels: dimensions (in pixels) of the image to retrieve. (Larger downloads take longer).
      :type pixels: int, optional (default 500)
      :param show_aperture: flag for whether to show an apertuer (rectangular slits and circular apertures supported). If this flag turned on, the following must be true.
                            For slits, your configuration must have properties `slit_width`, `slit_length`, and `PA`.
                            For circular apertures, your configuration must have a property `fiber_radius`.
      :type show_aperture: bool, optional (default True)
      :param \*\*implot_kwargs: arguments passed to the utility function `implot` to display the image. These include scale (images are scaled about their mean pixel value), colorbar flag, etc.
      :type \*\*implot_kwargs: optional

      :returns: **fig, ax** -- the fig and ax on which the dss image and possible aperture was plotted.
      :rtype: matplotlib figure and axes objects


   .. py:method:: add_custom_image(self, config_name, image_name, image, wcs=None)

      Add a custom image of your target. Allows for your image to be added to the observing plan along with, e.g., retrieved DSS imaging.

      :param config_name: configuration for which this image should apply. Can be a single configuration string, a list of configuration strings, or 'all'.
      :type config_name: str or list
      :param image_name: a name for the image (for later plotting and access).
      :type image_name: str
      :param image: the array containing the image
      :type image: array_like
      :param wcs: a wcs object defining the coordinates of the image. This must be provided for some functionality, like overplotting slits/apertures.
      :type wcs: astropy.WCS, optional (default None)


   .. py:method:: show_custom_image(self, config_name, image_name, show_aperture=True, **implot_kwargs)

      Display the custom image provided by user. If possible, show aperture (slit/fiber) over it.


   .. py:method:: list_configurations(self)


   .. py:method:: nudge_configuration(self, config_name, arcsec_east, arcsec_north)

      Nudge the coordinates of a configuration east or north in arcsec
      for better alignment.

      :param config_name: name of configuration to nudge
      :type config_name: str
      :param arcsec_east: amount to nudge east (west is negative) in arcsec
      :type arcsec_east: float
      :param arcsec_north: amount to nudge north (south is negative) in arcsec
      :type arcsec_north: float


   .. py:method:: configurations(self)
      :property:



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



.. py:class:: ObservingPlan(target_list, observatory, obsdates, utcoffset)

   .. py:method:: plot_visibility(self, date, target='all', view_range=12, plot_current_time=False, figsize=(30, 12), alt_min=10, alt_max=90)

      Produce a plot of altitude and airmass for targets on a given night of observing.

      :param date: date on which to make the plot. In form 'YYYY-MM-DD'. Must be a date specified within obsdate.
      :type date: str
      :param target: target to plot track for. deault is "all". Individual targets can be specified by name (str), via a list or, for single objects, a str.
      :type target: str or list, optional (default 'all')
      :param view_range: number of hours on either side of midnight over which to show the plot
      :type view_range: int, optional, default 12
      :param plot_current_time: show a vertical bar at the current time (when code executed). Useful when making plot interactively during a night.
      :type plot_current_time: bool, optional (default False)
      :param figsize: tuple specifying the figure size
      :type figsize: tuple, optional (default (15,12))


   .. py:method:: html_summary(self, date, save_dir, overwrite=True, view_range=12)

      Produce a 'beautiful' html output with the observing plan.

      :param date: date for which to construct the report. Must be a date present in the obsdates provided.
      :type date: str
      :param save_dir: location to save the observing plan. Spawns an image directory for relevant plots.
      :type save_dir: str
      :param overwrite: When set to true, constituant elements (like finder charts, airmass plots, etc) will be re-computed and saved to disk.
      :type overwrite: bool, optional (default True)



