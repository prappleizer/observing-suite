:py:mod:`observing_suite.observing_plan`
========================================

.. py:module:: observing_suite.observing_plan


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   observing_suite.observing_plan.ObservingPlan




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



