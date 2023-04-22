
from ast import Continue
from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import pandas as pd
from photutils.aperture import SkyRectangularAperture, SkyCircularAperture
from .imaging import implot 
from astroquery.skyview import SkyView
import legacystamps
__all__ = ['Target']


class Target():
  def __init__(self,name,parse_name=True,coordinates=None,coord_units=None):
    '''
    Initializes the Target object. 
    
    Parameters
    ----------
    name: str
      name to use for the target throughout the suite. 
    parse_name: bool, optional, default True
      If the name is that of a known object resolvable by Simbad, parse it to determine coordinates. 
    coordinates: SkyCoord or str, optional, default None
      If parse_name is False, supply coordinates for the target manually. Must be a SkyCoord object or string with coordinates. If string, coord_units must be supplied. 
    coord_units: tuple or str, optional, default None
      if supplying coordinates as a string, the units as accepted by SkyCoord must be provided, e.g., (u.hourangle,u.deg) or 'deg'. 
    
    Returns
    -------
    None 
    
    Sets
    ----
    configs: dict
      a dictionary containing configuration information.
    
    Notes
    -----
    It is not strictly necessary for the Target itself to have coordinates defined, but every configuration must. 
    '''
    self.name = name
    if parse_name:
      self.coordinates = SkyCoord.from_name(name)
    else:
      if coordinates is not None:
        if isinstance(coordinates,str):
          if coord_units is None:
            raise AssertionError('When providing string coordinates, a coordinate units accepted by SkyCoord are required to be passed to coord_units')
          else:
            self.coordinates =  SkyCoord(coordinates,unit=coord_units)
        elif isinstance(coordinates,SkyCoord):
          self.coordinates = coordinates           
    self.configs = {}
    
  
  def add_configuration(self,config_name,obstype=None,coordinates=None,coord_units=None,**kwargs):
    '''
    Add an observing configuration for this target, specifying as many fields as desired.
    
    Parameters
    ----------
    config_name: str
      Name for this configuration. As names are eventually used in the exporting of targetlists, it is worth keeping the name short-ish, as many observatories have character limits on this column
    obstype: str, optional, default None
      For now, either 'imaging' or 'spectroscopy'. Some features later on depend on this. 
    coordinates: str or SkyCoord, optional, default None
      If the coordinates of this configuration differ from the object coordinates or from other configurations, supply coordinates (either SkyCoord or string). If string, coord_units must be provided. 
    coord_units: tuple or str, optional, default None
      If coordinates are provided as a string, a unit (e.g., (u.hourangle, u.deg) or 'deg') as accepted by SkyCoord is required. 
    **kwargs: optional 
      Any desired fields for this configuration one wants displayed later, e.g., slit pa, slit width, etc., can be added as keyword arguments with values, and will be stored.
    
    Returns
    -------
    None
    
    Sets
    ----
    self.configs: dict
      dictionary of all configuration specifications.
    '''
    if config_name in self.configs.keys():
        cont = input(f'Config Name {config_name} already a configuration. Overwrite? [Enter yes, N for no]: ')
        if  cont.upper() == 'N':
            return 
            
    self.configs[config_name] = {}
    self.configs[config_name]['obstype']= obstype
    if coordinates is not None:
      if isinstance(coordinates,SkyCoord):
        self.configs[config_name]['coordinates'] = coordinates
      elif isinstance(coordinates,str):
        if coord_units is None:
          raise AssertionError('When providing string coordinates, a coordinate units accepted by SkyCoord are required to be passed to coord_units')
        else:
          self.configs[config_name]['coordinates'] = SkyCoord(coordinates,unit=coord_units)   
    elif self.coordinates is not None:
      self.configs[config_name]['coordinates'] = self.coordinates
    else:
      self.configs[config_name]['coordinates'] = None
    for i in kwargs.keys():
      self.configs[config_name][i] = kwargs[i]

  def remove_configuration(self,config_name):
    '''
    Remove a configuration from the list
    
    Parameters
    ----------
    config_name: str
      the configuration name to remove
    '''
    try:
      self.configs.pop(config_name)
    except KeyError:
      print('config not found')
      return
  def edit_configuration(self,config_name,quantity,value):
    '''
    Edit a configuration by changing the value in one of the columns.
    
    Parameters
    ----------
    config_name: str
      the name of the configuration to edit
    quantity: str
      the name of the quantity (e.g., 'obstype', or a quantity added via keyword argument) to edit
    value: Any
      updated value. As a note, we recommend only using this for simple string/display values. Editing, e.g., coordinates this way does not run the code to make a new SkyCoord. To change the coordinates associated with a configuration, we suggest re-adding it (with the same name) but new coords to overwrite it.
    '''
    try:
      self.configs[config_name][quantity] = value
    except KeyError:
      print('configuration name not found')
      return
  

  def add_offset_star(self,coordinate,coord_units=None,configurations='all'):
    '''
    Add an offset star to the configuration. Offset stars are used to execute blind offsets when a source is too faint to see in typical aquisition exposures.
    If an offset star is provided, the offsets between the star and the configurations coordinates (in arcsec east and north) is automatically calculated and added to the configuration.
    
    Parameters
    ----------
    coordinate: str or SkyCoord
      coordinates of the offset star. Either SkyCoord object or string. If string provided, must also provide coord_units for creation of SkyCoord object.
    coord_units: tuple or str, optional, default None
      if coordinates provided as a string, units acceptable by SkyCoord (e.g., (u.hourangle, u.deg) or 'deg') must be provided here. 
    configurations: str or list, optional, default 'all'
      Which configurations to apply this offset star to. Default is 'all', one can pass individual configuration names as strings, or a list of configuration names (as strings).
      
    Returns
    -------
    None
    
    Sets
    ----
    Sets the 'offset star' key for the chosen configuration(s) as the star coordinates and the 'offsets' key to the offsets, visible via view_configurations().
    '''
    if isinstance(coordinate,str):
      if coord_units is not None:
        coord = SkyCoord(coordinate,unit=coord_units)
      else:
        raise AssertionError('If string coordinate provided, units must be provided for SkyCoord creation')
    elif isinstance(coordinate,SkyCoord):
      coord = coordinate
    if configurations=='all':
      for i in self.configs.keys():
        os = coord.spherical_offsets_to(self.configs[i]['coordinates'])
        os = [os[0].to(u.arcsec).value,os[1].to(u.arcsec).value]
        add_str = f'''{os[0]:.3f}'' E, {os[1]:.3f}'' N'''
        self.configs[i]['offset star'] = coord
        self.configs[i]['offsets'] = add_str
    elif isinstance(configurations,str):
      os = coord.spherical_offsets_to(self.configs[configurations]['coordinates'])
      os = [os[0].to(u.arcsec).value,os[1].to(u.arcsec).value]
      add_str = f'''{os[0]:.3f}'' E, {os[1]:.3f}'' N'''
      self.configs[configurations]['offset star'] = coord
      self.configs[configurations]['offsets'] = add_str
    elif isinstance(configurations,list):
      for i in configurations:
        os = coord.spherical_offsets_to(self.configs[i]['coordinates'])
        os = [os[0].to(u.arcsec).value,os[1].to(u.arcsec).value]
        add_str = f'''{os[0]:.3f}'' E, {os[1]:.3f}'' N'''
        self.configs[i]['offset star'] = coord
        self.configs[i]['offsets'] = add_str
        
  def set_survey(self,survey_name):
    self.survey_name = survey_name
  
  def retrieve_finder_chart(self,config_name,size,pixels=500,show_aperture=True,**implot_kwargs):
    '''
    Retrieve a DSS image (finder chart) around the target. If obsmode is spectroscopy, optionally show the location of the slit or circular fiber on the image.
    
    Parameters
    ----------
    config_name: str
      name of the configuration to retrieve finder for 
    size: astropy Quantity
      dimensions of the finder box to use. Box is square. 
    pixels: int, optional (default 500)
      dimensions (in pixels) of the image to retrieve. (Larger downloads take longer).
    show_aperture: bool, optional (default True)
      flag for whether to show an apertuer (rectangular slits and circular apertures supported). If this flag turned on, the following must be true. 
      For slits, your configuration must have properties `slit_width`, `slit_length`, and `PA`. 
      For circular apertures, your configuration must have a property `fiber_radius`. 
    **implot_kwargs: optional
      arguments passed to the utility function `implot` to display the image. These include scale (images are scaled about their mean pixel value), colorbar flag, etc. 
    
    Returns
    -------
    fig, ax: matplotlib figure and axes objects
      the fig and ax on which the dss image and possible aperture was plotted.
    '''
    
    if hasattr(self,'survey_name'):
      survey=self.survey_name
    else:
      survey='SDSSdr7g'
    if survey != 'legacy':
      sv = SkyView()
      paths = sv.get_images(position=self.configs[config_name]['coordinates'],
                        survey=[survey],
                        coordinates='J2000',
                        width=size,
                        height=size,
                        grid=True,
                        gridlabels=True,
                        pixels=str(pixels))
      image = paths[0][0].data
      wcs = WCS(paths[0][0].header)
    elif survey == 'legacy':
      ra = self.configs[config_name]['coordinates'].ra.value
      dec = self.configs[config_name]['coordinates'].dec.value
      with fits.open(f"https://www.legacysurvey.org/viewer/cutout=fits?&ra={ra}&dec={dec}&layer=ls-dr9&zoom=14") as hdu:
        image = hdu[0].data[0]
        wcs = WCS(hdu[0].header)[0]
    fig, ax = implot(image,wcs=wcs,cmap='gray',**implot_kwargs)
    if show_aperture:
      if self.configs[config_name].keys() >= {'slit_width','slit_length','PA'}:
        slit = SkyRectangularAperture(self.configs[config_name]['coordinates'],
                                    w=self.configs[config_name]['slit_width'],
                                    h=self.configs[config_name]['slit_length'],
                                    theta=self.configs[config_name]['PA']+90*u.deg)
        slit.to_pixel(wcs).plot(color='r',lw=3)
      elif self.configs[config_name].keys() >= {'fiber_radius'}:
        fiber = SkyCircularAperture(self.configs[config_name]['coordinates'],
                                    r=self.configs[config_name]['fiber_radius'])
        fiber.to_pixel(wcs).plot(color='r',lw=3)
      else:
        raise KeyError('''show_slit set to true, but this configuration does not have 'slit_width','slit_length', and 'PA' set, which are needed for slit display, or 'fiber_radius' set, for circular aperture.''')
    return fig, ax          
    
  def add_custom_image(self,config_name,image_name,image,wcs=None):
    '''
    Add a custom image of your target. Allows for your image to be added to the observing plan along with, e.g., retrieved DSS imaging.
    
    Parameters
    ----------
    config_name: str or list
      configuration for which this image should apply. Can be a single configuration string, a list of configuration strings, or 'all'. 
    image_name: str
      a name for the image (for later plotting and access). 
    image: array_like
      the array containing the image 
    wcs: astropy.WCS, optional (default None)
      a wcs object defining the coordinates of the image. This must be provided for some functionality, like overplotting slits/apertures.
    '''
    self.configs[config_name]['user_images'] = {}
    self.configs[config_name]['user_images'][image_name] = {}
    self.configs[config_name]['user_images'][image_name]['image'] = image
    self.configs[config_name]['user_images'][image_name]['wcs'] = wcs
                              
  def show_custom_image(self,config_name,image_name,show_aperture=True,**implot_kwargs):
    '''
    Display the custom image provided by user. If possible, show aperture (slit/fiber) over it.
    '''
    image = self.configs[config_name]['user_images'][image_name]['image']
    wcs = self.configs[config_name]['user_images'][image_name]['wcs']
    fig, ax = implot(image,wcs=wcs,cmap='gray',**implot_kwargs)
    if show_aperture:
      if self.configs[config_name].keys() >= {'slit_width','slit_length','PA'}:
        slit = SkyRectangularAperture(self.configs[config_name]['coordinates'],
                                    w=self.configs[config_name]['slit_width'],
                                    h=self.configs[config_name]['slit_length'],
                                    theta=self.configs[config_name]['PA'])
        slit.to_pixel(wcs).plot(color='r',lw=3)
      elif self.configs[config_name].keys() >= {'fiber_radius'}:
        fiber = SkyCircularAperture(self.configs[config_name]['coordinates'],
                                    r=self.configs[config_name]['fiber_radius'])
        fiber.to_pixel(wcs).plot(color='r',lw=3)
      else:
        raise KeyError('''show_slit set to true, but this configuration does not have 'slit_width','slit_length', and 'PA' set, which are needed for slit display, or 'fiber_radius' set, for circular aperture.''')
    return fig, ax          
  
  def list_configurations(self):
    df = pd.DataFrame.from_dict(self.configs,orient='index')
    df['coordinates'] = [i.to_string() for i in df['coordinates'] if isinstance(i,SkyCoord)]
    if 'offset star' in df.columns:
      df['offset star'] = [i.to_string() if isinstance(i,SkyCoord) else np.nan for i in df['offset star']]
    if 'user_images' in df.columns:
      df['user_images'] = ['Y' if isinstance(i,dict) else np.nan for i in df['user_images']]
    df.index.name = 'configurations'
    df = df.replace({np.nan: '---'})
    return df
  
  def nudge_configuration(self,config_name,arcsec_east,arcsec_north):
    '''
    Nudge the coordinates of a configuration east or north in arcsec
    for better alignment.

    Parameters
    ----------
    config_name: str
      name of configuration to nudge
    arcsec_east: float
      amount to nudge east (west is negative) in arcsec
    arcsec_north: float
      amount to nudge north (south is negative) in arcsec
    '''
    new_coordinate = self.configs[config_name]['coordinates'].directional_offset_by(0,arcsec_north*u.arcsec)
    new_coordinate = new_coordinate.directional_offset_by(90,arcsec_east*u.arcsec)
    self.configs[config_name]['coordinates'] = new_coordinate


  @property
  def configurations(self):
    return self.list_configurations()