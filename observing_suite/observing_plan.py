
from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_moon
from matplotlib import dates
from matplotlib.dates import HourLocator,MinuteLocator
import pandas as pd
from photutils.aperture import SkyRectangularAperture, SkyCircularAperture
import os

__all__ = ['ObservingPlan']

class ObservingPlan():
  def __init__(self,target_list,observatory,obsdates,utcoffset):
    '''
    Initialize ObservingPlan object. 
    
    Parameters
    ----------
    target_list: list
      list of Target() objects that are part of this observing plan.
    observatory: str or EarthLocation
      observatory for run. Either a string with the name of an observatory in the Astropy list, or an astropy EarthLocation object.
    obsdates: list or str
      list of dates as strings being observed in UTC, in form ['YYYY-MM-DD','YYYY-MM-DD'], or for a single date, 'YYYY-MM-DD'
    utcoffset: int
      the offset from utc in which the observations are taking place (e.g., US/Pacific is -8). 
    '''
    self.target_list = target_list
    if isinstance(observatory,str):
      try:
        self.obsloc = EarthLocation.of_site(observatory)
      except:
        raise AssertionError('Observatory not in database; add EarthLocation manually.')
    elif not isinstance(observatory,EarthLocation):
        raise AssertionError('If a non string is provided, it must be an astropy EarthLocation')
    
    self.obs_info = {}
    if isinstance(obsdates,str):
      # Parse Date 
      midnight = Time(f'{obsdates} 00:00:00') - utcoffset*u.hour
      self.obs_info[obsdates] = {}
      self.obs_info[obsdates]['midnight'] = midnight
    elif isinstance(obsdates,list):
      for date in obsdates:
        self.obs_info[date] = {}
        midnight = Time(f'{date} 00:00:00') - utcoffset*u.hour
        self.obs_info[date]['midnight'] = midnight 
    else:
      raise AssertionError('obsdates must be str or list')
   

  
  def plot_visibility(self,date,target='all',view_range=12,plot_current_time=False,figsize=(30,12),alt_min=10,alt_max=90):
      '''
      Produce a plot of altitude and airmass for targets on a given night of observing. 

      Parameters
      ----------
      date: str
        date on which to make the plot. In form 'YYYY-MM-DD'. Must be a date specified within obsdate.
      target: str or list, optional (default 'all')
        target to plot track for. deault is "all". Individual targets can be specified by name (str), via a list or, for single objects, a str.
      view_range: int, optional, default 12
        number of hours on either side of midnight over which to show the plot
      plot_current_time: bool, optional (default False)
        show a vertical bar at the current time (when code executed). Useful when making plot interactively during a night.
      figsize: tuple, optional (default (15,12))
        tuple specifying the figure size
      '''
      
      midnight = self.obs_info[date]['midnight']
      delta_midnight = np.linspace(-view_range, view_range, 1000)*u.hour
      obs_times = midnight+delta_midnight
      frame = AltAz(obstime=obs_times,location=self.obsloc)
      fig, ax = plt.subplots(figsize=figsize)
      # Parse Target Dictionary 
      for i in self.target_list:
        name = i.name
        if target != 'all':
          if isinstance(target,str):
            if name != target:
              continue
          elif isinstance(target,list):
            if name not in target:
              continue
        
        if i.coordinates is not None:
          ax.plot(obs_times.plot_date,i.coordinates.transform_to(frame).alt,label=name,lw=4)
        else:
          for j in configs.keys():
            if not isinstance(configs[j]['coordinates'],SkyCoord):
              continue
            else:
              coord = i.configs[j]['coordinates']
              break
          ax.plot(obs_times.plot_date,coord.transform_to(frame).alt,label=name,lw=4)

      moon = get_moon(obs_times).transform_to(frame)
      sun = get_sun(obs_times).transform_to(frame)
      ax.plot(obs_times.plot_date,moon.alt,color='lightgray',ls='--',lw=3,label='Moon')

      ax.fill_between(obs_times.plot_date, 0*u.deg, 90*u.deg,
                   sun.alt < -0*u.deg, color='0.9', zorder=0)
      ax.fill_between(obs_times.plot_date, 0*u.deg, 90*u.deg,
                       sun.alt < -12*u.deg, color='0.7', zorder=0)
      ax.fill_between(obs_times.plot_date, 0*u.deg, 90*u.deg,
                       sun.alt < -18*u.deg, color='0.2', zorder=0)




      if plot_current_time:
        now = Time.now()
        ax.axvline(now.plot_date,color='r',lw=3)

      ax.set_ylim(alt_min,alt_max)
      ax.legend(prop={'size': figsize[0]})
      ax.set_xlim([obs_times[0].plot_date, obs_times[-1].plot_date])
      loc = HourLocator(interval=1)
      date_formatter = dates.DateFormatter("%H")
      ax.xaxis.set_major_locator(loc)
      ax.xaxis.set_minor_locator(MinuteLocator(interval=30))
      ax.xaxis.set_major_formatter(date_formatter)
      #plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
      airmass_ticks = np.array([1,1.065,1.155,1.55, 2, 2.9])
      altitude_ticks = 90 - np.degrees(np.arccos(1/airmass_ticks))

      ax2 = ax.twinx()
      # ax2.invert_yaxis()
      ax2.set_yticks(altitude_ticks)
      ax2.set_yticklabels([round(i,1) for i in airmass_ticks])
      ax2.set_ylim(ax.get_ylim())
      ax2.set_ylabel('Airmass',fontsize=figsize[0])
      ax.set_ylabel('Altitude [deg]',fontsize=figsize[0])
      ax.set_xlabel('Time [UTC]', fontsize=figsize[0])
      ax.grid()
      #ax2.grid()
      ax.tick_params(labelsize=figsize[0])
      ax2.tick_params(labelsize=figsize[0])
      ax.tick_params(which='minor',direction='in',length=6,width=2,color='gray',top=True)
      return fig, ax 
  
  def html_summary(self,date,save_dir):
    '''
    Produce a 'beautiful' html output with the observing plan. 
    '''
    if not os.path.exists(os.path.join(save_dir,f'ObservingPlan_{date}')):
      os.mkdir(os.path.join(save_dir,f'ObservingPlan_{date}'))
    if not os.path.exists(os.path.join(save_dir,f'ObservingPlan_{date}','img')):
      os.mkdir(os.path.join(save_dir,f'ObservingPlan_{date}','img'))                  
    save_airmass = os.path.join('img',f'visibility_{date}.png')
    title = f'Observing Plan for UTC {date}'
    if not os.path.exists(os.path.join(save_dir,f'ObservingPlan_{date}','img',f'visibility_{date}.png')):
      fig, ax = self.plot_visibility(date)
      fig.savefig(save_airmass)
    
    top = f'''
    <html>
        <head>
            <title>{title}</title>
           
        </head>
        <body>
           <div align='center'>
            <h1>{title}</h1>
            <p></p>
            
            <h2>Visibility (all targets)</h2>
            <img src={str(save_airmass)} width="1500"> 
            </div>
    '''
    for target in self.target_list:
      text = f'''
             <p></p>
             <div align='center'>
             <hr>
             <hr>
             <h1>Target: {target.name}</h1>
             <h4>Coordinates: {target.coordinates.to_string()}</h4>
             <h3>Configurations</h3>
             {target.configurations.to_html()}
             </div>
             <div align='center'>
             '''
      for key in target.configs.keys():
        if not os.path.exists(os.path.join(save_dir,f'ObservingPlan_{date}','img',f'{target.name}_config_{key}_{date}.png')):
          fig,ax=target.retrieve_finder_chart(key,size=10*u.arcmin)
          fig.savefig(os.path.join(save_dir,f'ObservingPlan_{date}','img',f'{target.name}_config_{key}_{date}.png'))
        ims = f'''
              <img src={os.path.join('img',f'{target.name}_config_{key}_{date}.png')} width='500'>
        '''
        text = text+ims
      top = top + text
    close = f'''
    </div>
    </body>
    </html>
    '''
    with open(os.path.join(save_dir,f'ObservingPlan_{date}','html_report.html'), 'w') as f:
      f.write(top)