
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
  def __init__(self,
            target_list: list,
            observatory: str,
            obsdates: list,
            utcoffset: float):
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
    self.observatory = observatory
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
    self.total_configs = 0
    for i in target_list:
      self.total_configs += len(i.configurations)

  
  def plot_visibility(self,
                    date: str,
                    target: str = 'all',
                    view_range: float = 12,
                    plot_current_time: bool = False,
                    figsize: tuple = (30,12),
                    alt_min: float = 10,
                    alt_max: float = 90,
                    legend=True):
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
      delta_midnight = np.linspace(-view_range, view_range, 400)*u.hour
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
          for j in i.configs.keys():
            if not isinstance(i.configs[j]['coordinates'],SkyCoord):
              continue
            else:
              coord = i.configs[j]['coordinates']
              break
          ax.plot(obs_times.plot_date,coord.transform_to(frame).alt,label=name,lw=4)

      moon = get_moon(obs_times).transform_to(frame)
      sun = get_sun(obs_times).transform_to(frame)

      ax.plot(obs_times.plot_date,moon.alt,color='lightgray',ls='--',lw=3,label='Moon')
      mask1 = sun.alt.to(u.deg) < -0*u.deg
      ax.fill_between(obs_times.plot_date[mask1],0,90,color='0.9',zorder=0)
      mask2 = sun.alt.to(u.deg) < -12*u.deg
      ax.fill_between(obs_times.plot_date[mask2],0,90,color='0.7',zorder=0)
      mask3 = sun.alt.to(u.deg) < -18*u.deg
      ax.fill_between(obs_times.plot_date[mask3],0,90,color='0.2',zorder=0)

      # ax.fill_between(obs_times.plot_date, 0*u.deg, 90*u.deg,
      #              sun.alt.to(u.deg) > -0*u.deg, color='0.9', zorder=0)
      # ax.fill_between(obs_times.plot_date, 0*u.deg, 90*u.deg,
      #                  sun.alt.to(u.deg) > -12*u.deg, color='0.7', zorder=0)
      # ax.fill_between(obs_times.plot_date, 0*u.deg, 90*u.deg,
      #                  sun.alt.to(u.deg) > -18*u.deg, color='0.2', zorder=0)




      if plot_current_time:
        now = Time.now()
        ax.axvline(now.plot_date,color='r',lw=3)

      ax.set_ylim(alt_min,alt_max)
      if legend:
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
  
  def export_targetlist(self,
                      include_extras: list =[],
                      include_offset_stars: bool = True,
                      save_path: str = './',
                      name: str = 'targetlist'):
    '''
    Export an observatory-compliant targetlist of all targets and configurations.
    If only one configuration exists, the name column will be the target name.
    Otherwise, it will be targetname_configname. It's worth noting that keeping
    both of these short is advantageous for many facsum ingesting tools.
    By default, only key info (name,ra,dec,equinox) are added. Extra info
    from the configs can be added and will be appended as commented lines 
    below each row.

    Parameters
    ----------
    include_extras: list, default: []
      extra keywords to include in comments. Code will try to add them if they
      exist for a given configuration. (e.g., PAs or offsets)
    include_offset_stars: bool, default: True
      whether to include offset stars as entries in the targetlist. Name format is <target>_<config>_os.
    save_path: str, default: './'
      path to save targetlist to. default is current directory.
    name: str, default: 'targetlist'
      name of the file. 


    Returns
    -------
    None
      but saves the relevant targetlist file.
    
    Notes
    -----
    CURRENT_SUPPORTED_OBS: Palomar Observatory
    '''
    if not save_path.endswith('/'):
      save_path = save_path + '/' + name +'.csv'
    else:
      save_path = save_path + name +'.csv'
    if isinstance(self.observatory,str):
      if self.observatory.upper() == 'PALOMAR':
        write_str = ''
        for target in self.target_list:
          df = target.configurations
          if len(df) == 1:
            config = list(target.configs.keys())[0]
            radec = target.configs[config]['coordinates'].to_string(style='hmsdms',sep=':')
            #ra = target.configs[config]['coordinates'].ra.value
            #dec = target.configs[config]['coordinates'].dec.value
            ra = radec.split(' ')[0].replace(':',' ')
            dec = radec.split(' ')[1].replace(':',' ')
            write_str += f'{target.name},{ra},{dec},J2000 \n'
            mini_str = '! ^^ '
            for j in include_extras:
              try:
                mini_str += f'{j}: {target.configs[config][j]},'
              except:
                continue 
            mini_str += '\n'
            if len(include_extras)>0:
              write_str += mini_str 
            if 'offset star' in list(target.configs[config].keys()):
              if include_offset_stars:
                write_name = target.name + '_os'
                radec = target.configs[config]['offset star'].to_string(style='hmsdms',sep=':')
                ra = radec.split(' ')[0].replace(':',' ')
                dec = radec.split(' ')[1].replace(':',' ')
                #ra = target.configs[config]['offset star'].ra.value
                #dec = target.configs[config]['offset star'].dec.value
                write_str += f'{write_name},{ra},{dec},J2000 \n'
          elif len(df)>1:
            for config in list(target.configs.keys()):
              write_name = target.name + '_' + config 
              radec = target.configs[config]['coordinates'].to_string(style='hmsdms',sep=':')
              #ra = target.configs[config]['coordinates'].ra.value
              #dec = target.configs[config]['coordinates'].dec.value
              ra = radec.split(' ')[0].replace(':',' ')
              dec = radec.split(' ')[1].replace(':',' ')
              write_str += f'{write_name},{ra},{dec},J2000 \n'
              mini_str = '! ^^ '
              for j in include_extras:
                try:
                  mini_str += f'{j}: {target.configs[config][j]},'
                except:
                  continue
              mini_str += '\n'
              if len(include_extras)>0:
                write_str += mini_str
              if 'offset star' in list(target.configs[config].keys()):
                if include_offset_stars:
                  write_name = target.name + '_' + config + '_os'
                  radec = target.configs[config]['offset star'].to_string(style='hmsdms',sep=':')
                  #ra = target.configs[config]['offset star'].ra.value
                  #dec = target.configs[config]['offset star'].dec.value
                  ra = radec.split(' ')[0].replace(':',' ')
                  dec = radec.split(' ')[1].replace(':',' ')
                  write_str = write_str + f'{write_name},{ra},{dec},J2000 \n' 
          else:
            print('?????')
        with open(save_path,'w') as f:
          f.write(write_str)
            


  def html_summary(self,
                date: str,
                save_dir: str,
                overwrite: bool = True,
                view_range: float = 12):
    '''
    Produce a 'beautiful' html output with the observing plan. 
    
    Parameters
    ----------
    date: str
      date for which to construct the report. Must be a date present in the obsdates provided.
    save_dir: str
      location to save the observing plan. Spawns an image directory for relevant plots. 
    overwrite: bool, optional (default True)
      When set to true, constituant elements (like finder charts, airmass plots, etc) will be re-computed and saved to disk. 
    '''
    
    if not os.path.exists(os.path.join(save_dir,f'ObservingPlan_{date}')):
      os.mkdir(os.path.join(save_dir,f'ObservingPlan_{date}'))
    if not os.path.exists(os.path.join(save_dir,f'ObservingPlan_{date}','img')):
      os.mkdir(os.path.join(save_dir,f'ObservingPlan_{date}','img'))                  
    save_airmass = os.path.join(f'ObservingPlan_{date}','img',f'visibility_{date}.jpg')
    title = f'Observing Plan for UTC {date}'
    if overwrite:
      fig, ax = self.plot_visibility(date,view_range=view_range)
      fig.savefig(save_airmass)
      plt.close('all')
    for target in self.target_list:
      if overwrite:
        fig,ax = self.plot_visibility(date,target=target.name,view_range=view_range)
        fig.savefig(os.path.join(save_dir,f'ObservingPlan_{date}','img',f'visibility_{target.name}_{date}.jpg'))
        plt.close('all')
    top = f'''
    <html>
        <head>
            <title>{title}</title>
               <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
         <link rel="stylesheet" href="style.css">
         <meta name="viewport" content="width=device-width, initial-scale=1">
         </head>
        <body>
           <div class='myheader'>
           <h2 class="display-3">Observing Plan for UTC {date}</h2>
           <p>{self.observatory} Observatory</p>
           <h4>Displaying {len(self.target_list)} Targets and {self.total_configs} Unique Configurations</h4>
           </div>
           
             <div class='visibility'>
            <h4 class='display-4'>Visibility (all targets)</h4>
            </div>
       
          <div align='center'>
            <img src={os.path.join('img',f'visibility_{date}.jpg')} width="1300"> 
            </div>
           
            
            
    '''
    text0 = "<div class='container'>"
    for target in self.target_list:     
      # Let's make the coordinates row specific links.
      df = target.configurations
      df['coordinates'] = [f"<a href='https://www.legacysurvey.org/viewer?ra={i.split(' ')[0]}&dec={i.split(' ')[1]}&layer=ls-dr9&zoom=12' target='_blank'>{i}</a>" for i in df.coordinates]
      
      text = f'''
              <hr>
              <hr>
             <div align='center'>
             <h3 class='display-4'>Target: {target.name}</h3>
             
             <p></p>
             <a href="https://www.legacysurvey.org/viewer?ra={target.coordinates.to_string().split(' ')[0]}&dec={target.coordinates.to_string().split(' ')[1]}&layer=ls-dr9&zoom=12" target="_blank"><button type="button" class="btn btn-outline-primary btn-lg">Coordinates: {target.coordinates.to_string()}</button></a>
             <p></p>
             </div>
             <img src={os.path.join('img',f'visibility_{target.name}_{date}.jpg')} class="img-fluid" width="1300"> 
             <div class='box'>
             <button class='btn btn-outline-secondary btn-lg'><h1>Observing Configurations</h1></button>
             <p></p>
             <p class="lead">
  Configurations for {target.name} are shown below. Coordinate links lead to Legacy Survey Viewer centered at coordinates.
</p>
             </div>
             
             

             {df.to_html(escape=False,classes=["table table-striped",'table table-hover'])}
             <div class='box'>             
             <button class='btn btn-outline-success btn-lg'><h1>Finder Charts</h1></button>
             <p></p>
             <p class="lead">
  Click on the images below to open a new tab with a larger finder image. 
</p>
</div>
             <div class="row">
             '''
      # Determine size to use for cutouts. 
      if 'finder_size' in list(target.configs.keys()):
        size = target.configs['finder_size']
      elif 'slit_length' in list(target.configs.keys()):
        size = 1.5*target.configs['slit_length']
      else:
        size = 200*u.arcsec
      
      for key in target.configs.keys():
        if overwrite:
          if 'image_scaling' in list(target.configs.keys()):
            scale = target.configs['image_scaling']
            fig,ax=target.retrieve_finder_chart(key,size=size,scale=scale)
            plt.tight_layout()
            fig.savefig(os.path.join(save_dir,f'ObservingPlan_{date}','img',f'{target.name}_config_{key}_{date}.jpg'))
          else:
            fig,ax=target.retrieve_finder_chart(key,size=size,scale=1.0)
            plt.tight_layout()
            fig.savefig(os.path.join(save_dir,f'ObservingPlan_{date}','img',f'{target.name}_config_{key}_{date}.jpg'))
        img_path = os.path.join('img',f'{target.name}_config_{key}_{date}.jpg')
        ims = f'''
              <div class="col-lg-4 col-md-4 col-xs-4 thumb">
              <figure class="text-center">
              <a href="{img_path}" target="_blank"><img src="{img_path}" class="figure-img img-fluid rounded" alt="..." width='550'></a>
              <figcaption class="figure-caption">Finder Chart for Configuration: {key}</figcaption>
              </figure>
              </div> 
        '''
      
        text = text+ims 
      text = text+'</div></div>'
      top = top+text0+text
    close = f'''
    </div>
    </div>
    <div class='myfooter'>
    <p>Generated by observing-suite Software. Copyright Imad Pasha 2022</p>
    </div>
    </body>
    </html>
    '''
    final = top + close
    with open(os.path.join(save_dir,f'ObservingPlan_{date}','html_report.html'), 'w') as f:
      f.write(final)
    css_string = '''

body{
background: #ffffff;
}

.myheader {
  padding: 80px;
  text-align: center;
  background: #0a2648;
  color: white;
  font-size: 30px;
}

.box {
  padding: 60px;
  text-align: center;
}

.visibility {
  padding-top: 150 px; 
  text-align: center;
}


.myfooter {
  padding: 10px;
  text-align: center;
  background: #0a2648;
  color: white;
  font-size: 15px;
}
    '''
    with open(os.path.join(save_dir,f'ObservingPlan_{date}','style.css'), 'w') as f:
      f.write(css_string)
    