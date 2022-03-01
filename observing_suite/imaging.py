import numpy as np 
import matplotlib.pyplot as plt 
from astropy.wcs import WCS 
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib.widgets import Slider, Button

def implot(image,figsize=(15,13),cmap='gray_r',scale=0.5,colorbar=False,header=None,wcs=None,**kwargs):
    '''
    Plot an astronomical image, setting default options and easy tweaking of parameters
    
    Parameters
    ----------
    image: array_like
        2D array containing an astronomical image to be plotted. Cutouts should be input as cutout.data.
    figsize: tuple, optional
        figure size to use. Default: (15,13)
    cmap: str, optional
        Colormap to use for the image. Default: 'gray_r'
    scale: float, optional
        By default, function will scale image to some number of standard deviations about the mean pixel value. Scale sets this number (or fraction). Default: 0.5.
    colorbar: bool, optional
        Whether to add a colorbar or not. Default: False
    header: dict, optional
        If input, function will attempt to create a WCS object from header and plot in celestial coordinates. Default: None
    wcs: WCS object
        If input, the function will plot using a projection set by the WCS. Default: None
    **kwargs
        Additional arguments are passed to matplotlib plotting commands. Currently supported: vmin, vmax.
        
    Returns
    -------
    fig, ax
        figure and axes objects containing currently plotted data.
    '''
    if (header==None) and (wcs==None):
        fig, ax = plt.subplots(figsize=figsize)
    elif wcs is not None:
        fig, ax = plt.subplots(figsize=figsize,subplot_kw={'projection':wcs})
        ax.set_xlabel('Right Ascension [hms]',fontsize=15)
        ax.set_ylabel('Declination [degrees]',fontsize=15)
        ax.coords.grid(color='gray', alpha=0.5, linestyle='solid')
    elif header is not None:
        wcs = WCS(header)
        fig, ax = plt.subplots(figsize=figsize,subplot_kw={'projection':wcs})
        ax.set_xlabel('Right Ascension [hms]',fontsize=15)
        ax.set_ylabel('Declination [degrees]',fontsize=15)
        ax.coords.grid(color='gray', alpha=0.5, linestyle='solid')
    mu = np.mean(image)
    s = np.std(image)
    dvmin = mu - scale*s
    dvmax = mu + scale*s
    if all(['vmin','vmax']) in kwargs.keys():
        im = ax.imshow(image,origin='lower',cmap=cmap,vmin=kwargs['vmin'],vmax=kwargs['vmax'])
    elif 'vmin' in kwargs.keys():
        im = ax.imshow(image,origin='lower',cmap=cmap,vmin=kwargs['vmin'],vmax=dvmax)
    elif 'vmax' in kwargs.keys():
        im = ax.imshow(image,origin='lower',cmap=cmap,vmin=dvmin,vmax=kwargs['vmax'])
    else:
        im = ax.imshow(image,origin='lower',cmap=cmap,vmin=dvmin,vmax=dvmax)
    if colorbar:
        cbar = plt.colorbar(im,ax=ax)
    ax.tick_params(direction='in',length=9,width=1.5,labelsize=15)
    
    return fig, ax

class ScalePlot():

    def __init__(self,image):
        self.image=image

    def plot(self,meanScalar=10,stdScalar=10,**kwargs):
        '''
        Function for plotting astronomical images, adding the convenience of scaling sliders and a colormap inverter.

        Parameters
        ----------
        image: array_like
            the image to be plotted. 
        kwargs: dict, optional
            any specific keywords for the `myutils.plotting.implot()` function 
        
        Returns
        -------
        fig, ax: the figure and axes objects
        '''
        fig, ax = implot(self.image,figsize=(10,10),**kwargs)
        aximage = ax.get_images()[0]
        divider = make_axes_locatable(ax)
        pad_fraction = 0.5
        width = axes_size.AxesY(ax, aspect=1./25)
        width2 = axes_size.AxesY(ax, aspect=1./10.)
        pad = axes_size.Fraction(pad_fraction, width)
        mean_ax = divider.append_axes("right", size=width, pad=pad)
        scale_ax = divider.append_axes("right",size=width,pad=pad+width)
        reverse_ax = divider.append_axes('top',size=width2,pad=pad)

        mean = Slider(mean_ax, 
                    'mean', 
                    (1./meanScalar)*np.mean(self.image), 
                    meanScalar*np.mean(self.image), 
                    valinit=np.mean(self.image), 
                    valstep=1,
                    orientation='vertical')
        scale = Slider(scale_ax, 
                    'scale',
                    (1./stdScalar)*np.mean(self.image), 
                    meanScalar*np.std(self.image), 
                    valinit=np.std(self.image),
                    orientation='vertical')
        
        reverse = Button(reverse_ax,'reverse colormap')
        def reverse_cmap(val):
            cmap = aximage.get_cmap().name
            if cmap.endswith('_r'):
                aximage.set_cmap(cmap.split('_')[0])
                fig.canvas.draw_idle()
            else:
                aximage.set_cmap(cmap+'_r')
                fig.canvas.draw_idle()

        def update(val):
            vmin = mean.val - scale.val 
            vmax = mean.val + scale.val
            aximage.set_clim(vmin,vmax)
            fig.canvas.draw_idle()

        reverse.on_clicked(reverse_cmap)
        mean.on_changed(update)
        scale.on_changed(update)

        plt.show()
        return reverse