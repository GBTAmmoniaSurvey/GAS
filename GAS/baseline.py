import numpy as np
import numpy.polynomial.legendre as legendre
from .utils import VlsrByCoord
from .first_look import trim_edge_cube
import catalogs
from scipy.optimize import least_squares as lsq
from spectral_cube import SpectralCube
import astropy.units as u
import astropy.utils.console as console
import pyspeckit.spectrum.models.ammonia_constants as acons


def ammoniaWindow(spectrum, spaxis, freqthrow=4.11 * u.MHz,
                  window=3, v0=8.5, line='oneone', outerwindow=None):
    """
    This defines a narrow window around the v0 value for ammonia hyperfines.

    Parameters
    ----------
    spectrum : np.array
        The one-dimensional spectrum
    spaxis : np.array
        The spectral axis in units of km/s for the input spectrum
    window : np.float
        Width, in km/s, around the input velocity to consider emission
    outerwindow : np.float
        Velocity separations larger than this value are ignored in the baseline
    v0 : np.float
        Central velocity in km/s for emission window
    freqthrow : astropy.Quantity
        frequency swith throw for the observations.
    line : str
        String name of the ammonia line to consider, e.g., 'oneone', 'twotwo'
    Returns
    -------
    mask : np.array
        Boolean array of regions to use in the baseline fitting.
    """

    mask = np.zeros_like(spectrum, dtype=np.bool)
    voffs = np.array([dv for dv in acons.voff_lines_dict[line]])
    for voff in voffs:
        mask[(spaxis > (v0 + voff - window)) *
             (spaxis < (v0 + voff + window))] = True

    deltachan = freqthrow / ((spaxis[1] - spaxis[0]) / 299792.458 *
                             acons.freq_dict[line] * u.Hz)
    deltachan = deltachan.to(u.dimensionless_unscaled).value
    deltachan = (np.floor(np.abs(deltachan))).astype(np.int)

    if (deltachan < spaxis.size):
        mask = np.logical_or(mask, np.r_[mask[deltachan:-1],
                                         np.zeros(deltachan + 1,
                                                  dtype=np.bool)])
        mask = np.logical_or(mask, np.r_[np.zeros(deltachan + 1,
                                                  dtype=np.bool),
                                         mask[0:(-deltachan - 1)]])

    if outerwindow is not None:
        mask[(spaxis > (v0 + outerwindow + voffs.max()))] = True
        mask[(spaxis < (v0 - outerwindow + voffs.min()))] = True
    return ~mask


def tightWindow(spectrum, spaxis,
                window=4,
                outerwindow=None,
                v0=8.5, freqthrow=4.11 * u.MHz):
    """
    This defines a narrow window around the v0 value for the emission region

    Parameters
    ----------
    spectrum : np.array
        The one-dimensional spectrum
    spaxis : np.array
        The spectral axis in units of km/s for the input spectrum
    window : np.float
        Width, in km/s, around the input velocity to consider emission
    outerwindow : np.float
        Velocity separations larger than this value are ignored in the baseline
    v0 : np.float
        Central velocity in km/s for emission window
    freqthrow : astropy.Quantity
        frequency swith throw for the observations.
    Returns
    -------
    mask : np.array
        Boolean array of regions to use in the baseline fitting.
    """

    mask = np.zeros_like(spectrum, dtype=np.bool)
    mask[(spaxis > (v0 - window)) * (spaxis < (v0 + window))] = True
    deltachan = freqthrow / ((spaxis[1] - spaxis[0]) /
                             299792.458 * 23.5 * u.GHz)
    deltachan = deltachan.to(u.dimensionless_unscaled).value
    deltachan = (np.floor(np.abs(deltachan))).astype(np.int)

    if (deltachan < spaxis.size):

        mask = np.logical_or(mask, np.r_[mask[deltachan:-1],
                                         np.zeros(deltachan + 1,
                                                  dtype=np.bool)])
        mask = np.logical_or(mask, np.r_[np.zeros(deltachan + 1,
                                                  dtype=np.bool),
                                         mask[0:(-deltachan - 1)]])
    if outerwindow is not None:
        mask[(spaxis > (v0 + outerwindow))] = True
        mask[(spaxis < (v0 - outerwindow))] = True
    return(~mask)


def mad1d(x):
    med0 = np.median(x)
    return np.median(np.abs(x - med0)) * 1.4826


def legendreLoss(coeffs, y, x, noise):
    return (y - legendre.legval(x, coeffs)) / noise


def robustBaseline(y, baselineIndex, blorder=1, noiserms=None):
    x = np.linspace(-1, 1, len(y))
    if noiserms is None:
        noiserms = mad1d((y - np.roll(y, -2))[baselineIndex])
    opts = lsq(legendreLoss, np.zeros(blorder + 1), args=(y[baselineIndex],
                                                          x[baselineIndex],
                                                          noiserms),
               loss='arctan')

    return y - legendre.legval(x, opts.x)


def rebaseline(filename, blorder=3, 
               baselineRegion=[slice(0, 262, 1), slice(-512, 0, 1)],
               windowFunction=None, blankBaseline=False,
               flagSpike=True, v0=None, trimEdge=True, **kwargs):
    """
    Rebaseline a data cube using robust regression of Legendre polynomials.

    Parameters
    ----------
    filename : string
        FITS filename of the data cube
    blorder : int
        Order of the polynomial to fit to the data
    baselineRegion : list
        List of slices defining the default region of the spectrum, in
        channels, to be used for the baseline fitting.
    windowFunction : function
        Name of function to be used that will accept spectrum data, and
        velocity axis and will return a binary mask of the channels to be used
        in the  baseline fitting.  Extra **kwargs are passed to windowFunction
        to do with as it must.
    blankBaseline : boolean
        Blank the baseline region on a per-spectrum basis
    trimEdge : boolean
        remove the map edges

    Returns
    -------
    Nothing.  A new FITS file is written out with the suffix 'rebaseN' where N
    is the baseline order

    """
    cube = SpectralCube.read(filename)
    originalUnit = cube.spectral_axis.unit
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')
    spaxis = cube.spectral_axis.to(u.km / u.s).value

    goodposition = np.isfinite(cube.apply_numpy_function(np.max, axis=0))
    y, x = np.where(goodposition)
    outcube = np.zeros(cube.shape) * np.nan
    RegionName = (filename.split('_'))[0]

    if hasattr(windowFunction, '__call__'):
        catalog = catalogs.GenerateRegions()

    nuindex = np.arange(cube.shape[0])
    runmin = nuindex[-1]
    runmax = nuindex[0]

    for thisy, thisx in console.ProgressBar(zip(y, x)):
        spectrum = cube[:, thisy, thisx].value

        if v0 is not None:
            baselineIndex = windowFunction(spectrum, spaxis,
                                           v0=v0, **kwargs)
        elif hasattr(windowFunction, '__call__'):
            _, Dec, RA = cube.world[0, thisy, thisx]
            # This determines a v0 appropriate for the region
            v0 = VlsrByCoord(RA.value, Dec.value, RegionName,
                             regionCatalog=catalog)
            baselineIndex = windowFunction(spectrum, spaxis,
                                           v0=v0, **kwargs)
        else:
            baselineIndex = np.zeros_like(spectrum,dtype=np.bool)
            for ss in baselineRegion:
                baselineIndex[ss] = True

        runmin = np.min([nuindex[baselineIndex].min(), runmin])
        runmax = np.max([nuindex[baselineIndex].max(), runmax])

        # Use channel-to-channel difference as the noise value.
        if flagSpike:
            jumps = (spectrum - np.roll(spectrum, -1))
            noise = mad1d(jumps) * 2**(-0.5)
            baselineIndex *= (np.abs(jumps) < 5 * noise)
            noise = mad1d((spectrum -
                           np.roll(spectrum, -2))[baselineIndex]) * 2**(-0.5)
        else:
            noise = mad1d((spectrum -
                           np.roll(spectrum, -2))[baselineIndex]) * 2**(-0.5)

        if blankBaseline:
            spectrum = robustBaseline(spectrum, baselineIndex,
                                      blorder=blorder,
                                      noiserms=noise)
            spectrum[baselineIndex] = np.nan
            outcube[:, thisy, thisx] = spectrum
        else:
            outcube[:, thisy, thisx] = robustBaseline(spectrum, baselineIndex,
                                                      blorder=blorder,
                                                      noiserms=noise)

    outsc = SpectralCube(outcube, cube.wcs, header=cube.header)
    outsc = outsc[runmin:runmax, :, :]  # cut beyond baseline edges
    if trimEdge: 
        trim_edge_cube(cube)
    # Return to original spectral unit
    outsc = outsc.with_spectral_unit(originalUnit)
    outsc.write(filename.replace('.fits', '_rebase{0}.fits'.format(blorder)),
                overwrite=True)
