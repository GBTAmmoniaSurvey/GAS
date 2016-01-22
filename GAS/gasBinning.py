import numpy as np
from GAS.gridregion import channelShift
import scipy.ndimage as nd
import astropy.units as u
from spectral_cube import SpectralCube
from GAS import voronoi_2d_binning as v2d

def BinByLabel(DataCube, CentroidMap, LabelMap,
               CentroidAggregator = np.nanmean, BackgroundLabels = [0]):
    """
    Bin a data cube by a label mask, aligning the data to a common centroid.

    Parameters
    ----------
    DataCube : SpectralCube
        The original spectral cube with spatial dimensions Nx, Ny and spectral dimension Nv
    CentroidMap : 2D numpy.ndarray
        A 2D map of the centroid velocities for the lines to stack of dimensions Nx, Ny.
        Note that DataCube and Centroid map must have the same spectral units (e.g., km/s)
    LabelMap : 2D numpy.ndarray
        A 2D map containing integer labels for each pixel into objects defining the stacking.
    CentroidAggregator : numpy.ufunc
        Operates on a vector of centroid data and returns the value summarizing that object.
    BackgroundLabels : list
        List of values in the label map that correspond to background objects and should not
        be processed with the stacking. 

    Returns
    -------
    OutputCube : SpectralCube
        A SpectralCube instance matching the input but with spectra aligned in velocity and averaged.

    """
    UniqLabels = np.unique(LabelMap)
    ChannelWidth = np.median(DataCube.spectral_axis-
                             np.roll(DataCube.spectral_axis,1))
    RawData = DataCube.unmasked_data[:].value
    for ThisLabel in UniqLabels:
        if ThisLabel not in BackgroundLabels:
            y,x= np.where(ThisLabel == LabelMap)
            CentroidValue = CentroidAggregator(CentroidMap[y,x])
            AccumSpec = np.zeros(DataCube.shape[0])
            for ThisX,ThisY in zip(x,y):
                DeltaV = CentroidMap[ThisY,ThisX] - CentroidValue
                # Note this assumes the units of the centroid map
                # are in same units as the spectral axis of the cube.
                DeltaChan = DeltaV/ChannelWidth.value
                AccumSpec += channelShift(DataCube[:,ThisY,ThisX].value,
                                          -DeltaChan)
            AccumSpec /= x.size
            AccumSpec.shape = AccumSpec.shape+ (1,)
            RawData[:,y,x] = AccumSpec    
    return (SpectralCube(data = RawData,wcs = DataCube.wcs))

def VoronoiBin(IntegratedIntensity, NoiseMap, TargetValue= 5,
               threshold = 0, mask = None, aggregator = v2d._sn_func):
    """
    This routine wraps the 2D Binning approach of Cappellari and Copin with altered
    functionality to make it more generally useful in arbitrary spectral line cases.
    The output of the labeling process should be used in `BinByLabel`.
    
    Parameters
    ----------
    IntegratedIntensity : 2D numpy.ndarray
        Map containing the data to be processed by the binning algorithm.
    NoiseMap : 2D numpy.ndarray
        Map containing the RMS at each point.  This must be the same dimensions
        as the IntegratedIntensity.
    TargetValue : number
        Value that the binning should try to achieve for the aggregator
        function (see below).  Since aggregation defaults to signal-to-noise,
        the value would be the target S/N.
    threshold : number
        Positions with IntegratedIntensity less than the threshold are not binned.
    mask : 2D numpy.ndarray
        Binary mask that tests False where aggregation should be ignored.
    aggregator : function
        Function with call signature f(Signal, Noise, index) which returns a single number
        for pixels of value index.  
    """
    
    SignalToNoise = IntegratedIntensity/NoiseMap
    if mask is None:
        mask = SignalToNoise > threshold
    x, y = np.where(mask)
    signal = IntegratedIntensity[x,y]
    noise = NoiseMap[x,y]
    labels = np.zeros(IntegratedIntensity.shape)
    SNmap = np.zeros(IntegratedIntensity.shape)
    voronoi_label, x0, y0, \
        xbar, ybar, sn, npix, scale = v2d.voronoi_2d_binning(
        x, y, signal, noise, TargetValue, plot = False, aggregator= aggregator)
    for ThisLabel,SNval in enumerate(sn):
        label_index = np.where(voronoi_label == ThisLabel)
        if SNval > TargetValue:
            labels[x[label_index],y[label_index]] = ThisLabel 
            SNmap[x[label_index],y[label_index]] = sn[ThisLabel]
    return(labels,SNmap)
