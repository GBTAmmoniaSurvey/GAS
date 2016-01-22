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
    DataCube is a SpectralCube
    """
    UniqLabels = np.unique(LabelMap)
    ChannelWidth = np.median(DataCube.spectral_axis-np.roll(DataCube.spectral_axis,1))
    RawData = DataCube.unmasked_data[:].value
    for ThisLabel in UniqLabels:
        if ThisLabel not in BackgroundLabels:
            y,x= np.where(ThisLabel == LabelMap)
            CentroidValue = CentroidAggregator(CentroidMap[y,x])
            AccumSpec = np.zeros(DataCube.shape[0])
            for ThisX,ThisY in zip(x,y):
                DeltaV = CentroidMap[ThisY,ThisX] - CentroidValue
                # Note this assumes the units of the centroid map are km/s
                DeltaChan = DeltaV/ChannelWidth.value
                AccumSpec += channelShift(DataCube[:,ThisY,ThisX].value,
                                          -DeltaChan)
            AccumSpec /= x.size
            AccumSpec.shape = AccumSpec.shape+ (1,)
            RawData[:,y,x] = AccumSpec    
    return (SpectralCube(data = RawData,wcs = DataCube.wcs))

def SignalBin(IntegratedIntensity, NoiseMap, targetSN = 5,
              threshold = 0, mask = None):
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
        x, y, signal, noise, targetSN, plot = False)
    for ThisLabel,SNval in enumerate(sn):
        label_index = np.where(voronoi_label == ThisLabel)
        if SNval > targetSN:
            labels[x[label_index],y[label_index]] = ThisLabel 
            SNmap[x[label_index],y[label_index]] = sn[ThisLabel]
    return(labels,SNmap)
