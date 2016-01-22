import numpy as np
from GAS.gridregion import channelShift
import scipy.ndimage as nd
import astropy.units as u
from spectral_cube import SpectralCube

def BinByLabel(DataCubeIn, CentroidMap, LabelMap,
               CentroidAggregator = np.nanmean, BackgroundLabels = [0]):
    """
    Bin a data cube by a label mask, aligning the data to a common centroid.
    DataCube is a SpectralCube
    """
    UniqLabels = np.unique(LabelMap)
    DataCube = DataCubeIn.with_spectral_unit(u.km/u.s,velocity_convention = 'radio')
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
                AccumSpec += channelShift(DataCube[:,ThisY,ThisX].value,-DeltaChan)
            AccumSpec /= x.size
            AccumSpec.shape = AccumSpec.shape+ (1,)
            RawData[:,y,x] = AccumSpec    
    return (SpectralCube(data = RawData,wcs = DataCube.wcs))
