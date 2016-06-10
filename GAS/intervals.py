import numpy as np
import astropy.units as u
import astropy.constants as con 

class VelocitySet(list):
    def __init__(self,input):
        if input is None:
            self.ranges = []
        else:
            self.ranges = []
            for elt in input:
                self.ranges.append(elt)

    def append(self,item):
        self.ranges.append(item)

    def __iter__(self):
        return iter(self.ranges)

    def __getitem__(self,index):
        return self.ranges[index]

    def __len__(self):
        return(len(self.ranges))

    def __repr__(self):
        return('Velocity set consisting of: \n\t'+''.join([elt.__repr__()+'\n\t' for elt in self.ranges]))

    def __add__(self,other):
        result = VelocitySet(None)
        if hasattr(other,'__iter__'):
            for otherelt in other.ranges:
                for selfelt in self.ranges:
                    result.ranges += selfelt+otherelt
        else:
            for selfelt in self.ranges:
                result.ranges += selfelt+other
        for idx, interval in enumerate(result.ranges):
            if interval.lower == interval.upper:
                result.ranges.pop(idx)
            
        return(result)

    def __sub__(self,other):
        result = VelocitySet(None)
        if hasattr(other,'__iter__'):
            for otherelt in other.ranges:
                for selfelt in self.ranges:
                    result.ranges += selfelt-otherelt
        else:
            for selfelt in self.ranges:
                result.ranges += selfelt-other
        for idx, interval in enumerate(result.ranges):
            if interval.lower == interval.upper:
                result.ranges.pop(idx)
        return(result)
    
    def toslice(self,cdelt = None, crpix = None, crval = None, vframe = 0.0*u.m/u.s, restfreq = None):
        slicelist = []
        for elt in self.ranges:
            ch1 =((restfreq*(1-((elt.lower+vframe)/con.c).to(u.dimensionless_unscaled))-crval)/cdelt)+crpix-1
            ch2 = ((restfreq*(1-((elt.upper+vframe)/con.c).to(u.dimensionless_unscaled))-crval)/cdelt)+crpix-1
            start = round(np.min([ch1,ch2]))
            stop = round(np.max([ch1,ch2]))
            slicelist += [slice(start,stop,1)]
        return(slicelist)

class VelocityInterval:
    def __init__(self,lower,upper):
        self.lower = lower
        self.upper = upper

    def __repr__(self):
        return('Velocity Interval from {0} to {1}'.format(self.lower,self.upper))
    
    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        boundaries = [self.lower,other.lower,self.upper,other.upper]
        casearray = (np.argsort(boundaries)) % 2
        if np.all(casearray == [0,0,1,1]):
            result = VelocitySet([self,other])
        if np.all(casearray == [1,1,0,0]):
            result = VelocitySet([other,self])
        if np.all(casearray == [0,1,1,0]):
            result = VelocitySet([self])
        if np.all(casearray == [1,0,0,1]):
            result = VelocitySet([other])
        if np.all(casearray == [0,1,0,1]):
            result = VelocitySet([VelocityInterval(self.lower,other.upper)])
        if np.all(casearray == [1,0,1,0]):
            result = VelocitySet([VelocityInterval(other.lower,self.upper)])
        return(result)

    def __sub__(self, other):
        boundaries = [self.lower,other.lower,self.upper,other.upper]
        casearray = (np.argsort(boundaries)) % 2
        if np.all(casearray == [0,0,1,1]):
            warnings.warn('Attempt to subtract non-overlapping intervals')
            result = VelocitySet(None)
        if np.all(casearray == [1,1,0,0]):
            warnings.warn('Attempt to subtract non-overlapping intervals')
            result = VelocitySet(None)
        if np.all(casearray == [0,1,1,0]):
            result = VelocitySet([VelocityInterval(self.lower,other.lower),VelocityInterval(other.upper,self.upper)])
        if np.all(casearray == [1,0,0,1]):
            result = VelocitySet([VelocityInterval(other.lower,self.lower),VelocityInterval(self.upper,other.upper)])
        if np.all(casearray == [0,1,0,1]):
            result = VelocitySet([VelocityInterval(self.lower,other.upper)])
        if np.all(casearray == [1,0,1,0]):
            result = VelocitySet([VelocityInterval(other.lower,self.upper)])
        return(result)
    
    def applyshift(self,offsetv):
        self.lower += offsetv
        self.upper += offsetv
