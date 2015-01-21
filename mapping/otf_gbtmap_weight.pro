;+
; Weight function used by OTF_GBTMap_Grid. It is based on Mangum et al. (2007) paper.
;
; <B> Contributed by: Jaime E. Pineda</B>
;
; <p> It assumes that the distance is in units of beam/(3 pixels), unless 
;  uniform weight is selected and the distance is in units of beam/(2 pixel)
;
; <p> Convolving with a Sinc function multiplied by a Gaussian (Mangum et al. 2007)
; The parameters used to create maps with cellsize = FWHM/3 are a=1.55 b=2.52
;
; <p> Convolving with a Bessel_1 function multiplied by a Gaussian (Mangum et al. 2007)
; The parameters used to create maps with cellsize = FWHM/3 are b=2.52 c=1.55
;
; @param dx {in}{required}{type=double} Array containing distance in x-axis.
; @param dy {in}{required}{type=double} Array containing distance in y-axis.
; @param weight {in}{required}{type=string} Type of weight used. Currently Sinc, 
;   Bessel and Uniform weights are implemented.
;
; @version $Id: otf_gbtmap_weight.pro,v 1.0 2009/01/28 16:11:03 jpineda Exp $
;-
FUNCTION OTF_GBTMap_Weight, dx, dy, weight
compile_opt IDL2
Case strupcase(Weight) of
  'BESSEL' : begin
    b=2.52 & c=1.55 & Rsup=3.
    b2 = 1.d/(b*b)
    distance = sqrt( dx^2+dy^2)
    wt = distance*0.0
    ind = where( distance LE Rsup, cc)
    if cc GT 0 then begin
      dmin=0.0001d
      distance = distance>dmin ;; this removes both singularities at d=0.0 !!
      wt[ind]= exp(-distance[ind]^2*b2)*BESELJ(distance[ind]/c, 1, /double)/(distance[ind]/c)
    endif
  end
  'SINC' : begin
    a  =1.55 & b=2.52 & Rsup=3.
    pia=!dpi/a
    b2 = 1.d/(b*b)
    distance = sqrt( dx^2+dy^2)
    wt = distance*0.0
    ind = where( distance LE Rsup, cc)
    if cc GT 0 then begin
      dmin=0.0001d
      distance = distance>dmin ;; this removes both singularities at d=0.0 !!
      wt[ind]= exp(-distance[ind]^2*b2)*Sin(pia*distance[ind])/(pia*distance[ind])
    endif
  end
  'UNIFORM' : begin
    distance = max([[abs(dx)],[abs(dy)]], dimension=2) 
    wt       = distance*0.0
    ind      = where( distance LE 1.0, cc)
    if cc GT 0 then wt[ind] = 1.0
  end
  ELSE : begin
    message, 'Weight function not defined', /con
    return, -1
  end
ENDCASE
return, wt
END
