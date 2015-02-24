;+
; Regrids GBT data in a file or directory into a FITS cube.
; 
; <B> Contributed by: Jaime E. Pineda</B>
;
; <p> It removes baseline and align spectra before regridding. Proper tcals 
; and opacity corrections should have been applied previously.
;
; <p> It uses several functions from IDLAstro library to make sure the header
; is AOK.
; The default projection used is TAN, but it can be easily changed.
;
; <p> Map_DBeam can reduce the raw GBT data (it applies Tcal and tau corrections) 
; and produce files ready for OTF_GBTMap_Grid to process.
;
; @param file_in {in}{required}{type=string} Input file with data reduced
; @param file_out {in}{required}{type=string} Output file
; @keyword dir {in}{optional}{type=boolean} File_In is a directory
;
; @keyword source_name {in}{optional}{type=string} Source name that will 
;   be added to the FITS header. If none given, it uses the source in !g.s[0].source
; @keyword blwin {in}{optional}{type=long} An array setting the window(s) 
;   for baseline fitting. It is passed to nregion function. If no region 
;   is given, then the whole spectrum is fitted.
; @keyword blorder {in}{optional}{type=integer} Order of polynomio fitted 
;   in the baseline removal. It's is passed to the baseline function. (Default=1)
;
; @keyword save_file {in}{optional}{type=string} Output filename for IDL 
;   save function. Just in case you want to check some intermediate steps (not tested).
; @keyword cont {in}{optional}{type=boolean} It tries to use a previously 
;   saved file (not tested).
; @keyword weight {in}{optional}{type=string} Weight function to be used: 
;   'Besse', 'Sinc' or 'Uniform'. 'Uniform' should be used for raster masp.
;   Sinc is used by default. Parameter is passed to OTF_GBTMap_Weight.
;
; @keyword pol {in}{optional}{type=integer} Polarizations that will be stored. 
;   If not provided it uses pl=[0,1].
; @keyword fdnum {in}{optional}{type=integer} Beams that will be stored. 
;   If not provided it uses fdnum=[0,1].
;
; @keyword beam {in}{optional}{type=double} Beam size (in degrees). It is used 
;   in the convolving function, and also stored in the FITS header. 
;   Default is 31.32 arcsec.
; @keyword pix {in}{optional}{type=double} Pixel size (in degrees) of the 
;   datacube. Default is beam/3.
; @keyword naxis {in}{optional}{type=long} Size of the final datacube. By 
;   default the header is set up to store all data.
; @keyword crpix {in}{optional}{type=double} Array with reference pixels 
;   for RA-DEC projection (velocity has not been implemented).
;   By default it is set at the center of the map.
; @keyword crval {in}{optional}{type=double} Array with reference RA-DEC values
;   to be used in the header.
;
; @keyword nu {in}{optional}{type=double} Rest frequency of line observed. 
;   Needed to compute the velocity of the line.
;   If none is given then a search within a VERY small list of frequencies
;   is done. It is added to the FITS header.
; @keyword mol {in}{optional}{type=string} Molecule observed. It 
;   will be added to the FITS header.
; @keyword transition {in}{optional}{type=string} Transition observed. It
;   will be added to the FITS header.
;
; @keyword freq {in}{optional}{type=boolean} Output file has frequency as 
;   third axis. (This is the default behaviour)
; @keyword vel {in}{optional}{type=boolean} Output file has velocity as 
;   third axis.
; @keyword nwin {in}{optional}{type=Long} The range in channels kept in 
;   the final cube. By default it uses the whole IF.
;
; @examples
; OTF_GBTMap_Grid, 'TMC-1C/NH3-11', 'cubes/TMC-1C_NH3_11.fits', /dir
;
; OTF_GBTMap_Grid, 'reduced/mar23_09.0.fits', 'cubes/NH3_11.fits', /vel
;
;; This is an example that shows many of the features in the program
;; It creates a cube of NH3 (1,1) emission line (several windows have to 
;; be defined to avoid all components)
; beam=31.d/60./60.
; pix=0.5*beam            ; pixel size is beam/2, where beam=31.arcsec
; ra0=15*ten(03,47,36.0)  ; center of map
; dec0=ten(32,51,00.0)    ;
;
; crval=[ra0,dec0]
; naxis=[45,56]
; crpix=[22.5, 29.5]
; pol=[0,1]
; fdnum=[0,1]
; nwin=[2000, 3500]
;
; blorder=1
; blwin=[1500, 2160, 2266, 2463, 2572, 2650, 2785, 2857, 2977, 3170, 3285, 4100]
; filein ='test_red/B5-11'
; fileout='test_sou/B5_11_f.fits'
; nu =23694.4955d6 ; NH3 (1,1) Lovas (2004)
; mol='NH3'
; transition='(1,1)'
; OTF_GBTMap_Grid, filein, fileout, /dir, source_name='B5', pix=pix, beam=beam, $
;         crval=crval, naxis=naxis, crpix=crpix, pol=pol, weight='BESSEL', $
;         save_file='B5-11.save', /freq, nwin=nwin, fdnum=fdnum, $
;         blwin=blwin, blorder=blorder, nu=nu, mol=mol, transition=transition
;
; @version $Id: otf_gbtmap_grid.pro,v 1.1 2011/04/18 18:02:04 jpineda Exp $
;               Updated to use bult-in handling of multi-beam header information.
;               Data is now process as float (to reduce data volume used), this is 
;               also consistent with native data format in GBTIDL (i.e., we are 
;               degrading the data).
; @version $Id: otf_gbtmap_grid.pro,v 1.0 2009/01/28 16:10:11 jpineda Exp $
;-
PRO OTF_GBTMap_Grid, file_in, file_out, source_name=source_name, $
           save_file=save_file, $
           blwin=blwin, blorder=blorder, $
           cont = cont, dir = dir, weight=weight, $
           pix=pix, beam=beam, $
           pol=pol, fdnum=fdnum, crpix=crpix, naxis=naxis, crval=crval, $
           nu=nu, freq = freq, vel = vel, nwin=nwin, $
           mol=mol, transition=transition

Compile_opt idl2
On_Error,2
if N_Params() LT 2 then begin
  message, 'Not enough parameters', /con
  return
endif

if (not keyword_set(save_file) and keyword_set(cont)) then begin
  message, 'can not continue a previous run without save_file file', /con
  return
endif
t0=systime(/sec)

if not keyword_set(weight) then weight='SINC'

if otf_gbtmap_weight( [0], [0], strupcase(weight)) eq -1 then begin
  message, 'Error! Weight function is not recognized', /con
  return
endif
; pixel and beam size
if n_elements(beam) eq 0 then $
  beam=0.522/60.   ;; beam size is 31.32 arcsec
if n_elements(pix) EQ 0 then $
  pix=beam/3.      ;; 3 pixels per beam

if keyword_set(save_file) and keyword_set(cont) then $
  goto, previous

if N_elements(crval) eq 2 then begin
  ra0  = crval[0]
  dec0 = crval[1]
endif
;
freeze
if keyword_set(dir) then $
  dirin, file_in $
else $
  filein,  file_in
;
; List of indexes is created and then read
; Is there a way to do this without creating a file?
;
if n_elements(fdnum) EQ 0 then fdnum=[0,1]
if n_elements(pol)   EQ 0 then pol=[0,1]

emptystack
select, plnum=pol, fdnum=fdnum, /quiet
n=astack()
;;n=n[where(n GT 0)]
getrec,n[0]
aux =getdata()
if keyword_set(nwin) then begin
  nmin=min(nwin)
  nmax=max(nwin)
endif else begin
  nmin=0L
  nmax=n_elements( aux)-1L
endelse
if not keyword_set(nu) then begin
; 
; Set of line frequencies that we'll use.
; These frequencies are in MHz. More lines can be easily added.
;
  nu0=[23694.4955d6, $ ; NH3 (1,1) Lovas (2004)
       23722.6333d6, $ ; NH3 (2,2) Lovas (2004)
       23870.1300d6, $ ; NH3 (3,3) Lovas (2004)
       24137.4200d6, $ ; NH3 (4,4) Lovas (2004)
       22344.0330d6, $ ; CCS (2_1 - 1_0) Yamamoto et al. (1990)
       23963.9010d6  $ ; HC5N (9-8) Lovas (2004)   this transition has not been observed in the lab
        ]
  mol0=['NH3', $
        'NH3', $
        'NH3', $
        'NH3', $
        'CCS', $
        'HC5N' $
        ]
  transition0=['(1,1)', $
        '(2,2)', $
        '(3,3)', $
        '(4,4)', $
        '2_1 - 1_0', $
        '9 - 8' $
        ]
;
  aux=MIN( nu0-!g.s[0].center_frequency, index, /absolute) ; search for the closest one.
  nu=nu0[index]                                            ; and assign it as the rest freq.
  if not keyword_set(mol)   then mol=mol0[index]           ; also save the molecule 
  if not keyword_set(transition) then transition=transition0[index] ; and transition at that frequency
  
endif else begin ;; if frequency is given then mol and transition must be given
  if not keyword_set(mol)   then mol='UNDEF'
  if not keyword_set(trans) then transition='UNDEF'
endelse
!g.s[0].line_rest_frequency=nu                             ; Rest frequency is updated
xfreq=frequency_axis(!g.s[0],frame='LSR')                  ; Get the frequency axis
xfreq=xfreq[nmin:nmax]                                     ; and select the desired range
;
velm=velocity_axis(!g.s[0],frame='LSR',veldef='RADIO')     ; Get velocity axis
velm=velm[nmin:nmax]                                       ; and select the desired range
; Create data and rms matrices
sz   = N_ELEMENTS(n)
nvel = N_ELEMENTS(velm)
data = FLTARR(sz,nvel)
radec= FLTARR(sz,2)
rms2 = FLTARR(sz)*!values.f_NaN
; load data into accum to align spectra later on
sclear
getrec, n[0]
accum
if not keyword_set(blorder) then blorder=1
if not keyword_set(blwin) then blwin=[nmin, nmax]
; Loop over records: baseline removal, 

FOR i=0L, sz-1L DO BEGIN
  getrec, n[i]                                  ; get spectrum
  print, i, sz-1, string(13B), format='("i=",i," / ",i,a1,$)'
  if total( finite(getdata()) ) eq 0 then begin  ; if all data is blanked
    rms2[i] = !values.f_nan                      ; do this... is there a
    radec[i,*]=[!values.f_nan, !values.f_nan]    ; status keyword in getrec?
  endif else begin
    nregion, blwin
    baseline, nfit=blorder
    print,fshift()
    gshift, fshift()                              ; align spectrum
    radec[i,*]=GetRaDec(!g.s[0])                  ; coordinates are stored
    aux=GetData()
    rms2[i]  =(robust_sigma(aux[nmin:nmax])>0)^2  ; rms is calculated in the respective window
    data[i,*]=aux[nmin:nmax]                      ; spectrum is stored
    if rms2[i] eq 0 then print, n[i]
  endelse
ENDFOR
print,''
t1=systime(/sec)
; 
; Remove all invalid scans
gd=where(finite(rms2) and rms2 ne 0)
data =data[gd,*]
rms2 =rms2[gd]
radec=radec[gd,*]
; save file with data if desired
if keyword_set(save_file) then $
  save, file=save_file, data, radec, rms2, nvel, velm, xfreq, nu
previous:
if keyword_set(save_file) and keyword_set(cont) then $
  restore, file=save_file
;
radec0=[MAX(radec[*,0], /nan), MIN(radec[*,1], /nan)]
mycd  =[[-1,0],[0,1]]*pix  ;; square pixels
; number of extra pixels for the map 
epsilon=1
if N_Elements(naxis) ne 2 then begin
  dra =(radec[*,0]-radec0[0])*cos(!dtor*radec0[1])
  ddec=(radec[*,1]-radec0[1])
  nax1=CEIL(MAX( ABS(dra/mycd[0,0] ), /nan ))
  nax2=CEIL(MAX( ABS(ddec/mycd[1,1]), /nan ))
  nax1=nax1+2*epsilon
  nax2=nax2+2*epsilon
endif else begin
  nax1=naxis[0]
  nax2=naxis[1]
endelse

; prepare to create final image and cube
fdata=FLTARR(nax1,nax2,nvel)*!VALUES.F_NAN
mkhdr, hd, fdata
aux=!g.s.units
nv0=1.0
v0 =velm[0]
dv =velm[1]-velm[0]
xfreq0 =xfreq[0]
dxfreq =xfreq[1]-xfreq[0]
if n_elements(ra0) eq 0 then begin
  ra0 =median(radec[*,0])
  dec0=median(radec[*,1])
endif
radec0=[ra0, dec0]
if N_Elements(crpix) ne 2 then begin
;
; crpix=[nax1, nax2]*0.5-0.5
;
; crpix is set to the center of map
  crpix=round( ( [ra0,dec0]-[MAX(radec[*,0],/nan), MIN(radec[*,1],/nan)] )$
         *[cos(!dtor*dec0), 1 ]$
         /[mycd[0,0],mycd[1,1]]*2)*0.5+[1,1]+[epsilon,epsilon]
endif
; set the FITS header
equinox=!g.s[0].equinox
putast, hd, mycd, double(crpix), radec0, ['RA---TAN','DEC--TAN'], EQUINOX=equinox, cd_type=1
if aux[0] EQ 'Jy' then $
  sxaddpar, hd, 'BUNIT', 'Jy' $
else $
  sxaddpar, hd, 'BUNIT', 'K', aux[0]
if keyword_set(vel) then begin
  sxaddpar, hd, 'CTYPE3', 'VRAD', after='ctype2'
  sxaddpar, hd, 'CUNIT3',  'm/s', after='cunit2'
  sxaddpar, hd, 'CDELT3',     dv, after='cdelt2'
  sxaddpar, hd, 'CRPIX3',    nv0, after='crpix2'
  sxaddpar, hd, 'CRVAL3',     v0, after='crval2'
endif else begin
  sxaddpar, hd, 'WCSNAME', 'LSRK-Freq'
  sxaddpar, hd, 'CTYPE3', 'FREQ', after='ctype2'
  sxaddpar, hd, 'CUNIT3',   'Hz', after='cunit2'
  sxaddpar, hd, 'CDELT3', dxfreq, after='cdelt2'
  sxaddpar, hd, 'CRPIX3',    nv0, after='crpix2'
  sxaddpar, hd, 'CRVAL3', xfreq0, after='crval2'
  sxaddpar, hd, 'SPECSYS', 'LSRK'
endelse
sxaddpar, hd, 'NAXIS3',   nvel, after='naxis2'
sxaddpar, hd, 'RESTFRQ', nu, 'Hz'
sxaddpar, hd, 'MOLECULE', mol
sxaddpar, hd, 'TRANSITI', transition
sxaddpar, hd, 'RADESYS', 'FK5'
; 
; Not sure if this will really help
;
;sxaddpar, hd, 'RESTFRQR', nu, 'Hz'
;sxaddpar, hd, 'CTYPE3R', 'VRAD'
;sxaddpar, hd, 'CNAME3R', 'LSRK Radio Velocity'
;sxaddpar, hd, 'CUNIT3R',  'm/s'
;sxaddpar, hd, 'CDELT3R',     dv
;sxaddpar, hd, 'CRPIX3R',    nv0
;sxaddpar, hd, 'CRVAL3R',     v0
;sxaddpar, hd, 'SPECSYSR', 'LSRK'
;sxaddpar, hd, 'SSYSOBSR', 'LSRK'
;
; Beam is defined and its circular
sxaddpar, hd, 'BMAJ', beam
sxaddpar, hd, 'BMIN', beam
sxaddpar, hd, 'BPA',  0.0d
  
if keyword_set(source_name) then $  ; Store source name if provided
  sxaddpar, hd, 'OBJECT', source_name, after='NAXIS3' $
else $                              ; or it is obtained from the first scan
  sxaddpar, hd, 'OBJECT', !g.s[0].source, after='NAXIS3'
adxy,hd, radec[*,0], radec[*,1], ii, jj
xind=findgen(nax1)#replicate(1,nax2)
yind=replicate(1,nax1)#findgen(nax2)
;
; i index loops over final image pixel
; j index loops over spectra
;
b3pix=beam/3./pix  ; (FWHM/3) in pixel units (the same units as distance) for weight function
if weight eq 'UNIFORM' then b3pix=0.5
FOR i=0L, nax1*nax2-1L DO BEGIN
  wt  = OTF_GBTMap_Weight( (ii-xind[i])/b3pix, (jj-yind[i])/b3pix, weight)
  ;iwt = where( wt GT 0.0, cwt)
  iwt = where( wt GT 0.0 and rms2 ne 0., cwt)
  if cwt GT 0 then begin
    wt[iwt] /=rms2[iwt]
    fdata[xind[i],yind[i],*]=0.0
    FOR j=0L, cwt-1L DO BEGIN
      fdata[ xind[i], yind[i], *] += data[iwt[j],*]*wt[iwt[j]]
    ENDFOR
    fdata[ xind[i] , yind[i] ,*] /= total(wt[iwt])
  endif
ENDFOR
unfreeze
spawn, '\rm '+file_out
fileout, file_out
writefits, file_out, fdata, hd
t2=systime(/sec)
print, 'it took '+STRING( (t2-t0)/60., format='(F6.2)')+' minutes to create the map'
print, '        '+STRING( (t1-t0)/60., format='(F6.2)')+' minutes to load the data'
print, '        '+STRING( (t2-t1)/60., format='(F6.2)')+' minutes to regrid'
;stop
END
