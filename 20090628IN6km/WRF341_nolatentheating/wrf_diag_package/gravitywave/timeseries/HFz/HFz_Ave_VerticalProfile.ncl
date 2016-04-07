; a script to demonstrate deplaning, high-pass filtering, 
; computing vertical energy flux, and low-passed vertical energy flux

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./funcs.ncl"

begin

  setvalues NhlGetWorkspaceObjectId()
    "wsMaximumSize": 1000000000
  end setvalues

lev = 5000.0 ; m

nlon = 329
nlat = 313


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; start of looping

fList = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01*")
fList_basename = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01* | xargs -n1 basename")
nFiles = dimsizes(fList)

levs = lev;
nlevs = dimsizes(levs)

hfzAve = new((/nFiles,nlevs/),double)

do iFile = 0, nFiles-1
  filename = stringtocharacter(fList_basename(iFile))
  a = addfile(fList(iFile),"r")
 
 ; first read 3-D variables
  tk = wrf_user_getvar(a,"tk",0)  ; K
  qv = wrf_user_getvar(a,"QVAPOR",0) ;kg/kg
  p = wrf_user_getvar(a,"p",0)  ; Pa
  u = wrf_user_getvar(a,"ua",0) ; m/s
  w = wrf_user_getvar(a,"wa",0) ; m/s
  z = wrf_user_getvar(a,"z",0)  ; m

  Tv = (1.+(0.622*qv))*tk
  rho = p/(287*Tv)


  do ilev = 0, nlevs-1
    ; interpolate to constant height
    w_plane = wrf_user_intrp3d(w,z,"h",levs(ilev),0,False)
    t_plane = wrf_user_intrp3d(tk,z,"h",levs(ilev),0,False)
    rho_plane = wrf_user_intrp3d(rho,z,"h",levs(ilev),0,False)
    w_plane@description = "Original w"
    t_plane@description = "Original T"
    rho_plane@description = "Original air density"

    ; have 2-D fields at a constant height
    ; now, deplane these fields
    wdp = deplane(w_plane)
    tdp = deplane(t_plane)
    wdp@description = "Deplaned w"
     tdp@description = "Deplaned T"
    
    ; now, high-pass filter these fields to get the perturbation parts
    L = 400.0 ; cut-off length scale in km
    dx = (a@DX)*0.001   ; grid resolution in km
    wp = hpfilter(wdp,dx,L)
    tp = hpfilter(tdp,dx,L)
    wp@description = "Perturbation w"
    tp@description = "Perturbation T"
    
    tmp = wp*tp
    hfzAve(iFile,ilev) = avg(tmp(111:130,174:210)) ; tmp(lat,lon)
  end do
 
   delete(a)
   delete(filename)
end do 
  
;;;; Output
 outputFileName = "HFz_Ave_VerticalProfile.nc"

 system("rm -f " + outputFileName)   ; remove any pre-existing file
 ncout = addfile(outputFileName,"c")  ; open output netCDF file
 filedimdef(ncout,"time",-1,True)
 hfzAve!0    = "time"                     ; assign named dimensions
 hfzAve!1    = "lev"
 
 hfzAve@long_name = "Zonal Heat Flux"   ; assign attributes 
 hfzAve@units     = "W/m2"
 ncout->hfzAve = hfzAve
;;;; End of Output


end

