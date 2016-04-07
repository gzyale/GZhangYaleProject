; a script to demonstrate deplaning, high-pass filtering, 
; computing vertical energy flux, and low-passed vertical energy flux

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./funcs.ncl"

begin

  setvalues NhlGetWorkspaceObjectId()
    "wsMaximumSize": 1000000000
  end setvalues

  ; level of interest
;  lev = 16000.0 ; m

nlon = 329
nlat = 313
nlev = 21

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; start of looping

fList = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01*")
fList_basename = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01* | xargs -n1 basename")
nFiles = dimsizes(fList)

do iFile = 0, nFiles-1
  filename = stringtocharacter(fList_basename(iFile))
  a = addfile(fList(iFile),"r")
 
 ; first read 3-D variables
  p = wrf_user_getvar(a,"p",0)  ; Pa
  u = wrf_user_getvar(a,"ua",0) ; m/s
  z = wrf_user_getvar(a,"z",0)  ; m

  efx = new((/nFiles, nlev, nlat, nlon/), double);
  efx@description = "Zonal Energy Flux (W/m2)"

  do lev = 0, 1000, 20000
  ; interpolate to constant height
  p_plane = wrf_user_intrp3d(p,z,"h",lev,0,False)
  u_plane = wrf_user_intrp3d(u,z,"h",lev,0,False)
  p_plane@description = "Original Pressure"
  u_plane@description = "Original u"

  ; have 2-D fields at a constant height
  ; now, deplane these fields
  pdp = deplane(p_plane)
  udp = deplane(u_plane)
  pdp@description = "Deplaned Pressure (Pa)"
  udp@description = "Deplaned u (m/s)"

  ; now, high-pass filter these fields to get the perturbation parts
  L = 400.0 ; cut-off length scale in km
  dx = (a@DX)*0.001   ; grid resolution in km
  pp = hpfilter(pdp,dx,L)
  up = hpfilter(udp,dx,L)
  pp@description = "Perturbation Pressure (Pa)"
  up@description = "Perturbation u (m/s)"

  efx(iFile,lev,:,:) = pp*up

  end do
;;;; Output
 outputFileName = "EFx"+filename(7:23)+".nc"

 xlat = a->XLAT(0,:,:)
 xlong= a->XLONG(0,:,:) 

 system("rm -f " + outputFileName)   ; remove any pre-existing file
 ncout = addfile(outputFileName,"c")  ; open output netCDF file
 filedimdef(ncout,"time",-1,True)
 efx!0    = "lat"                     ; assign named dimensions
 efx!1    = "lon"
 ncout->XLAT = xlat
 ncout->XLONG = xlong 
 
 efx@long_name = "Zonal Energy Flux"   ; assign attributes 
 efx@units     = "W/m2"
 ncout->efx = efx
;;;; End of Output

  delete(a)
  delete(filename)
end do

end

