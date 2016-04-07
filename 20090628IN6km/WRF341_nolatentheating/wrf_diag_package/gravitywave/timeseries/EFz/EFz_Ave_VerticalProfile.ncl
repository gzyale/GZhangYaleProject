; a script to demonstrate deplaning, high-pass filtering, 
; computing vertical energy flux, and low-passed vertical energy flux

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./funcs.ncl"

begin

  setvalues NhlGetWorkspaceObjectId()
    "wsMaximumSize": 1000000000
  end setvalues

nlon = 329
nlat = 313


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; start of looping

fList = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01*")
fList_basename = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01* | xargs -n1 basename")
nFiles = dimsizes(fList)

levs = ispan(2000,20000,1000)*1.0;
nlevs = dimsizes(levs)

efzAve = new((/nFiles,nlevs/),double)

do iFile = 0, nFiles-1
  filename = stringtocharacter(fList_basename(iFile))
  a = addfile(fList(iFile),"r")
 
 ; first read 3-D variables
  p = wrf_user_getvar(a,"p",0)  ; Pa
  w = wrf_user_getvar(a,"wa",0) ; m/s
  z = wrf_user_getvar(a,"z",0)  ; m

  do ilev = 0, nlevs-1
    ; interpolate to constant height
    p_plane = wrf_user_intrp3d(p,z,"h",levs(ilev),0,False)
    w_plane = wrf_user_intrp3d(w,z,"h",levs(ilev),0,False)
    p_plane@description = "Original Pressure"
    w_plane@description = "Original w"

    ; have 2-D fields at a constant height
    ; now, deplane these fields
    pdp = deplane(p_plane)
    wdp = deplane(w_plane)
    pdp@description = "Deplaned Pressure (Pa)"
    wdp@description = "Deplaned w (m/s)"

    ; now, high-pass filter these fields to get the perturbation parts
    L = 400.0 ; cut-off length scale in km
    dx = (a@DX)*0.001   ; grid resolution in km
    pp = hpfilter(pdp,dx,L)
    wp = hpfilter(wdp,dx,L)
    pp@description = "Perturbation Pressure (Pa)"
    wp@description = "Perturbation w (m/s)"
    
    tmp = pp*wp
    efzAve(iFile,ilev) = avg(tmp(111:130,174:210)) ; tmp(lat,lon)
  end do
 
   delete(a)
   delete(filename)
end do 
  
;;;; Output
 outputFileName = "EFz_Ave_VerticalProfile.nc"

 system("rm -f " + outputFileName)   ; remove any pre-existing file
 ncout = addfile(outputFileName,"c")  ; open output netCDF file
 filedimdef(ncout,"time",-1,True)
 efzAve!0    = "time"                     ; assign named dimensions
 efzAve!1    = "lev"
 
 efzAve@long_name = "Zonal Energy Flux"   ; assign attributes 
 efzAve@units     = "W/m2"
 ncout->efzAve = efzAve
;;;; End of Output


end

