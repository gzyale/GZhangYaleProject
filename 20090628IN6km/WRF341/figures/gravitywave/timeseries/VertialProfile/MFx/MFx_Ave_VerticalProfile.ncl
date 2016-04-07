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

;;;generate a mask;;;;;;;;;;;;;;;;;;;;;;
aa = addfile("../../../../WRF/test/em_real/wrfout_d01_2009-06-28_00:00:00","r")
xlong = wrf_user_getvar(aa,"XLONG",0)
xlat = wrf_user_getvar(aa,"XLAT",0)
latpoly = (/8,8,20,20/)
lonpoly = (/77,79,76,74/)
boolmask = gc_inout(xlat,xlong,latpoly,lonpoly) ; returns array of True/False
maskpoly = new(dimsizes(boolmask),integer)
maskpoly = where(boolmask,1,0) ; where True = 1, False = 0  ; Gives mask array of 0s and 1s
delete(aa)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; start of looping

fList = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01*")
fList_basename = systemfunc("ls -1 ../../../../WRF/test/em_real/wrfout_d01* | xargs -n1 basename")
nFiles = dimsizes(fList)

levs = ispan(2000,20000,1000)*1.0;
nlevs = dimsizes(levs)

mfxAve = new((/nFiles,nlevs/),double)

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
    u_plane = wrf_user_intrp3d(u,z,"h",levs(ilev),0,False)
    rho_plane = wrf_user_intrp3d(rho,z,"h",levs(ilev),0,False)
    w_plane@description = "Original w"
    u_plane@description = "Original u"
    rho_plane@description = "Original air density"
    
    ; have 2-D fields at a constant height
    ; now, deplane these fields
    wdp = deplane(w_plane)
    udp = deplane(u_plane)
    wdp@description = "Deplaned w"
    udp@description = "Deplaned u (m/s)"

    ; now, high-pass filter these fields to get the perturbation parts
    L = 400.0 ; cut-off length scale in km
    dx = (a@DX)*0.001   ; grid resolution in km
    wp = hpfilter(wdp,dx,L)
    up = hpfilter(udp,dx,L)
    wp@description = "Perturbation w"
    up@description = "Perturbation u (m/s)"
    
    mfxAve(iFile,ilev) = sum(wp*up*maskpoly)/sum(maskpoly)
    ;print(pp*up*maskpoly)
  end do
 
   delete(a)
   delete(filename)
end do 
  
;;;; Output
 outputFileName = "MFx_Ave_VerticalProfile.nc"

 system("rm -f " + outputFileName)   ; remove any pre-existing file
 ncout = addfile(outputFileName,"c")  ; open output netCDF file
 filedimdef(ncout,"time",-1,True)
 mfxAve!0    = "time"                     ; assign named dimensions
 mfxAve!1    = "lev"
 
 mfxAve@long_name = "Zonal Energy Flux"   ; assign attributes 
 mfxAve@units     = "W/m2"
 ncout->mfxAve = mfxAve
;;;; End of Output


end

