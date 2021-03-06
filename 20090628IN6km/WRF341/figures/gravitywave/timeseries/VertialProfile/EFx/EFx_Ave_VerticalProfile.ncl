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

efxAve = new((/nFiles,nlevs/),double)

do iFile = 0, nFiles-1
  filename = stringtocharacter(fList_basename(iFile))
  a = addfile(fList(iFile),"r")
 
 ; first read 3-D variables
  p = wrf_user_getvar(a,"p",0)  ; Pa
  u = wrf_user_getvar(a,"ua",0) ; m/s
  z = wrf_user_getvar(a,"z",0)  ; m

  do ilev = 0, nlevs-1
    ; interpolate to constant height
    p_plane = wrf_user_intrp3d(p,z,"h",levs(ilev),0,False)
    u_plane = wrf_user_intrp3d(u,z,"h",levs(ilev),0,False)
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

    efxAve(iFile,ilev) = sum(pp*up*maskpoly)/sum(maskpoly)
    ;print(pp*up*maskpoly)
  end do
 
   delete(a)
   delete(filename)
end do 
  
;;;; Output
 outputFileName = "EFx_Ave_VerticalProfile.nc"

 system("rm -f " + outputFileName)   ; remove any pre-existing file
 ncout = addfile(outputFileName,"c")  ; open output netCDF file
 filedimdef(ncout,"time",-1,True)
 efxAve!0    = "time"                     ; assign named dimensions
 efxAve!1    = "lev"
 
 efxAve@long_name = "Zonal Energy Flux"   ; assign attributes 
 efxAve@units     = "W/m2"
 ncout->efxAve = efxAve
;;;; End of Output


end

