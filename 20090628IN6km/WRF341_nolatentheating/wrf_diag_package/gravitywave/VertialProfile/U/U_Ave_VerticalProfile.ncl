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

uAve = new((/nFiles,nlevs/),double)

do iFile = 0, nFiles-1
  filename = stringtocharacter(fList_basename(iFile))
  a = addfile(fList(iFile),"r")
 
 ; first read 3-D variables
  p = wrf_user_getvar(a,"p",0)  ; Pa
  u = wrf_user_getvar(a,"ua",0) ; m/s
  z = wrf_user_getvar(a,"z",0)  ; m

  do ilev = 0, nlevs-1
    ; interpolate to constant height
    u_plane = wrf_user_intrp3d(u,z,"h",levs(ilev),0,False)
    u_plane@description = "Original u"

    uAve(iFile,ilev) = sum(u_plane*maskpoly)/sum(maskpoly)
  end do
 
   delete(a)
   delete(filename)
end do 
  
;;;; Output
 outputFileName = "U_Ave_VerticalProfile.nc"

 system("rm -f " + outputFileName)   ; remove any pre-existing file
 ncout = addfile(outputFileName,"c")  ; open output netCDF file
 filedimdef(ncout,"time",-1,True)
 uAve!0    = "time"                     ; assign named dimensions
 uAve!1    = "lev"
 
 uAve@long_name = "Zonal Wind"   ; assign attributes 
 uAve@units     = "m/s"
 ncout->uAve = uAve
;;;; End of Output


end

