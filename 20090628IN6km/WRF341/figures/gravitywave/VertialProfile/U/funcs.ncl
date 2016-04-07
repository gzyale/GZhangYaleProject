
; functions contained in this file
; "calc_fourierposx"
; "calc_fourierposy"
; "lpfilter"
; "hpfilter"
; "make1d"
; "deplane"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

undef("calc_fourierposx")
function calc_fourierposx(var:numeric,dx:numeric)

begin

  ; setup the wavenumber magnitude matrix.
  coef     = fft2df(var)    ; 2d fft to get dims of fft
  coefdims = dimsizes(coef) ; dimensions of fft
  nl = coefdims(1)          ; number coeffs in y direction
  nk = coefdims(2)          ; number coeffs in x direction
  NY = coef@M              ; original number of points in y dir
  NX = coef@N              ; original number of points in x dir

  ; span of wave numbers in x and y directions (from 0 to Nyquist)
  dy = dx 
  k = ispan(0,nk-1,1)/todouble(NX)/dx
  l = ispan(0,nl-1,1)/todouble(NY)/dy
  inds = ind(l.gt.1.0/2.0/dy)
  l(inds) = l(inds)-2.0/2.0/dy
;  if NY%2.eq.0 then ; if even
;    l1 = ispan(0,nl/2-1,1)/dy/todouble((NY))
;    l2 = ispan(-nl/2,-1,1)/dy/todouble((NY))
;    l  = array_append_record(l1,l2,0)
;  end if 
;  if NY%2.eq.1 then ; if odd
;    l1  = ispan(0,(nl-1)/2,1)/dy/todouble(NY-1)
;    l2  = ispan(-(nl-1)/2,-1,1)/dy/todouble(NY-1)
;    l   = array_append_record(l1,l2,0)
;  end if 

  kk = new((/nl,nk/),double) ; 2-D east-west wavenumber array
  ;ll = new((/nm,nn/),double)
  ;do ki = 0,nk-1
  ;  ll(:,ki) = l
  ;end do
  do li = 0,nl-1
    kk(li,:) = k
  end do

  return(kk)
end

undef("calc_fourierposy")
function calc_fourierposy(var:numeric,dx:numeric)

begin

  ; setup the wavenumber magnitude matrix. will be valid at
  ; every level and every time since the horizontal grid does
  ; not change with height and time. Do this once outside of
  ; the time/height loop to increase efficiency
  coef     = fft2df(var)    ; 2d fft to get dims of fft
  coefdims = dimsizes(coef) ; dimensions of fft
  nl = coefdims(1)          ; number coeffs in y direction
  nk = coefdims(2)          ; number coeffs in x direction
  NY = coef@M              ; original number of points in y dir
  NX = coef@N              ; original number of points in x dir

  ; span of wave numbers in x and y directions (from 0 to Nyquist)
  dy = dx 
  k = ispan(0,nk-1,1)/todouble(NX)/dx
  l = ispan(0,nl-1,1)/todouble(NY)/dx
  inds = ind(l.gt.1.0/2.0/dy)
  l(inds) = l(inds) - 2.0/2.0/dy
;  if NY%2.eq.0 then ; if even
;    l1 = ispan(0,nl/2-1,1)/dy/NY
;    l2 = ispan(-nl/2,-1,1)/dy/NY
;    l  = array_append_record(l1,l2,0)
;  end if 
;  if NY%2.eq.1 then ; if odd
;    l1  = ispan(0,(nl-1)/2,1)/dy/NY
;    l2  = ispan(-(nl-1)/2,-1,1)/dy/NY
;    l   = array_append_record(l1,l2,0)
;  end if 

  ;kk = new((/nl,nk/),double) ; 2-D east-west wavenumber array
  ll = new((/nl,nk/),double)
  do ki = 0,nk-1
    ll(:,ki) = l
  end do
  ;do li = 0,nl-1
  ;  kk(li,:) = k
  ;end do

  return(ll)
end

undef("lpfilter")
function lpfilter(varu:numeric,dx:numeric,L:numeric)

begin

  var = todouble(varu)

  dims  = dimsizes(var)
  varlp = new(dims,double)

  kk = calc_fourierposx(var,todouble(dx))
  ll = calc_fourierposy(var,todouble(dx))
  kappa2 = kk*kk + ll*ll

  ; Now, do the filtering
  fc = 1./L ; spatial "frequency"

  ; With wave number magnitude matrix (kappa), can calculate 
  ; the low pass 2-D response function in fourier
  ; space. Not a function of time or height! Should greatly
  ; increase efficiency. Can subtract the lp field from the full field
  ; to get the hp field

  ; calculate the 2-D response function in fourier space
  rlp = exp(-1.0*kappa2/fc^2)
  dims = dimsizes(rlp)

  ; stack this response function array for simple point-wise multiplication of
  ; fourier coefficient array. Takes care of real and imag coeffs.
  dimsrlp = dimsizes(rlp)
  rlp3d = new((/2,dimsrlp(0),dimsrlp(1)/),double)
  rlp3d(0,:,:) = rlp
  rlp3d(1,:,:) = rlp

  ; perform the filtering. This single line includes three sets: 1. 2-d FFT. 
  ; 2. Multiply 2-D FFT with the response function. 3. Perform inverse of appropriately
  ; reduced Fourier coefficients.
  coef  = fft2df(var)  ; 2d DFT
  coef  = coef*rlp3d   ; apply response function
  varlp = fft2db(coef) ; perform 2-D DFT

  ;print("Finished LP filtering")
  return(varlp)
end

undef("hpfilter")
function hpfilter(var:numeric,dx:numeric,L:numeric)

begin

  dims  = dimsizes(var)
  varhp = new(dims,double)

  kk = calc_fourierposx(var,todouble(dx))
  ll = calc_fourierposy(var,todouble(dx))
  kappa2 = kk*kk + ll*ll

  ; Now, do the filtering
  fc = 1./L ; spatial "frequency"

  ; With wave number magnitude matrix (kappa), can calculate 
  ; the low pass 2-D response function in fourier
  ; space. Not a function of time or height! Should greatly
  ; increase efficiency. Can subtract the lp field from the full field
  ; to get the hp field

  ; calculate the 2-D response function in fourier space
  rhp = 1.0 - exp(-1.0*kappa2/fc^2)
  dims = dimsizes(rhp)

  ; stack this response function array for simple point-wise multiplication of
  ; fourier coefficient array. Takes care of real and imag coeffs.
  dimsrhp = dimsizes(rhp)
  rhp3d = new((/2,dimsrhp(0),dimsrhp(1)/),double)
  rhp3d(0,:,:) = rhp
  rhp3d(1,:,:) = rhp

  ; perform the filtering. This single line includes three sets: 1. 2-d FFT. 
  ; 2. Multiply 2-D FFT with the response function. 3. Perform inverse of appropriately
  ; reduced Fourier coefficients.
  coef  = fft2df(todouble(var))  ; 2d DFT
  coef  = coef*rhp3d   ; apply response function
  varhp = fft2db(coef) ; perform 2-D DFT

  ;print("Finished LP filtering")
  return(varhp)
end

undef("make1d")
function make1d(var:numeric)

begin

  ; number of points in the 2-D array
  dims = dimsizes(var)
  ny   = dims(0)
  nx   = dims(1)
  n    = nx*ny

  ; define the 1d array
  oned = new(n,double)

  ; loop through and fill up 1-D array with slices of the 2-D array
  do j=0,ny-1
    oned(j*nx:((j+1)*nx-1)) = var(j,:)
  end do

  return(oned)
end

undef("deplane")
function deplane(var:numeric)

begin

  ; make a coordinate system and position arrays
  dims = dimsizes(var)
  ny   = dims(0)
  nx   = dims(1)

  skip = 2 ; number of points to skip eastward and northward. Choice of one
           ; means every other point is retain.

  indx = ispan(0,nx-1,skip+1)
  indy = ispan(0,ny-1,skip+1)

  nxs = dimsizes(indx)
  nys = dimsizes(indy)

  x  = new((/ny,nx/),double)
  y  = new((/ny,nx/),double)
  xs = new((/nys,nxs/),double)
  ys = new((/nys,nxs/),double)
  vars = new((/nys,nxs/),double)

  ; fill out the full x arrays
  do j=0,ny-1
    x(j,:) = ispan(0,nx-1,1)
  end do
  do i=0,nx-1
    y(:,i) = ispan(0,ny-1,1)
  end do

  ; fill out the subsampled x arrays
  do j=0,nys-1
    xs(j,:) = indx
  end do
  do i=0,nxs-1
    ys(:,i) = indy
  end do

  ; fill out the subsampled var array
  do j=0,nys-1
    ;print(nys + " " + ny + " " + j + " " + j*(skip+1))
    vars(j,:) = var(j*(skip+1),indx)
  end do

  ; Now have positions and values. Need to make these arrays 1-D
  var1d = make1d(vars)
  x1d = make1d(xs)
  y1d = make1d(ys)

  n = dimsizes(var1d)
  ones = new(n,double)
  ones(:) = 1.0

  A = transpose((/ones,x1d,y1d/))
  beta = inverse_matrix(transpose(A)#A)#transpose(A)#var1d

  planar = beta(0) + beta(1)*x + beta(2)*y

  return(var-planar)
end
