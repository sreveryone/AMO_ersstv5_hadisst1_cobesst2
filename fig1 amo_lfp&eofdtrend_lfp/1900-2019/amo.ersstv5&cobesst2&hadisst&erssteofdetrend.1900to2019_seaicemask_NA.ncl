load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;---------------------------------------------------------------------------
;                         User input
;---------------------------------------------------------------------------
gg                 = 9.806
pi                 = 3.14159265359
;pi                = 4.0*ATAN(1.0)
d2rad              = 0.01745329
;R0                = 6371393.0
;RADIUS            = 6371000.0
R0                 = 6371000.0
rad                = (pi/180.)
cp                 = 3996.0
p0                 = 1026.0
;-----------------------oissst----------------------
file_path1 = "/data05/zhaobw/atlantic/netcdf/SST/oisstv2/oisstv2.mon.198201to201912.subclim1982to2017.r360x180.nc"
f1 = addfile(file_path1, "r")
;-------- mask data-----------------
file_path3 = "/data05/zhaobw/atlantic/netcdf/SST/oisstv2/oisstv2.landseamask.r360x180.nc"
f3= addfile(file_path3, "r")
landseamask = f3->mask

;begin
time = f1->time
a = cd_calendar(time, -1)
ssta_init = f1->sst(ind(a.ge.198201.and.a.le.201912),:,:)
oisst = short2flt(ssta_init)
; ssta := dtrend_leftdim(ssta, 0)
lat = f1->lat
lon = f1->lon
oisst!0 = "time"
oisst!1 = "lat"
oisst!2 = "lon"
oisst&time = time
oisst&lat = lat
oisst&lon = lon
; sst anomaly -> sst anomaly mask(only ocean)
oisst_mask = oisst
ntime_oisst = dimsizes(oisst&time)
; mask land area
do t = 0, dimsizes(time)-1
    ; write loop content
    oisst_mask(t,:,:) = where(landseamask(0,:,:).eq.0, oisst@_FillValue, oisst(t,:,:))
end do
delete(oisst)
oisst = oisst_mask
; printVarSummary(oisst)
delete([/lat, lon, a, time, f1, oisst_mask/])
;---------------------- ERSSTv5-------------------------------------
    fil_ssta      = "/data05/zhaobw/atlantic/netcdf/SST/ersstv5/ersstv5.monmean.190001to201912.subclim1971.r360x180.nc"
    f1   = addfile(fil_ssta, "r")
    TIME   = f1->time
    lat = f1->lat
    lon = f1->lon
    nlat = dimsizes(lat)
    nlon = dimsizes(lon)
  yrStrt  = 1900          ; manually specify for convenience
  yrLast  = 2019          ; 20th century ends 2010
  TIME   = f1->time
  YYYYMM = cd_calendar(TIME, -1)
  ; print(YYYYMM)
  tStrt  = ind(YYYYMM .eq. (yrStrt*100 +  1))       ; indices 20th century
  tLast  = ind(YYYYMM .eq. (yrLast*100 + 12))       ; relative to entire dataset
    ersst = f1->sst(tStrt:tLast,:,:) ; (time, lat, lon)
    ; printVarSummary(ersst)
    ntime_ersst = dimsizes(ersst&time)
    delete([/fil_ssta, f1,TIME/])

;---------------------- ersstv5 after EOF1 detrend -------------------------------------
    fil_ssta      = "/data05/zhaobw/atlantic/netcdf/eof_ersstv5_lfca/ersstv5_detrendeof1_ssta.mon.1900to2019.nc"
    f1   = addfile(fil_ssta, "r")

    ersst_eof1dtrend = f1->ssta;(tStrt:tLast,:,:) ; (time, lat, lon)
    copy_VarCoords(ersst, ersst_eof1dtrend)
    copy_VarAtts(ersst, ersst_eof1dtrend)
    printVarSummary(ersst_eof1dtrend)
    ntime_ersst_eof1dtrend = dimsizes(ersst_eof1dtrend&time)
    ; printVarSummary(cobesst)
    delete([/fil_ssta, f1/])

;---------------------- COBESST -------------------------------------
    fil_ssta      = "/data05/zhaobw/atlantic/netcdf/SST/cobesst/cobesst.190001to201912.monmean.subclim1971.r360x180.nc"
    f1   = addfile(fil_ssta, "r")

  yrStrt  = 1900          ; manually specify for convenience
  yrLast  = 2019          ; 20th century ends 2010
  TIME   = f1->time
  YYYYMM = cd_calendar(TIME, -1)
  ; print(YYYYMM)
  tStrt  = ind(YYYYMM .eq. (yrStrt*100 +  1))       ; indices 20th century
  tLast  = ind(YYYYMM .eq. (yrLast*100 + 12))       ; relative to entire dataset

    cobesst = f1->sst(tStrt:tLast,:,:) ; (time, lat, lon)
    ntime_cobesst = dimsizes(cobesst&time)
    ; printVarSummary(cobesst)
    delete([/fil_ssta, f1, TIME/])
;---------------------- COBESST2 -------------------------------------
    fil_ssta      = "/data05/zhaobw/atlantic/netcdf/SST/cobesst/cobesst2.190001to201912.monmean.subclim1971.r360x180.nc"
    f1   = addfile(fil_ssta, "r")

  yrStrt  = 1900          ; manually specify for convenience
  yrLast  = 2019          ; 20th century ends 2010
  TIME   = f1->time
  YYYYMM = cd_calendar(TIME, -1)
  ; print(YYYYMM)
  tStrt  = ind(YYYYMM .eq. (yrStrt*100 +  1))       ; indices 20th century
  tLast  = ind(YYYYMM .eq. (yrLast*100 + 12))       ; relative to entire dataset

    cobesst2 = f1->sst(tStrt:tLast,:,:) ; (time, lat, lon)
    ntime_cobesst2 = dimsizes(cobesst2&time)
    ; printVarSummary(cobesst)
    delete([/fil_ssta, f1, TIME/])
; ;---------------------- HadISST -------------------------------------
;     fil_ssta      = "/data05/zhaobw/atlantic/netcdf/SST/hadisst/HadISST.mon.190001to201912.subclim1971.r360x180.nc"
;     f1   = addfile(fil_ssta, "r")
;     fil_ice       = "/data05/zhaobw/atlantic/netcdf/SST/hadisst/HadISST.ice.mon.190001to201912.r360x180.nc"
;     f2   = addfile(fil_ice, "r")

;   yrStrt  = 1900          ; manually specify for convenience
;   yrLast  = 2019          ; 20th century ends 2010
;   TIME   = f1->time
;   YYYYMM = cd_calendar(TIME, -1)
;   ; print(YYYYMM)
;   tStrt  = ind(YYYYMM .eq. (yrStrt*100 +  1))       ; indices 20th century
;   tLast  = ind(YYYYMM .eq. (yrLast*100 + 12))       ; relative to entire dataset

;     hadisst_tempt = f1->sst(tStrt:tLast,:,:) ; (time, lat, lon)
;     ice = f2->sic(tStrt:tLast,:,:)
;     hadisst       = where( ice.gt.0.0, hadisst_tempt@_FillValue,hadisst_tempt)
;     copy_VarCoords(hadisst_tempt, hadisst)
;     ntime_hadisst = dimsizes(hadisst&time)
;     printVarSummary(hadisst)
;     delete([/fil_ssta, f1, TIME/])
;---------------------- HadISST -------------------------------------
    fil_ssta      = "/data05/zhaobw/atlantic/netcdf/SST/hadisst/hadisst.ssta.mon.1900to2019.rmSeaIce.nc"
    f1   = addfile(fil_ssta, "r")
    fil_ice       = "/data05/zhaobw/atlantic/netcdf/SST/hadisst/HadISST.ice.mon.190001to201912.r360x180.tt.nc"
    f2   = addfile(fil_ice, "r")

  ; yrStrt  = 1900          ; manually specify for convenience
  ; yrLast  = 2019          ; 20th century ends 2010
  ; TIME   = f1->time
  ; YYYYMM = cd_calendar(TIME, -1)
  ; ; print(YYYYMM)
  ; tStrt  = ind(YYYYMM .eq. (yrStrt*100 +  1))       ; indices 20th century
  ; tLast  = ind(YYYYMM .eq. (yrLast*100 + 12))       ; relative to entire dataset

    hadisst = f1->ssta ; (time, lat, lon)
    ice = f2->sic;(tStrt:tLast,:,:)
    copy_VarCoords(ersst, hadisst)
    copy_VarAtts(ersst, hadisst)
    hadisst       = where( ice.gt.0.0, hadisst@_FillValue,hadisst)
    ntime_hadisst = dimsizes(hadisst&time)
    printVarSummary(hadisst)
    delete([/fil_ssta, f1/])
;----------------------  SODA  -------------------------------------
    fil_ssta      = "/data05/zhaobw/atlantic/netcdf/SST/soda/SODAsi.3_sst_1900to2013_tt.subclim1971.r360x180.nc"
    f1   = addfile(fil_ssta, "r")

  yrStrt  = 1900          ; manually specify for convenience
  yrLast  = 2013          ; 20th century ends 2010
  TIME   = f1->time
  ; YYYYMM = cd_calendar(TIME, -1)
    soda = f1->TEMP(:,0,:,:) ; (time, lat, lon)
    ntime_soda = dimsizes(soda&time)
    ; printVarSummary(hadisst)
    delete([/fil_ssta, f1, TIME, YYYYMM/])



;-------------------------------------topography for identify Atlantic Basin----------------------------------------
finland1            = addfile("/data05/linpf/FGOALS-g3/data/basin.r360x180.nc","r")
land41              = finland1->ind(:,:)
land_na             = where( (land41.ge.1.5.and.land41.le.2.8), 1, 0)
copy_VarCoords(land41, land_na)
printVarSummary(land41)
Atlantic_down = 1.5
Atlantic_up = 3
;----------------------------------------remove Mediterranean-----------------------------------------
  do j = 0, nlat-1
    do i= 0, nlon-1
      land_na(j,i)   = where( (i.ge.0.and.i.le.50).and.(j.ge.110.and.j.le.140),0,land_na(j,i))
    end do
  end do
;-------------------------------------- remove south Atlantic ------------------------------------
  do j = 0, nlat-1
    do i= 0, nlon-1
      land_na(j,i)   = where( j.ge.0.and.j.le.90,0,land_na(j,i))
    end do
  end do
;--------------------- remove strange place ------------------
do j = 0, nlat-1
    do i= 0, nlon-1
      land_na(j,i)   = where( i.ge.30.and.i.le.260,0,land_na(j,i))
    end do
  end do
; plot basin when switch.eq.1
switch = 0
if (switch.eq.1) then
  ; write branch
wks1 = gsn_open_wks("x11", get_script_prefix_name())
resa = True
resa@cnFillOn = True
resa@lbOrientation                        = "vertical"
; PLOT = gsn_csm_contour_map(wks1, land41, resa)
PLOT = gsn_csm_contour_map(wks1, land_na, resa)
end if

ersst_eof1dtrend_na = ersst_eof1dtrend
ersst_na = ersst
cobesst_na = cobesst
cobesst2_na = cobesst2
hadisst_na = hadisst
; kaplan_na = kaplan
soda_na = soda

;---- extract SSTA only over Atlantic region--------------------------
do t = 0, ntime_ersst_eof1dtrend-1
  ersst_eof1dtrend_na(t,:,:) = where( (land41.ge.1.5.and.land41.le.2.8), ersst_eof1dtrend(t,:,:), ersst_eof1dtrend@_FillValue)
end do
;---- extract SSTA only over Atlantic region--------------------------
do t = 0, ntime_ersst-1
  ersst_na(t,:,:) = where( (land41.ge.1.5.and.land41.le.2.8), ersst(t,:,:), ersst@_FillValue)
end do
;---- extract SSTA only over Atlantic region--------------------------
do t = 0, ntime_cobesst-1
  cobesst_na(t,:,:) = where( (land41.ge.1.5.and.land41.le.2.8) , cobesst(t,:,:), cobesst@_FillValue)
end do
;---- extract SSTA only over Atlantic region--------------------------
do t = 0, ntime_cobesst-1
  cobesst2_na(t,:,:) = where( (land41.ge.1.5.and.land41.le.2.8) , cobesst2(t,:,:), cobesst2@_FillValue)
end do
;---- extract SSTA only over Atlantic region--------------------------
do t = 0, ntime_hadisst-1
  ; hadisst_na(t,:,:) = where( (land41.ge.0.5.and.land41.le.2.8).or.(land41.ge.5.and.land41.le.5.5) , hadisst(t,:,:), hadisst@_FillValue)
  hadisst_na(t,:,:) = where( (land41.ge.1.5.and.land41.le.2.8) , hadisst(t,:,:), hadisst@_FillValue)
end do
;---- extract SSTA only over Atlantic region--------------------------
; do t = 0, ntime_soda-1
;   ; hadisst_na(t,:,:) = where( (land41.ge.0.5.and.land41.le.2.8).or.(land41.ge.5.and.land41.le.5.5) , hadisst(t,:,:), hadisst@_FillValue)
;   soda_na(t,:,:) = where( (land41.ge.1.5.and.land41.le.2.8) , soda(t,:,:), soda@_FillValue)
; end do
;-------------------------设置lanczos滤波器---------------------
    ihp=0       ; ihp=0 : low-pass filter
                ; ihp=1 : hi gh-pass filter
                ; ihp=2 : bandpass filter
    nWgt=61    ;
    nsigma=1.0  ;nsigma=1.0 is common
    fca=1./121.       ;因为数据的时间步长是年平均，要做10年的低通滤波:10
    fcb=-999.
    wgt=filwgts_lanczos(nWgt,ihp,fca,fcb,nsigma) ;120个月的滑动平均
;-------------------------------- load AMO from LFCA detrend method ---------------------------------
    dir          = ""
; 0-75n
    fil_S1       = "/data05/zhaobw/atlantic/netcdf/lfca_ersstv5/Atlantic_ersstv5_1900_2019_month_0_75N_local.nc"
    fil_S2       = "/data05/zhaobw/atlantic/netcdf/lfca_cobesst/Atlantic_cobesst_1900_2019_month_0_75N_local.nc"
    fil_S3       = "/data05/zhaobw/atlantic/netcdf/lfca_cobesst2/Atlantic_cobesst2_1900_2019_month_0_75N_local.nc"
    fil_S4       = "/data05/zhaobw/atlantic/netcdf/lfca_hadisst/Atlantic_hadisst_1900_2019_month_0_75N_local_seaice.nc"
    ; fil_S4       = "/data05/zhaobw/atlantic/netcdf/lfca_hadisst_nomask/Atlantic_hadisst_1900_2019_month_0_75N_local.nc"
    ; ersst eof_detrend
    fil_S5       = "/data05/zhaobw/atlantic/netcdf/eof_ersstv5_lfca/Atlantic_ersstv5_1900_2019_month_0_75N_local_detrendeof1_lfca_oriregress.nc" 

    f_s6       = addfile("seaIceMask_region.nc", "r")

    f_s1         = addfile(dir+fil_S1, "r")
    f_s2         = addfile(dir+fil_S2, "r")
    f_s3         = addfile(dir+fil_S3, "r")
    f_s4         = addfile(dir+fil_S4, "r")
    f_s5         = addfile(dir+fil_S5, "r")

    ssta_s1_detrend     = f_s1->LFC2
    ssta_s2_detrend     = f_s2->LFC2
    ssta_s3_detrend     = f_s3->LFC2
    ssta_s4_detrend     = f_s4->LFC3 ; hadisst
    ssta_s5_detrend     = f_s5->LFC1 ; hadisst

    topo                = f_s6->topo
    copy_VarCoords(land_na, topo)
    printVarSummary(topo)
; switch = 0
; if (switch.eq.1) then
;   ; write branch
; wks1 = gsn_open_wks("x11", get_script_prefix_name())
; resa = True
; resa@cnFillOn = True
; resa@lbOrientation                        = "vertical"
; ; PLOT = gsn_csm_contour_map(wks1, land41, resa)
; PLOT = gsn_csm_contour_map(wks1, topo, resa)
; end if
    topo_conform  =conform(ersst_na, topo, (/1,2/))
    ; mask topography
    ersst_na = where(topo_conform.ne.0, ersst_na, ersst_na@_FillValue)
    cobesst_na = where(topo_conform.ne.0, cobesst_na, cobesst_na@_FillValue)
    cobesst2_na = where(topo_conform.ne.0, cobesst2_na, cobesst2_na@_FillValue)
    hadisst_na = where(topo_conform.ne.0, hadisst_na, hadisst_na@_FillValue)
    ersst_eof1dtrend_na = where(topo_conform.ne.0, ersst_eof1dtrend_na, ersst_eof1dtrend_na@_FillValue)


    hadi_lfp3 = f_s4->LFP3
    copy_VarCoords(ersst_na(0,:,:), hadi_lfp3)
    hadi_lfp3&lon = fspan(0, 360, 360)
    hadi_lfp3@_FillValue = -9.99e36
    hadi_lfp3 = where( (land41.ge.1.5.and.land41.le.2.8) , hadi_lfp3, hadi_lfp3@_FillValue )


    lfp1 = f_s5->LFP1
    copy_VarCoords(ersst_na(0,:,:), lfp1)
    lfp1&lon = fspan(0, 360, 360)
    lfp1@_FillValue = -9.99e36

    ; lfp1       = 

    lfp1_eofdtrend = where( (land41.ge.1.5.and.land41.le.2.8) , lfp1, lfp1@_FillValue )
    copy_VarCoords(ersst_na(0,:,:), lfp1_eofdtrend)
    lfp1_eofdtrend&lon = fspan(0, 360, 360)
;------------------ filtering with 10 years ----------------------------------
    ssta_s1_detrend_filt = wgt_runave(ssta_s1_detrend, wgt, 0)
    ssta_s2_detrend_filt = wgt_runave(ssta_s2_detrend, wgt, 0)
    ssta_s3_detrend_filt = wgt_runave(ssta_s3_detrend, wgt, 0)
    ssta_s4_detrend_filt = wgt_runave(ssta_s4_detrend, wgt, 0)
    ssta_s5_detrend_filt = wgt_runave(ssta_s5_detrend, wgt, 0)
;------------------ standard ----------------------------------
    amo_s1 = dim_standardize(ssta_s1_detrend_filt,0)
    amo_s2 = dim_standardize(ssta_s2_detrend_filt,0)
    amo_s3 = dim_standardize(ssta_s3_detrend_filt,0)
    amo_s4 = dim_standardize(ssta_s4_detrend_filt,0)
    amo_s5 = dim_standardize(ssta_s5_detrend_filt,0)
;---------------------- n1---------------------------------------------------------------------
    amo_s1_rc = regCoef_n(amo_s1, dim_standardize_n_Wrap(ersst_na, 0, 0), 0, 0)
; 5 applications of the smth9 function
    ; amo_s1_rc  = (/ smth9(amo_s1_rc_tempt,0.5,0,True) /)
    ; printVarSummary(amo_n1_rc) ; smth9 destroy attribute of amo_rc
    copy_VarCoords(ersst_na(0,:,:), amo_s1_rc)
    amo_s1_rc&lon  = fspan(0, 360, 360)

   df          = onedtond(amo_s1_rc@nptxy,dimsizes(amo_s1_rc))-2   ; degrees of freedom
   tval        = onedtond(amo_s1_rc@tval,dimsizes(amo_s1_rc))      ; t-statistic
   b           = tval
   b           = 0.5          ; b must be same size as tval (and df)
   prob        = betainc(df/(df+tval^2),df/2.0,b)
   prob_amo_s1_rc  = (/prob/)
   copy_VarCoords(amo_s1_rc,prob_amo_s1_rc)
   printVarSummary(prob_amo_s1_rc)
   delete([/df,tval,b,prob/])
;---------------------- n2---------------------------------------------------------------------
    ; printVarSummary(amo_s2)
    ; printVarSummary(cobesst_na)
    amo_s2_rc = regCoef_n(amo_s2, dim_standardize_n_Wrap(cobesst_na, 0, 0), 0, 0)
    ; printVarSummary(amo_n1_rc) ; smth9 destroy attribute of amo_rc
    copy_VarCoords(cobesst_na(0,:,:), amo_s2_rc)
    amo_s2_rc&lon  = fspan(0, 360, 360)
    df          = onedtond(amo_s2_rc@nptxy,dimsizes(amo_s2_rc))-2   ; degrees of freedom
   tval        = onedtond(amo_s2_rc@tval,dimsizes(amo_s2_rc))      ; t-statistic
   b           = tval
   b           = 0.5          ; b must be same size as tval (and df)
   prob        = betainc(df/(df+tval^2),df/2.0,b)
   prob_amo_s2_rc  = (/prob/)
   copy_VarCoords(amo_s2_rc,prob_amo_s2_rc)
   printVarSummary(prob_amo_s2_rc)
   delete([/df,tval,b,prob/])
;---------------------- n3---------------------------------------------------------------------
    amo_s3_rc = regCoef_n(amo_s3, dim_standardize_n_Wrap(cobesst2_na, 0, 0), 0, 0)
    ; printVarSummary(amo_n1_rc) ; smth9 destroy attribute of amo_rc
    copy_VarCoords(cobesst2_na(0,:,:), amo_s3_rc)
    amo_s3_rc&lon  = fspan(0, 360, 360)

    df          = onedtond(amo_s3_rc@nptxy,dimsizes(amo_s3_rc))-2   ; degrees of freedom
   tval        = onedtond(amo_s3_rc@tval,dimsizes(amo_s3_rc))      ; t-statistic
   b           = tval
   b           = 0.5          ; b must be same size as tval (and df)
   prob        = betainc(df/(df+tval^2),df/2.0,b)
   prob_amo_s3_rc  = (/prob/)
   copy_VarCoords(amo_s3_rc,prob_amo_s3_rc)
   printVarSummary(prob_amo_s3_rc)
   delete([/df,tval,b,prob/])
;---------------------- n4---------------------------------------------------------------------
    amo_s4_rc = regCoef_n(amo_s4, dim_standardize_n_Wrap(hadisst_na, 0, 0), 0, 0)
    ; printVarSummary(amo_n1_rc) ; smth9 destroy attribute of amo_rc
    copy_VarCoords(hadisst_na(0,:,:), amo_s4_rc)
    amo_s4_rc&lon  = fspan(0, 360, 360)

    df         = onedtond(amo_s4_rc@nptxy,dimsizes(amo_s4_rc))-2   ; degrees of freedom
   tval        = onedtond(amo_s4_rc@tval,dimsizes(amo_s4_rc))      ; t-statistic
   b           = tval
   b           = 0.5          ; b must be same size as tval (and df)
   prob        = betainc(df/(df+tval^2),df/2.0,b)
   prob_amo_s4_rc  = (/prob/)
   copy_VarCoords(amo_s4_rc,prob_amo_s4_rc)
   copy_VarAtts(amo_s4_rc, prob_amo_s4_rc)
   ; prob_amo_s4_rc       = where( ice.gt.0.0, prob_amo_s4_rc@_FillValue,prob_amo_s4_rc)
   printVarSummary(prob_amo_s4_rc)
   delete([/df,tval,b,prob/])
;---------------------- n5---------------------------------------------------------------------
    amo_s5_rc = regCoef_n(amo_s5, dim_standardize_n_Wrap(ersst_eof1dtrend_na, 0, 0), 0, 0)
    ; printVarSummary(amo_n1_rc) ; smth9 destroy attribute of amo_rc
    copy_VarCoords(ersst(0,:,:), amo_s5_rc)
    copy_VarAtts(ersst(0,:,:), amo_s5_rc)
    amo_s5_rc&lon  = fspan(0, 360, 360)

   df          = onedtond(amo_s5_rc@nptxy,dimsizes(amo_s5_rc))-2   ; degrees of freedom
   tval        = onedtond(amo_s5_rc@tval,dimsizes(amo_s5_rc))      ; t-statistic
   b           = tval
   b           = 0.5          ; b must be same size as tval (and df)
   prob        = betainc(df/(df+tval^2),df/2.0,b)
   prob_amo_s5_rc  = (/prob/)
   copy_VarCoords(amo_s5_rc,prob_amo_s5_rc)
   printVarSummary(prob_amo_s5_rc)
   delete([/df,tval,b,prob/])
; ;---------------------- n5---------------------------------------------------------------------
    ; amo_s5_rc = regCoef_n(amo_s5, dim_standardize_n_Wrap(kaplan_na, 0, 0), 0, 0)
    ; ; printVarSummary(amo_n1_rc) ; smth9 destroy attribute of amo_rc
    ; copy_VarCoords(kaplan_na(0,:,:), amo_s5_rc)

;------------- pattern correlation -------------------------------------------------------------------
coswgt = cos(0.01745329*lat)
; ; erv5 & cobe/cobe2/hadi
; corr_erv5_cobe = pattern_cor(amo_s1_rc, amo_s2_rc, coswgt, 0)
; corr_erv5_cobe2 = pattern_cor(amo_s1_rc, amo_s3_rc, coswgt, 0)
; corr_erv5_hadi = pattern_cor(amo_s1_rc, amo_s4_rc, coswgt, 0)
; ;hadi & cobe/cobe2
; corr_hadi_cobe = pattern_cor(amo_s4_rc, amo_s2_rc, coswgt, 0)
; corr_hadi_cobe2 = pattern_cor(amo_s4_rc, amo_s3_rc, coswgt, 0)
; ; cobe & cobe2
; corr_cobe_cobe2 = pattern_cor(amo_s2_rc, amo_s3_rc, coswgt, 0)
amo_s1_rc_na = where(land_na.eq.1, amo_s1_rc, amo_s1_rc@_FillValue)
amo_s2_rc_na = where(land_na.eq.1, amo_s2_rc, amo_s2_rc@_FillValue)
amo_s3_rc_na = where(land_na.eq.1, amo_s3_rc, amo_s3_rc@_FillValue)
amo_s4_rc_na = where(land_na.eq.1, hadi_lfp3, hadi_lfp3@_FillValue)
; erv5 & cobe/cobe2/hadi
corr_erv5_cobe = pattern_cor(amo_s1_rc_na, amo_s2_rc_na, coswgt, 1); 1: uncentered patetrn correlation
corr_erv5_cobe2 = pattern_cor(amo_s1_rc_na, amo_s3_rc_na, coswgt, 1)
corr_erv5_hadi = pattern_cor(amo_s1_rc_na, amo_s4_rc_na, coswgt, 1)
;hadi & cobe/cobe2
corr_hadi_cobe = pattern_cor(amo_s4_rc_na, amo_s2_rc_na, coswgt, 1)
corr_hadi_cobe2 = pattern_cor(amo_s4_rc_na, amo_s3_rc_na, coswgt, 1)
; cobe & cobe2
corr_cobe_cobe2 = pattern_cor(amo_s2_rc_na, amo_s3_rc_na, coswgt, 1)
;--------------------- significance test for PCC -----------------------------------------------------------
amo_s1_rc_spaceseries = ndtooned(amo_s1_rc_na); ERv5
amo_s2_rc_spaceseries = ndtooned(amo_s2_rc_na); COBE
amo_s3_rc_spaceseries = ndtooned(amo_s3_rc_na); COBE2
amo_s4_rc_spaceseries = ndtooned(amo_s4_rc_na); HadI
;                    N_eff = N*(1-rx*ry)/(1+rx*ry), following Bretherton et al. [1999]
;                    rx: lag-1 autocorrelation of variable x
;                    ry: lag-1 autocorrelation of variable y
tempt = esacr(amo_s1_rc_spaceseries, 12)
r1 = tempt(1)
tempt = esacr(amo_s2_rc_spaceseries, 12)
r2 = tempt(1)
tempt = esacr(amo_s3_rc_spaceseries, 12)
r3 = tempt(1)
tempt = esacr(amo_s4_rc_spaceseries, 12)
r4 = tempt(1)
n = dimsizes(amo_s1_rc_spaceseries)
N_eff_r1_r2 = n*(1-r1*r2)/(1+r1*r2)
N_eff_r1_r3 = n*(1-r1*r3)/(1+r1*r3)
N_eff_r1_r4 = n*(1-r1*r4)/(1+r1*r4)

print(tempt)
print(N_eff_r1_r2)
print("effective degree of freedom: N(ERv5,COBE): "+N_eff_r1_r2+", N(ERv5,COBE2): "+N_eff_r1_r3+", N(ERv5,HadI): "+N_eff_r1_r4)

print(" PCC: r(ERv5,COBE)="+sprintf("%4.2f", corr_erv5_cobe)+", r(ERv5,COBE2)="+sprintf("%4.2f", corr_erv5_cobe2)+", r(ERv5,HadI)="+sprintf("%4.2f", corr_erv5_hadi)+", r(HadI,COBE)="+sprintf("%4.2f", corr_hadi_cobe)+", r(HadI,COBE2)="+sprintf("%4.2f", corr_hadi_cobe2)+", r(COBE,COBE2)="+sprintf("%4.2f", corr_cobe_cobe2))

switch = 1
if ( switch.eq.1) then
;*********************************
; plot figures
;*********************************
  plot   = new ( 5, "graphic")
  wks = gsn_open_wks("eps", get_script_prefix_name)
    gsn_define_colormap(wks, "significance_level_zbw") ; add my personal colorbar only used for plot significant pass colorbar
     
; - create the plots
   res                                = True
   res@gsnDraw                        = False
   res@gsnFrame                       = False
   res@gsnAddCyclic                   = False
   res@gsnMaximize                    = False
   ;res@mpProjection                   = "WinkelTripel"
   res@mpGeophysicalLineColor         = "gray42"
   res@mpPerimOn                      = False
   ; res@mpGridLatSpacingF              = 45             ; change latitude  line spacing
   ; res@mpGridLonSpacingF              = 60.            ; change longitude line spacing
   res@mpGridLineColor                = "transparent"  ; trick ncl into drawing perimeter
   res@mpGridAndLimbOn                = True           ; turn on lat/lon lines
   res@mpFillOn                       = True
   res@mpCenterLonF                   = -60
   ;res@mpMaxLonF                      = 30
   ;res@mpMinLonF                      = -150

   res@mpOutlineOn                    = True ;drawing of MapPlot area outlines. If True, the MapPlot enables drawing of map area outlines
   ; res@mpLandFillColor                = "White"
   ; res@mpDataSetName                  = "Earth..4"
   res@mpFillDrawOrder                = "PreDraw"
   ; res@cnLevelSelectionMode           = "ExplicitLevels"
   ; res@cnLevels                       = fspan(-0.36,0.36,24+1)
   res@cnLevelSelectionMode           = "ExplicitLevels"
   res@cnLevels                       =  (/-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5/)
   ; res@lbLabelStrings                 = (/"-0.5","-0.4","-0.3","-0.2","-0.1","0.1","0.2","0.3","0.4","0.5"/)
   res@cnLineLabelsOn                 = False
   res@cnFillOn                       = True
   res@cnLinesOn                      = False
   res@cnFillDrawOrder                = "PreDraw"
   res@lbLabelBarOn                   = False
   ;res@gsnRightStringParallelPosF     = 0.94
   ;res@gsnLeftStringOrthogonalPosF    = 0.0
   ;res@gsnLeftStringParallelPosF      = 0.01

   ;res@gsnStringFont                  = 0
   res@gsnStringFontHeightF           = 0.020
   res@tiMainFont                     = 0
   res@tiMainFontHeightF              = 0.030
   res@gsnRightString                 = ""

   res@tiMainString                   = ""
   res@gsnLeftString                  = ""
   res@gsnStringFont                  = "times-roman"
   res@tmYLLabelFont                  = "times-roman"
   res@tmXBLabelFont                  = "times-roman"  
   res@mpMaxLatF                      = 80
   res@mpMinLatF                      = 0  
   res@mpMaxLonF                      = 30
   res@mpMinLonF                      = -90
   res@mpCenterLonF                   = -30
   ; res@mpMaxLatF                      =  
;------------- significance level resource -----------------------
   res1                               = True
   res1@gsnDraw                       = False                     ; do not draw
   res1@gsnFrame                      = False                     ; do not advance frame
   res1@gsnAddCyclic                  = False
   res1@gsnMaximize                   = False
   res1@cnMonoFillPattern             = False
   res1@cnLevelSelectionMode          = "ManualLevels"
   res1@cnMinLevelValF                =  0              ; set min contour level
   res1@cnMaxLevelValF                =  0.05           ; set max contour level
   res1@cnLevelSpacingF               =  0.05           ; set contour interval
   res1@cnFillDrawOrder               =  "PreDraw"
   res1@gsnLeftString                 =  ""
   res1@cnLineLabelsOn                = False
   res1@cnLinesOn                     = False
   res1@cnInfoLabelOn                 = False           ; turn off contour label

;------------- significance level resource version2 -----------------------
    res2                              =True
    res2@gsnFrame=False
    res2@gsnDraw=False
    res2@gsnAddCyclic                 = False
    res2@cnFillOn                     = True
    res2@cnLinesOn                    =False
    res2@cnLineLabelsOn               =False;不要等值线上的标签
    res2@cnLevelSelectionMode         = "ExplicitLevels"
    res2@cnLevels= (/0.05/)
    res2@cnMonoFillPattern            =False;取消默认填充形状
    res2@cnMonoFillScale=False;取消默认填充密度
    res2@cnFillPatterns= (/17,-1/) ; the patterns
    res2@cnFillColors=(/"black","white"/)
    res2@lbLabelBarOn=False
    res2@cnInfoLabelOn= False
    res2@cnFillScales= (/0.8/) ;数值越大越稀疏
    res2@tiMainString=""
    res2@cnFillOpacityF=1.0
;---------------- shade resource --------------------------
   opt                                =   True
   opt@gsnShadeFillType               =    "pattern"
   opt@gsnShadeMid                    =    17

sig = 0.01
        res@gsnLeftStringFontHeightF  = 0.03
        res@gsnLeftStringOrthogonalPosF = 0.035
        res@gsnCenterStringFontHeightF= 0.03
        res@gsnRightStringFontHeightF = 0.03
    ; res@tmXBLabels = (/"","","180~S~o~N~W","","","90~S~o~N~W","","","0","","","90~S~o~N~E",""/)
    res@tmXBLabelFontHeightF          = 0.03

        ; res@tmYLLabels                = (/"","60~S~o~N~S","","0","","60~S~o~N~N",""/)
    res@tmYLLabelFontHeightF           = 0.03
;--------------- resource for only plot positive/negative contour line for boundary
    res_bound                          = True
    res_bound@gsnDraw                  = False
    res_bound@gsnFrame                 = False
    res_bound@cnFillOn                 = False
    res_bound@cnLineLabelsOn           = True
    ; res_bound@cnLabelMasking           = True
    res_bound@cnMonoLineThickness      = False
    res_bound@cnLinesOn                = True
    res_bound@cnLevelSelectionMode     = "ExplicitLevels"
    res_bound@cnLevels                 = (/-0.1,0.1/)
    res_bound@cnLineThicknesses        = (/2,2/)
    res_bound@cnInfoLabelOn            = False

    res@tiYAxisOffsetXF               = -0.05
    ; res@tiYAxisOffsetYF               = 0.1
    res@tiYAxisFontHeightF = 0.05
; ;-------------- amo_n1 ----------------
; res@tiYAxisString   = "0-75~S~o~N~N"
    res@gsnLeftString = "(a) ERv5"
    res@gsnRightString = "69.2%"
      plot(0) = gsn_csm_contour_map(wks,amo_s1_rc,res)
      plot1 = gsn_csm_contour(wks, prob_amo_s1_rc, res2)
      overlay(plot(0), plot1)
res@tiYAxisString   = ""
; ;-------------- amo_n1 ----------------
    res@gsnLeftString = "(b) COBE"
    res@gsnRightString = "98.9%"
      ; plot(1) = gsn_csm_contour_map(wks,amo_s2_rc,res)
      ; plot2 = gsn_csm_contour(wks, prob_amo_s2_rc, res2)
      ; overlay(plot(1), plot2)
; ;-------------- amo_n1 ----------------
    res@gsnRightString = "94%";"91.6%"
    res@gsnLeftString = "(c) COBE2"
      plot(2) = gsn_csm_contour_map(wks,amo_s3_rc,res)
      plot3 = gsn_csm_contour(wks, prob_amo_s3_rc, res2)
      overlay(plot(2), plot3)
; ;-------------- amo_n1 ----------------
    res@gsnLeftString = "(b) HadI1"
    res@gsnRightString = "93.5%"
      plot(1) = gsn_csm_contour_map(wks,hadi_lfp3,res); pay attention to hadi_lfp3 rather than amo_s4_rc
      plot4 = gsn_csm_contour(wks, prob_amo_s4_rc, res2)
      overlay(plot(1), plot4)
; ;-------------- amo_n5 ----------------
    res@gsnRightString = "41.6%"
    res@gsnLeftString = "(d) EOF1_dtrend"
      plot(3) = gsn_csm_contour_map(wks,amo_s5_rc,res);lfp1_eofdtrend,res)
      plot5 = gsn_csm_contour(wks, prob_amo_s5_rc, res2)
      overlay(plot(3), plot5)


;---------------------------------------------------------------------------------------
; - create the panel
   resp                               = True
   resp@gsnPaperOrientation           = "portrait"
   ; resp@gsnPanelRight                 = 0.9 ; 0<=   <=1.0
   ; resp@gsnPanelLeft    = 0.1 ; 0<=   <=1.0
  ;resp@gsnPanelDebug = True
  resp@gsnPanelMainString             = ""
  resp@gsnPanelMainFontHeightF        = 0.025
  resp@gsnPanelMainFont               = "times-roman"
   resp@gsnMaximize                   = True
   ; resp@gsnPanelLeft                  = 0.1
   ; resp@gsnPanelRight                 = 0.95
   ;resp@gsnPanelBottom                = 0.05
   resp@gsnPanelLabelBar              = True
   resp@lbOrientation                 = "Vertical"
   ; resp@pmLabelBarHeightF             = 0.4
   ; resp@pmLabelBarWidthF              = 0.05
   resp@lbLabelFont                   = "times-roman"
   resp@lbLabelFontThicknessF         = 1.5
   resp@lbLabelFontHeightF            = 0.016
   resp@pmLabelBarWidthF              = 0.035
   resp@pmLabelBarHeightF             = 0.45

   ; resp@lbLabelStrings                = (/"-0.36","","-0.30","","-0.24","","-0.18","","-0.12","","-0.06","","0","","0.06","","0.12","","0.18","","0.24","","0.30","","0.36"/)
   ; resp@lbLabelStrings                = (/"-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0.1","0.2","0.3","0.4","0.5","0.6"/)


   resp@lbTitleOn                     = False
   ;resp@lbBoxLinesOn                  = True
   ;resp@lbBoxLineColor                = "azure4"
   ;resp@lbBoxSeparatorLinesOn         = False
   ;resp@gsnPanelYF                    = (/0.95,0.74,0.53,0.32/)
   resp@lbBoxEndCapStyle                     = "TriangleBothEnds"

   gsn_panel(wks,plot,(/2,2/),resp)

   ; frame(wks)

   print("end of the script")

end if
end