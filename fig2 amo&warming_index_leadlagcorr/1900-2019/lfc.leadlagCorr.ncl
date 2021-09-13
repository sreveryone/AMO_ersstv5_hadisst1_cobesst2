load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************** Main code begin ***************************************************************************
begin

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
    ; fil_S5       = "/data05/zhaobw/atlantic/netcdf/lfca_kaplanv2/Atlantic_kaplanv2_1900_2019_month_eof30_cutoff10_0_75N_local.nc"
    ; fil_S6       = "/data05/zhaobw/atlantic/netcdf/lfca_soda/Atlantic_soda_1900_2013_month_eof30_cutoff10_0_75N_local.nc"

    fil_ersst    = "/data05/zhaobw/atlantic/netcdf/SST/ersstv5/ersstv5.monmean.190001to201912.subclim1971.r360x180.nc"
    ; fil_soda     = "/data05/zhaobw/atlantic/netcdf/SST/soda/SODAsi.3_sst_1900to2013_tt.subclim1971.r360x180.nc"

    f_s1         = addfile(dir+fil_S1, "r")
    f_s2         = addfile(dir+fil_S2, "r")
    f_s3         = addfile(dir+fil_S3, "r")
    f_s4         = addfile(dir+fil_S4, "r")
    ; f_s5         = addfile(dir+fil_S5, "r")
    ; f_s6         = addfile(dir+fil_S6, "r")
    f_ersst      = addfile(dir+fil_ersst, "r")
    ; f_soda      = addfile(dir+fil_soda, "r")


    ssta_s1_detrend     = f_s1->LFC2
    ssta_s2_detrend     = f_s2->LFC2
    ssta_s3_detrend     = f_s3->LFC2
    ssta_s4_detrend     = f_s4->LFC3
    ; ssta_s5_detrend     = f_s5->LFC2
    ; ssta_s6_detrend     = f_s6->LFC2

    ersst               = f_ersst->sst
    ; soda               = f_soda->TEMP(:,0,:,:)
    ; printVarSummary(ersst)
; print(ssta_s1_detrend)
;------------------ filtering with 10 years ----------------------------------
    ssta_s1_detrend_filt = wgt_runave(ssta_s1_detrend, wgt, 0)
    ssta_s2_detrend_filt = wgt_runave(ssta_s2_detrend, wgt, 0)
    ssta_s3_detrend_filt = wgt_runave(ssta_s3_detrend, wgt, 0)
    ssta_s4_detrend_filt = wgt_runave(ssta_s4_detrend, wgt, 0)
    ; ssta_s5_detrend_filt = wgt_runave(ssta_s5_detrend, wgt, 0)
    ; ssta_s6_detrend_filt = wgt_runave(ssta_s6_detrend, wgt, 0)

;------------------ standard ----------------------------------
    amo_s1 = dim_standardize(ssta_s1_detrend_filt,0)
    amo_s2 = dim_standardize(ssta_s2_detrend_filt,0)
    amo_s3 = dim_standardize(ssta_s3_detrend_filt,0)
    amo_s4 = dim_standardize(ssta_s4_detrend_filt,0)
    ; amo_s5 = dim_standardize(ssta_s5_detrend_filt,0)
    ; amo_s6 = dim_standardize(ssta_s6_detrend_filt,0)

amo_s1!0 = "time"
amo_s1&time = fspan(1900, 2019, (2019-1900+1)*12)
; printVarSummary(amo_s1)
copy_VarCoords(amo_s1, amo_s2)
copy_VarCoords(amo_s1, amo_s3)
copy_VarCoords(amo_s1, amo_s4)
; copy_VarCoords(amo_s1, amo_s5)
; amo_s6!0 = "time"
; amo_s6&time = fspan(1900, 2013, (2013-1900+1)*12)


switch = 1
if (switch.eq.1) then
;--------------------------------------------------------------------------------
;                   cross-correlation at 12-month lead-lag
;--------------------------------------------------------------------------------
     mxlag         = 25*12
; lead-lag correlation of s2
     amo_s1_Lead_amo_s2 = esccr(amo_s1,amo_s2,mxlag)
     amo_s2_Lead_amo_s1 = esccr(amo_s2,amo_s1,mxlag)    ; switch the order of the series

     leadlag_amo_s1_s2 = new ( 2*mxlag+1, float)    
     leadlag_amo_s1_s2(0:mxlag-1) = amo_s1_Lead_amo_s2(1:mxlag:-1)  ; "negative lag", -1 reverses order
     leadlag_amo_s1_s2(mxlag:)    = amo_s2_Lead_amo_s1(0:mxlag)     ; "positive lag"

; lead-lag correlation of s3
     amo_s1_Lead_amo_s3 = esccr(amo_s1,amo_s3,mxlag)
     amo_s3_Lead_amo_s1 = esccr(amo_s3,amo_s1,mxlag)    ; switch the order of the series

     leadlag_amo_s1_s3 = new ( 2*mxlag+1, float)    
     leadlag_amo_s1_s3(0:mxlag-1) = amo_s1_Lead_amo_s3(1:mxlag:-1)  ; "negative lag", -1 reverses order
     leadlag_amo_s1_s3(mxlag:)    = amo_s3_Lead_amo_s1(0:mxlag)     ; "positive lag"
; lead-lag correlation of s4
     amo_s1_Lead_amo_s4 = esccr(amo_s1,amo_s4,mxlag)
     amo_s4_Lead_amo_s1 = esccr(amo_s4,amo_s1,mxlag)    ; switch the order of the series

     leadlag_amo_s1_s4 = new ( 2*mxlag+1, float)    
     leadlag_amo_s1_s4(0:mxlag-1) = amo_s1_Lead_amo_s4(1:mxlag:-1)  ; "negative lag", -1 reverses order
     leadlag_amo_s1_s4(mxlag:)    = amo_s4_Lead_amo_s1(0:mxlag)     ; "positive lag"
; lead-lag correlation of s5
;      amo_s1_Lead_amo_s5 = esccr(amo_s1,amo_s5,mxlag)
;      amo_s5_Lead_amo_s1 = esccr(amo_s5,amo_s1,mxlag)    ; switch the order of the series

;      leadlag_amo_s1_s5 = new ( 2*mxlag+1, float)    
;      leadlag_amo_s1_s5(0:mxlag-1) = amo_s1_Lead_amo_s5(1:mxlag:-1)  ; "negative lag", -1 reverses order
;      leadlag_amo_s1_s5(mxlag:)    = amo_s5_Lead_amo_s1(0:mxlag)     ; "positive lag"
; ; lead-lag correlation of s6
;      amo_s1_Lead_amo_s6 = esccr(amo_s1(0:(2013-1900+1)*12-1),amo_s6,mxlag)
;      amo_s6_Lead_amo_s1 = esccr(amo_s6,amo_s1(0:(2013-1900+1)*12-1),mxlag)    ; switch the order of the series

;      leadlag_amo_s1_s6 = new ( 2*mxlag+1, float)    
;      leadlag_amo_s1_s6(0:mxlag-1) = amo_s1_Lead_amo_s6(1:mxlag:-1)  ; "negative lag", -1 reverses order
;      leadlag_amo_s1_s6(mxlag:)    = amo_s6_Lead_amo_s1(0:mxlag)     ; "positive lag"
;----------------------------------------------------------------------------------------;
; 寻找max并且标记
;----------------------------------------------
max_amo_s1_s2 = max(leadlag_amo_s1_s2)
do i = 0, dimsizes(leadlag_amo_s1_s2)-1
  if (leadlag_amo_s1_s2(i).ge.max_amo_s1_s2) then
    max_amo_s1_s2_position = i  
  end if
end do
max_amo_s1_s3 = max(leadlag_amo_s1_s3)
do i = 0, dimsizes(leadlag_amo_s1_s3)-1
  if (leadlag_amo_s1_s3(i).ge.max_amo_s1_s3) then
    max_amo_s1_s3_position = i  
  end if
end do
max_amo_s1_s4 = max(leadlag_amo_s1_s4)
do i = 0, dimsizes(leadlag_amo_s1_s4)-1
  if (leadlag_amo_s1_s4(i).ge.max_amo_s1_s4) then
    max_amo_s1_s4_position = i  
  end if
end do
; max_amo_s1_s5 = max(leadlag_amo_s1_s5)
; do i = 0, dimsizes(leadlag_amo_s1_s5)-1
;   if (leadlag_amo_s1_s5(i).ge.max_amo_s1_s5) then
;     max_amo_s1_s5_position = i  
;   end if
; end do
; max_amo_s1_s6 = max(leadlag_amo_s1_s6)
; do i = 0, dimsizes(leadlag_amo_s1_s6)-1
;   if (leadlag_amo_s1_s6(i).ge.max_amo_s1_s6) then
;     max_amo_s1_s6_position = i  
;   end if
; end do

; Calculate effective number of freedom (Bretherton et al., 1999)
; Algorithm:
; N = dimsizes(lfc1)
; N_eff = N*(1-rx*ry)/(1+rx*ry)
; rx = 1-lag autocorrelation of timeseries x
; ry = 1-lag autocorrelation of timeseries y
; N: length of timeseries
;-------------------------------------------------------------------------------------
mxlag = 1
N = dimsizes(amo_s1)
; 1-lag autocorrelation of TNH
ry_au = esccr(amo_s1, amo_s1, 1)
ry    = ry_au(1)

; 1-lag autocorrelation of blob_warm
rx_au = esccr(amo_s2, amo_s2, 1)
rx    = rx_au(1)
; effective degree of freedom
N_eff_blob_warm = N*(1-rx*ry)/(1+rx*ry)
; print("effective degree of freedom of AMO-ersst is "+N_eff_blob_warm+", "+" origional is "+N)
delete([/rx_au,rx/])

;----------------------------------------------------------------------------------;
; test leadlag correlation significance 
;-----------------------------------------------------------------------------------;
  ; siglvl = 0.05
  ; n      = new(dimsizes(leadlag_amo_s1_s2), "integer", floattointeger(N_eff_blob_warm))
  ; ; print(n)
  ; r      = leadlag_amo_s1_s2
  ; tval   = r/sqrt((1.0-r^2)/(n-2))
  ; prob   = rtest(r,n,0)
  ; yes_no = where(prob.lt.siglvl, True, False)
  ; ; print(r+"   "+n+"   "+tval+"   "+prob+"   "+yes_no )
  ; delete([/r, tval, prob, yes_no/])
  ; sig_lfc1_blob_warm = new(dimsizes(leadlag_blob_warm_lfc1), "float", -0.26)
  ; print(sig_lfc1_blob_warm)
  ; print(r+"   "+n+"   "+tval+"   "+prob+"   "+yes_no ) ; output to examine the correlation at significance level
  ; sig_lfc3_blob_warm = new(dimsizes(leadlag_blob_warm_lfc1), "float", 0.26)
;-------------------------------------------------------------------------;
;   Save data 
;---------------------------------------------------------------------------;
mxlag         = 25
Lines1      = new((/3,dimsizes(leadlag_amo_s1_s2)/), typeof(leadlag_amo_s1_s2))
Lines1(0,:) = (/leadlag_amo_s1_s2/)
Lines1(1,:) = (/leadlag_amo_s1_s3/)
Lines1(2,:) = (/leadlag_amo_s1_s4/)
; Lines1(3,:) = (/leadlag_amo_s1_s5/)
; Lines1(4,:) = (/leadlag_amo_s1_s6/)


Lines1!1    = "time"
Lines1&time = fspan(-mxlag, mxlag, 24*mxlag+1)

plot       = new(5, "graphic")
wks        = gsn_open_wks("x11", get_script_prefix_name())

;*****************************************
;          Plot    time series
;*****************************************
  sres                      = True      ; plot mods desired
  sres@gsnFrame             = False     ; don't advance frame
  sres@gsnDraw              = False

sres@tmXTOn                          = False  ; Setting True turns on top tick marks; False turns them off.
sres@tmYROn                          = False; Setting True turns on right tick marks; False turns them off.
sres@gsnStringFont                     = "times-roman"
   sres@tmYLLabelFont                  = "times-roman"
   sres@tmXBLabelFont                  = "times-roman"  
   sres@tiXAxisFont                    = "times-roman"  
   sres@tiYAxisFont                    = "times-roman"  
; ;add background grid line 
; sres@tmXMajorGrid                    = True
; sres@tmYMajorGrid                    = True
; sres@tmXMajorGridThicknessF          = 1.0
; sres@tmYMajorGridThicknessF          = 1.0
; sres@tmXMajorGridLineDashPattern     = 2
; sres@tmYMajorGridLineDashPattern     = 2
; sres@tmXMajorGridLineColor           = "gray40"
; sres@tmYMajorGridLineColor           = "gray40"

  ; sres@trXMinF              = 1950     ; set X-axis max and mins
  ; sres@trXMaxF              = 2020
  ; sres@trYMinF              = -3.      ; set Y-axis max and mins
  ; sres@trYMaxF              =  7.      

  sres@tmYLLabelFontHeightF = 0.015  ; font height
  sres@tmXBLabelFontHeightF = 0.015 
  sres@gsnTickMarksOn       = True        ; turn off tickmarks

  sres@gsnYRefLine          = 0.          ; Y-value for reference line

  sres@gsnCenterString      = ""       ; no center string
  sres@tiXAxisOn            = True   ; turn off X-axis title 
; this controls the size and location of the second plot
  ; sres@vpXF                 = 0.2   
  ; sres@vpYF                 = 0.7

  sres@vpWidthF             = 0.6
  sres@vpHeightF            = 0.2
  sres@tiYAxisOn            = True
  sres@tiYAxisString        = "STD"    ; y-axis label 
  sres@tiXAxisString        = "yr"
  ;sres@tiYAxisFontHeightF   = 0.1     
  ;sres@tiMainString     = "LFC1" ; title
  sres@gsnLeftStringFontHeightF      = 0.015
  ; sres@gsnLeftString        = "(a) Index" 
  ; sres@gsnRightString        = "r(blob_warm,IPO)= "+sprintf("%4.2f",ccr_blob_warm_ipo)+" "+"r(blob_warm,AMO)= "+sprintf("%4.2f",ccr_blob_warm_amo)
  sres@gsnRightStringOrthogonalPosF  = 0.01
;************画第一个序列*************************************************************************************
  polyres                  = sres
  ;xyMonoDashPattern ,if true: all the curves in the associated data object will be drawn with the 
  ;dash pattern specified by the xyDashPattern resource
  polyres@xyMonoDashPattern    = True         ; 
  polyres@xyDashPattern        = 0        ; set dash pattern
  polyres@xyLineThicknesses      = (/6,6,6,6,6/)      ; Define line thicknesses 
  polyres@xyMonoDashPattern      = False
  polyres@pmLegendDisplayMode    = "Always"            ; turn on legend
  polyres@pmLegendSide           = "Top"               ; Change location of 
  polyres@pmLegendParallelPosF   = .77                 ; move units right
  polyres@pmLegendOrthogonalPosF = -0.47                ; move units down
  polyres@pmLegendWidthF         = 0.08                ; Change width and
  polyres@pmLegendHeightF        = 0.07                ; height of legend.
  polyres@lgPerimOn              = False               ; turn off box around
  polyres@lgLabelFontHeightF     = .015                 ; label font height

  polyres@tmXBMode             = "Explicit"
  ; polyres@tmXBValues           = (/-20,-15/)
  ; polyres@tmXBLabels           = (/"-18","-15","-12","-9","-6","-3","0","3","6","9","12","15","18"/)
  ; polyres@tmXBValues           = (/-240,-180,-120,-60,0,60,120,180,240/)
  ; polyres@tmXBLabels           = (/"-20","-15","-10","-5","0","5","10","15","20"/)

  polyres@trXMinF              = -20     ; set X-axis max and mins
  polyres@trXMaxF              = 20

  polyres@trYMinF              = 0;-.3      ; set Y-axis max and mins
  polyres@trYMaxF              =  1.      
  polyres@xyLineColors         = (/"forestgreen","blue","red","purple","gray11"/)
  polyres@lgLabelFont            = "times-roman"
  polyres@xyExplicitLegendLabels = (/" ERv5 & COBE"," ERv5 & COBE2"," ERv5 & HadI1"," ERv5 & Kaplanv2", " ERSST & SODA"/)
  ; delete(polyres@xyDashPatterns)
  ; polyres@xyDashPatterns         = (/0,0,0,0,0/)
  polyres@tiYAxisString        = "Corr-Coeff"    ; y-axis label
  polyres@tiXAxisString        = "year"
  polyres@gsnLeftString        = "(a) Lead-lag correlation of AMO index" 
  polyres@gsnRightString       = ""
; polyres@xyExplicitLegendLabels = (/"","", ""/)         ; create explicit labels
  plot(0) = gsn_csm_xy(wks, Lines1&time, Lines1, polyres) 
  txres                               = True
  txres@txFontHeightF                 = 0.012
  txres@txFont                        = "times-roman"
  txres@txFontThicknessF                    = 3  
;add 0.05 significance level
pres                                = True
; pres@gsLineDashPattern              = (/1/)
pres@gsLineThicknessF               = 3.0
pres@gsLineColor                    = "gray11"
;pres@gsLineDashPattern              = (/3/)
; dum1                                = gsn_add_polyline(wks,plot(0),(/-40, 40/),(/0.108,0.108/),pres)
; dum2                                = gsn_add_polyline(wks,plot(0),(/-40,40/),(/-0.108,-0.108/),pres)
dum3                                = gsn_add_polyline(wks,plot(0),(/0,0/),(/-10,10/),pres)
markeres = True
markeres@gsMarkerIndex = (/1,1,1,1,1/)
markeres@gsMarkerSizeF  = 0.08
markeres@gsMarkerColor = (/"forestgreen"/);(/"forestgreen","red","blue","purple","gray"/)
marker = gsn_add_polymarker(wks, plot(0), max_amo_s1_s2_position/12-mxlag, max_amo_s1_s2, markeres)
markeres@gsMarkerColor = (/"blue"/)
marker2 = gsn_add_polymarker(wks, plot(0), max_amo_s1_s3_position/12-mxlag, max_amo_s1_s3, markeres)
markeres@gsMarkerColor = (/"red"/)
marker3 = gsn_add_polymarker(wks, plot(0), max_amo_s1_s4_position/12-mxlag, max_amo_s1_s4, markeres)
; markeres@gsMarkerColor = (/"purple"/)
; marker4 = gsn_add_polymarker(wks, plot(0), max_amo_s1_s5_position-mxlag, max_amo_s1_s5, markeres)
; markeres@gsMarkerColor = (/"gray11"/)
; marker5 = gsn_add_polymarker(wks, plot(0), max_amo_s1_s6_position-mxlag, max_amo_s1_s6, markeres)

;------------------------------- add text in plot ---------------------------------------------
 txres               = True                            
 txres@txFontHeightF = 0.015             ; Set the font height

 text1  = gsn_add_text(wks,plot(0),"ERv5 lead",-9, 0.6,txres)
 text2  = gsn_add_text(wks,plot(0),"ERv5 leag", 9, 0.6,txres)


; - create the panel
   resp                               = True
   ; resp@lbOrientation                 = "Vertical"
   ; resp@pmLabelBarHeightF             = 0.25
   ; resp@pmLabelBarWidthF              = 0.05
   ; resp@lbLabelFont                   = 0
   ; resp@lbLabelFontHeightF            = 0.02
   resp@lbBoxEndCapStyle                     = "TriangleBothEnds"
   resp@gsnPanelMainString                 = ""
   gsn_panel(wks,plot,(/1,1/),resp)
   ; frame(wks)
   print("end of the script")
end if

end