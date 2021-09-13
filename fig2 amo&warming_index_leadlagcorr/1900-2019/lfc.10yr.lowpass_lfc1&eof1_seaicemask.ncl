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
    ; fil_S1       = "/data05/zhaobw/atlantic/netcdf/lfca_ersstv5/Atlantic_ersstv5_1900_2019_month_0_75N_local.nc"
    ; fil_S2       = "/data05/zhaobw/atlantic/netcdf/cobesst_lfca/Atlantic_cobesst_1900_2019_month_eof30_cutoff10_0_75N_local.nc"
    ; fil_S3       = "/data05/zhaobw/atlantic/netcdf/lfca_cobesst2/Atlantic_cobesst2_1900_2019_month_eof30_cutoff10_0_75N_local.nc"
    ; fil_S4       = "/data05/zhaobw/atlantic/netcdf/lfca_hadisst/Atlantic_hadisst_1900_2019_month_eof30_cutoff10_0_75N_local.nc"

    fil_S1       = "/data05/zhaobw/atlantic/netcdf/lfca_ersstv5/Atlantic_ersstv5_1900_2019_month_0_75N_local.nc"
    fil_S2       = "/data05/zhaobw/atlantic/netcdf/lfca_cobesst/Atlantic_cobesst_1900_2019_month_0_75N_local.nc"
    fil_S3       = "/data05/zhaobw/atlantic/netcdf/lfca_cobesst2/Atlantic_cobesst2_1900_2019_month_0_75N_local.nc"
    fil_S4       = "/data05/zhaobw/atlantic/netcdf/lfca_hadisst/Atlantic_hadisst_1900_2019_month_0_75N_local_seaice.nc"
    ; detrend EOF of ERSST, then lfc1 obtained from LFCA
    fil_S5       = "/data05/zhaobw/atlantic/netcdf/eof_ersstv5/Atlantic_ersstv5_1900_2019_month_0_75N_local_eof.nc"

    f_s1         = addfile(dir+fil_S1, "r")
    f_s2         = addfile(dir+fil_S2, "r")
    f_s3         = addfile(dir+fil_S3, "r")
    f_s4         = addfile(dir+fil_S4, "r")
    ; eof1 of ERSSTv5
    f_s5         = addfile(dir+fil_S5, "r")
    ; f_s6         = addfile(dir+fil_S6, "r")


    ssta_s1_detrend     = f_s1->LFC1
    ssta_s2_detrend     = f_s2->LFC1
    ssta_s3_detrend     = f_s3->LFC1
    ssta_s4_detrend     = f_s4->LFC1
    ssta_s5_detrend     = f_s5->PC1

    ; explained variance
    variance_ratio_erv5 = f_s1->variance_ratio(0)
    variance_ratio_cobe = f_s2->variance_ratio(0)
    variance_ratio_cobe2= f_s3->variance_ratio(0)
    variance_ratio_hadi = f_s4->variance_ratio(0)
    variance_ratio_eof = f_s5->variance_ratio(0)

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
year_start1 = 1900
year_end1   = 2019
; simple pearson correlation of ERSSTv5 with COBE/COBE2/HadISST1
corr_erv5_cobe  = escorc(amo_s1, amo_s2)
corr_erv5_cobe2 = escorc(amo_s1, amo_s3)
corr_erv5_hadi  = escorc(amo_s1, amo_s4)
corr_erv5_eof1  = escorc(amo_s1, amo_s5)

corr_hadi_cobe = escorc(amo_s4, amo_s2)
corr_hadi_cobe2 = escorc(amo_s4, amo_s3)

corr_cobe_cobe2 = escorc(amo_s2, amo_s3)
;***************创建数组保存四条序列*************************************************
lines              = new((/5,(year_end1-year_start1+1)*12/), "float")
lines(0,:)         = (/amo_s1/)
lines(1,:)         = (/amo_s2/)
lines(2,:)         = (/amo_s3/)
lines(3,:)         = (/amo_s4/)
lines(4,:)         = (/amo_s5/)
; lines(4,0:(2013-1900+1)*12-1)         = (/amo_s6/)
; lines(5,0:(2017-1900+1)*12-1)         = (/amo_s6/)
; lines(6,:)         = (/amo_s7/)
; lines(7,(1982-year_start1)*12:(2019-year_start1+1)*12-1)          = (/amo_s8/)

lines!1        = "time"
lines&time     = fspan(year_start1,year_end1,(year_end1-year_start1+1)*12)
switch = 1
if (switch.eq.1) then
  print("plot begin")
plot       = new(1, "graphic")
wks        = gsn_open_wks("eps", get_script_prefix_name())

;*****************************************
;          Plot    time series
;*****************************************
  sres                      = True      ; plot mods desired
  sres@gsnFrame             = False     ; don't advance frame
  sres@gsnDraw              = False

sres@tmXTOn                          = False  ; Setting True turns on top tick marks; False turns them off.
sres@tmYROn                          = False; Setting True turns on right tick marks; False turns them off.
sres@gsnStringFont                  = "times-roman"
   sres@tmYLLabelFont                  = "times-roman"
   sres@tmXBLabelFont                  = "times-roman" 
;add background grid line 
; sres@tmXMajorGrid                    = True
; sres@tmYMajorGrid                    = True
; sres@tmXMajorGridThicknessF          = 1.0
; sres@tmYMajorGridThicknessF          = 1.0
; sres@tmXMajorGridLineDashPattern     = 2
; sres@tmYMajorGridLineDashPattern     = 2
; sres@tmXMajorGridLineColor           = "gray40"
; sres@tmYMajorGridLineColor           = "gray40"

  sres@trXMinF              = 1900     ; set X-axis max and mins
  sres@trXMaxF              = 2020
  sres@trYMinF              = -3.      ; set Y-axis max and mins
  sres@trYMaxF              =  3.      

  sres@tmYLLabelFontHeightF = 0.015  ; font height
  sres@tmXBLabelFontHeightF = 0.015 
  sres@gsnTickMarksOn       = True        ; turn off tickmarks

  sres@gsnYRefLine          = 0.          ; Y-value for reference line

  sres@gsnCenterString      = ""       ; no center string
  sres@tiXAxisOn            = False    ; turn off X-axis title 
; this controls the size and location of the second plot
  ; sres@vpXF                 = 0.2   
  ; sres@vpYF                 = 0.7

  sres@vpWidthF             = 0.6
  sres@vpHeightF            = 0.2
  sres@tiYAxisOn            = True
  sres@tiYAxisFont          = "times-roman"
  sres@tiYAxisString        = "STD"    ; y-axis label 
  sres@tiXAxisString        = "yr"
  ;sres@tiYAxisFontHeightF   = 0.1     
  ;sres@tiMainString     = "LFC1" ; title
  sres@gsnLeftStringFontHeightF      = 0.015
  sres@gsnLeftString        = "(c) NA warming timeseries" 
  sres@gsnRightString       = ""
  sres@gsnRightStringOrthogonalPosF  = 0.01
;************画第一个序列*************************************************************************************
  polyres                  = sres
  ;xyMonoDashPattern ,if true: all the curves in the associated data object will be drawn with the 
  ;dash pattern specified by the xyDashPattern resource
  ; polyres@xyMonoDashPattern    = True         ; 
  polyres@xyDashPatterns        = (/0,0,0,0,1,1,1,1/)        ; set dash pattern
  polyres@xyLineThicknesses    = (/4,4,4,4,4,4/)      ; Define line thicknesses 
  polyres@xyLineColors         = (/"gray11","forestgreen","blue","red","gray11"/)  ; "purple2","violet","gray" "white", Define line color
  polyres@pmLegendDisplayMode    = "Always"            ; turn on legend
  polyres@pmLegendSide           = "Top"               ; Change location of 
  polyres@pmLegendParallelPosF   = .17                 ; move units right
  polyres@pmLegendOrthogonalPosF = -0.58                ; move units down
  polyres@pmLegendWidthF         = 0.05                ; Change width and
  polyres@pmLegendHeightF        = 0.09              ; height of legend.
  polyres@lgPerimOn              = False               ; turn off box around
  polyres@lgLabelFontHeightF     = .013                 ; label font height
  ; polyres@lgOrientation          = "horizontal"
  polyres@lgLabelFont            = "times-roman" 
  polyres@xyExplicitLegendLabels = (/" ERv5 ("+sprintf("%4.2f", variance_ratio_erv5)+"%)"," COBE ("+sprintf("%4.2f", variance_ratio_cobe)+"%)"," COBE2 ("+sprintf("%4.2f", variance_ratio_cobe2)+"%)"," HadI1 ("+sprintf("%4.2f", variance_ratio_hadi)+"%)"," EOF1 ("+sprintf("%4.2f", variance_ratio_eof)+"%)"/)         ; create explicit labels
  ; polyres@gsnLeftString        = "" 
  plot(0) = gsn_csm_xy(wks, lines&time, lines, polyres)
;------------------------------- add text in plot ---------------------------------------------
 txres               = True                            
 txres@txFontHeightF = 0.013             ; Set the font height
 txres@txFont        = "times-roman"
 label1 = "r(ERv5,COBE) = "+sprintf("%4.2f", corr_erv5_cobe)
 label2 = "  r(ERv5,COBE2) = "+sprintf("%4.2f", corr_erv5_cobe2)
 label3 = " r(ERv5,HadI1) = "+sprintf("%4.2f", corr_erv5_hadi)
 label4 = " r(ERv5,EOF1) = "+sprintf("%4.2f", corr_erv5_eof1)

 label5 = " r(HadI,COBE) = "+sprintf("%4.2f", corr_hadi_cobe)
 label6 = " r(HadI,COBE2) = "+sprintf("%4.2f", corr_hadi_cobe2)
 label7 = " r(COBE,COBE2) = "+sprintf("%4.2f", corr_cobe_cobe2)

 text1  = gsn_add_text(wks,plot(0),label1,1999, -1.,txres)
 text2  = gsn_add_text(wks,plot(0),label2,1999, -1.55,txres)
 text3  = gsn_add_text(wks,plot(0),label3,1998.2, -2.1,txres)
 text4  = gsn_add_text(wks,plot(0),label4,1998.2, -2.6,txres)

 text5  = gsn_add_text(wks,plot(0),label5,1970, 2.5,txres)
 text6  = gsn_add_text(wks,plot(0),label6,1970, 2.0,txres)
 text7  = gsn_add_text(wks,plot(0),label7,1970, 1.5,txres)

; - create the panel
   resp                               = True
   resp@lbOrientation                 = "Vertical"
   resp@pmLabelBarHeightF             = 0.25
   resp@pmLabelBarWidthF              = 0.05
   ; resp@lbLabelFont                   = 0
   resp@lbLabelFontHeightF            = 0.02
   resp@lbBoxEndCapStyle                     = "TriangleBothEnds"
   resp@gsnPanelMainString                 = ""
   gsn_panel(wks,plot,(/1,1/),resp)
   ; frame(wks)
   print("end of the script")

end if
end
