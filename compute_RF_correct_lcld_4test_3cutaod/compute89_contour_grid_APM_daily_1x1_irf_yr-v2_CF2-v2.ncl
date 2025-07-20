;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;*************************************************
begin

   nlat = 180 
   nlon = 360 
   setfileoption("bin","ReadByteOrder","LittleEndian")

;  cmap = (/"white", "black", \
;           "navy", "(/0.00, 0.00, 1.00/)", "(/0.00, 0.3, 1.00/)", \
;           "(/0.00, 0.5, 1.00/)", "(/0.00, 0.7, 1.00/)", "(/0.00, 0.9, 1.00/)", \
;           "limegreen", "chartreuse", "yellow", "orange", "orangered", "red", "red3", \
;           "darkorchid4" /)

   wks = gsn_open_wks("x11","compute89_contour_6dt8_APM_1x1_irf_2010-v2-CF2-v2")
;   gsn_define_colormap(wks,"BlAqGrYeOrReVi200")  ; choose colormap
  gsn_define_colormap(wks,"BlWhRe");cmap)

nmon = 12 
glv_mon = new((/nmon/),float)
var_mon = new((/nmon,nlat,nlon/),float)
do im = 0, nmon-1
;
 fils = "./2010"+sprinti("%0.2i",im+1)+"CF/compute89_APM_1by1_irf_2010"+sprinti("%0.2i",im+1)+"_meme.dat"
 print(fils)
 var = fbindirread(fils,0, (/nlat,nlon/), "float")
 var@_FillValue = -999.9
 var = -var
; print(var)
 var_mon(im,:,:) = var(:,:)
;=============================
;  create lat and long coordinate variables
;============================
  dlat  = 180./nlat
  dlon  = 360./nlon
;
  lat = ispan(-90,89,1)*dlat
  lon = ispan(-180,179,1)*dlon
;
  jlat  = nlat
  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad
  dlon2   = abs(lon(2)-lon(1))*rr
  dx     = dlon2*cos(lat*rad)
;                                     lat can have variable spacing
  dy     = new ( jlat, typeof(dx))
                                                            ; close enough
  dy(0)  = abs(lat(2)-lat(1))*rr
  dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5
  dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr

  area   = dx*dy
;print("area="+area)
;  dx     = rr*cos(lat*rad)
; print(dx)
;
  glv_area = wgt_areaave(var, dx, 1.0, 1)
;  print("GlbAve_area="+glv_area)
;;;;;
;;;;; calculate weight average
;;;;;
  res_lat = 1.0
  res_lon = 1.0
;  rlat1 = -90.
;  rlat2 =  90.
  rlat1 = -60.
  rlat2 =  60.
  rlon1 = -180.
  rlon2 =  180.
  rlat11 = rlat1 + 90.
  rlat22 = rlat2 + 89.
  rlon11 = rlon1 + 180.
  rlon22 = rlon2 + 179.
  lat1 = floattoint(rlat11/res_lat)
  lat2 = floattoint(rlat22/res_lat)
  lon1 = floattoint(rlon11/res_lon)
  lon2 = floattoint(rlon22/res_lon)

print("lat1="+lat1+"   lat2="+lat2)
print("lon1="+lon1+"   lon2="+lon2)
  ave_var = wgt_areaave(var(lat1:lat2,lon1:lon2),area(lat1:lat2), 1.0, 0)
 print("ave_IRF="+ave_var)
;=============================
;  create lat and long coordinate variables
;============================
   var!0          = "lat"
   var&lat        = ispan(-90,89,1)*dlat
   var&lat@units  = "degrees_north"

   var!1          = "lon"
   var&lon        = ispan(-180,179,1)*dlon
   var&lon@units  = "degrees_east"
;**************************************************
; create plot
;**************************************************
  res                  = True                      ; plot mods desired
  res@cnFillOn         = True                      ; turn on color
;  res@gsnSpreadColors  = True                  ; use full range of colors
  res@cnLinesOn        = False                ; no contour lines
;  res@gsnSpreadColors       = True     ; use full range of color map
;  res@gsnSpreadColorStart  = 2
;  res@gsnSpreadColorEnd    = 14

  res@gsnPaperOrientation  = "portrait"
  res@mpShapeMode  = "FreeAspect"
  res@vpWidthF      = 0.7
  res@vpHeightF     = 0.5
  res@vpXF          = 0.1
  res@vpYF          = 0.8

  res@gsnMajorLonSpacing   = 60
  res@gsnMinorLonSpacing   = 10
  res@gsnMajorLatSpacing   = 30
  res@gsnMinorLatSpacing   = 10

  res@lbOrientation        = "vertical"  ; vertical label bar
;  res@pmTickMarkDisplayMode = "Always"; use NCL default lat/lon labels
;  res@gsnAddCyclic          = False    ; data already has cyclic point
                                       ; this must also be set for any zoom
  res@mpMinLatF            = -60.      ; range to zoom in on
  res@mpMaxLatF            =  60.
  res@mpMinLonF            = -180.
  res@mpMaxLonF            =  180.

;  res@vpHeightF= 0.6                    ; change aspect ratio of plot
;  res@vpWidthF = 0.7

  res@gsnAddCyclic = True     ;
  res@gsnRightStringParallelPosF = 1.09
  res@pmLabelBarWidthF     = 0.04

  xx  = flt2string(decimalPlaces(ave_var,2,True));floattoint(100.*ave_var-0.5)/100.
  res@gsnLeftString        = "mean: " + xx
;  res@gsnRightString        = "W m~S~-2"
  res@cnFillColors    = (/2,15,25,30,35,40,43,46,51,54,57,62,69,80,90/)
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels            = (/-15,-10,-5,-2.5,-1.5,-1.,-0.5,-0.1,0.1,0.5,1,1.5,2.5,5/)
;  res@gsnCenterString     = "Satellite-based estimates of IRF"                     ; title

  res@tiMainString="2010"+sprinti("%0.2i",im+1)

;  plot = gsn_csm_contour_map_ce(wks,var, res)

glv_mon(im) = ave_var
;
end do
;
var_ann = dim_avg_n_Wrap(var_mon,0)
ave_ann = wgt_areaave(var_ann(lat1:lat2,lon1:lon2),area(lat1:lat2), 1.0, 0)
 print("ave_ann_IRF="+ave_ann)
;=============================
;  create lat and long coordinate variables
;============================
   var_ann!0          = "lat"
   var_ann&lat        = ispan(-90,89,1)*dlat
   var_ann&lat@units  = "degrees_north"

   var_ann!1          = "lon"
   var_ann&lon        = ispan(-180,179,1)*dlon
   var_ann&lon@units  = "degrees_east"

  xx  = flt2string(decimalPlaces(ave_ann,2,True));floattoint(100.*ave_ann-0.5)/100.
  res@gsnLeftString        = "mean: " + xx
;  res@gsnRightString        = "W m~S~-2"
  res@lbBoxEndCapStyle     = "TriangleBothEnds"
  res@lbTitleString      ="RF~B~aci~N~ (W m~S~-2~N~)"
  res@lbTitleFontHeightF = 0.02

  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels            = (/-15,-10,-5,-2.5,-1.5,-1.,-0.5,-0.1,0.1,0.5,1,1.5,2.5,5/)
;  res@gsnCenterString     = "Satellite-based estimates of IRF"                     ; title

  plot = gsn_csm_contour_map_ce(wks,var_ann, res)
;
lines=sprintf("%10.4f",glv_mon)
asciiwrite("compute89_glv_irf_6dt8_2010-v2-CF.txt",lines)
end
