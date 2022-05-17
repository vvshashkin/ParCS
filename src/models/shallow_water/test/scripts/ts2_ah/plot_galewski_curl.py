import numpy as np
import re
import Ngl
from sys import argv

path = argv[1]
wktype = "eps"
schemes = ["Ah21","Ah42","Ah43","Ah63"]
cn_res = Ngl.Resources()
cn_res.cnFillOn = True
cn_res.cnLinesOn = True
cn_res.cnLineThicknessF = 0.5
#cn_res.cnLineColor = "darkgray"
cn_res.mpGridAndLimbOn = True
cn_res.mpGridLineThicknessF=0.3
cn_res.mpGridLineDashPattern = 1
cn_res.cnLineLabelsOn = False
cn_res.mpCenterLonF = 180.0
cn_res.lbLabelBarOn = False
cn_res.lbOrientation="Horizontal"
cn_res.lbLabelFontHeightF = 0.01
cn_res.tiMainFontHeightF = 0.01
cn_res.mpLimitMode = "LatLon"
cn_res.mpMaxLatF = 90.0
cn_res.mpMinLatF = 20.0
cn_res.cnLevelSelectionMode = "ExplicitLevels"
cn_res.cnLevels = [-12,-10,-8,-6,-4,-2,-0.5,0.5,2,4,6,8,10,12]
cn_res.cnFillColors = [2,4,6,7,8,9,10,0,12,13,14,15,17,18,19]

cn_res.nglDraw = False
cn_res.nglFrame = False

wkres = Ngl.Resources()
wkres.nglMaximize = True
wkres.wkColorMap = "BlueDarkRed18"

scheme_conv = "Ah42"
wks_conv = Ngl.open_wks(wktype, "curl_"+scheme_conv,wkres)
plots_conv = []

#for N in [32, 64, 128, 256]:#[20,40,80,160]:
for N in [96, 128, 192, 256]:#[20,40,80,160]:
    print "N=",N
    Nlon = 4*N
    Nlat = 2*N+1
    cn_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    cn_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
    wks = Ngl.open_wks(wktype, "curl_N{:03d}".format(N),wkres)

    plots = []
    for scheme in schemes:
        print scheme
        fd = open("curl_N"+"{:03d}_".format(N)+scheme+".dat","rb")
        fd.seek(4*Nlon*Nlat*6,0)
        z = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape(Nlat,Nlon)
        fd.close()
        cn_res.tiMainString = scheme
        print "plot"
        plots.append(Ngl.contour_map(wks,z*1e5,cn_res))
        if scheme == scheme_conv:
            cn_res.tiMainString = scheme+" N~B~c~N~="+str(N)
            plots_conv.append(Ngl.contour_map(wks_conv,z*1e5,cn_res))

    pres = Ngl.Resources()
    pres.nglPanelLabelBar = True
    Ngl.panel(wks,plots,(4,1),pres)
    Ngl.delete_wks(wks)
Ngl.panel(wks_conv,plots_conv,(4,1),pres)
Ngl.delete_wks(wks_conv)
