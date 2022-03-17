import numpy as np
import re
import Ngl
from sys import argv

path = argv[1]
schemes = ["Ah21","Ah42","Ah43","Ah63"]
t1 = 401
t2 = 1400
Nc = 96
wktype = "pdf"

def plot_Eldred(scheme, N, path,t1,t2):

    Nlon = 4*N
    Nlat = 2*N+1
    suffix = "_N{:03d}_".format(N)+scheme
    
    cn_res = Ngl.Resources()
    cn_res.cnFillOn = True
    cn_res.cnLinesOn = False
    #cn_res.cnLineThicknessF = 0.5
    cn_res.cnLineLabelsOn = False
    cn_res.mpCenterLonF = 180.0
    cn_res.mpGridAndLimbOn = False
    cn_res.lbLabelBarOn = True
    cn_res.lbOrientation="Horizontal"
    cn_res.lbLabelFontHeightF = 0.01
    cn_res.tiMainFontHeightF = 0.01
    cn_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    cn_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
   
   
    wkres = Ngl.Resources()
    wkres.nglMaximize = True
    wkres.wkColorMap = "cmp_b2r"
    wks = Ngl.open_wks(wktype, "Eldred_curl_N{:03d}_".format(N)+scheme,wkres)

    fd = open(path+"/curl"+suffix+".dat","rb")
    fd.seek(4*t2*Nlon*Nlat,0)
    z = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape((Nlat,Nlon))
    fd.close()
    plot = Ngl.contour_map(wks,z*1e5,cn_res)
    Ngl.delete_wks(wks)

    cn_res.nglDraw = False
    cn_res.nglFrame = False

    def get_mean_and_std(fname, t1, t2):
        Nt = t2-t1+1.0
        f_mean  = np.zeros((Nlat,Nlon),dtype=np.float64)
        f_mean2 = np.zeros((Nlat,Nlon),dtype=np.float64)
        fd = open(fname,"rb")
        fd.seek(4*t1*Nlon*Nlat,0)
        #lon=174
        #lat = 63
        for i in range(t2-t1+1):
            f = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape((Nlat,Nlon))
            f_mean  = f_mean+f/Nt
            f_mean2 = f_mean2+f**2/Nt
            #print "f", f[lat,lon], Nt


        fd.close()
        return f_mean, np.sqrt(f_mean2-f_mean**2)

    z_mean, z_std     = get_mean_and_std(path+"/curl"+suffix+".dat", t1, t2)
    div_mean, div_std = get_mean_and_std(path+"/div"+suffix+".dat", t1, t2)
    h_mean, h_std     = get_mean_and_std(path+"/h"+suffix+".dat", t1, t2)

    wks = Ngl.open_wks(wktype, "Eldred_N{:03d}_".format(N)+scheme,wkres)
    plot1 = Ngl.contour_map(wks,z_mean*1e5,cn_res)
    plot2 = Ngl.contour_map(wks,z_std*1e5,cn_res)
    plot3 = Ngl.contour_map(wks,div_mean*1e8,cn_res)
    plot4 = Ngl.contour_map(wks,div_std*1e8,cn_res)
    plot5 = Ngl.contour_map(wks,h_mean,cn_res)
    plot6 = Ngl.contour_map(wks,h_std,cn_res)
    pres = Ngl.Resources()
    Ngl.panel(wks,(plot5,plot6,plot1,plot2,plot3,plot4),(3,2),pres)
    Ngl.delete_wks(wks)

    wks = Ngl.open_wks(wktype, "xy_N{:03d}_".format(N)+scheme,wkres)
    xy_res = Ngl.Resources()
    xy_res.nglDraw = False
    xy_res.nglFrame = False
    plot1 = Ngl.xy(wks,cn_res.sfYArray,np.mean(h_std,axis=1),xy_res)
    plot2 = Ngl.xy(wks,cn_res.sfYArray,np.mean(z_std*1e5,axis=1),xy_res)
    Ngl.panel(wks,(plot1,plot2),(1,2),pres)
    Ngl.delete_wks(wks)

for scheme in schemes:
    plot_Eldred(scheme,Nc,path,t1,t2)

