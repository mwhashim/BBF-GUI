from __future__ import division
import os, sys
sys.setrecursionlimit(50000) # to solve maximum recursion depth exceeded error !!

import shutil
import copy

import logging
import glob
from threading import Thread
import Queue
#from KThread import *

import time

from numpy import *
from scipy import ndimage

import decimal
from collections import *

#import paramiko
#paramiko.util.log_to_file('/tmp/paramiko.log')

#----------------------------------
import matplotlib
matplotlib.use('TkAgg')

#from six.moves import tkinter_filedialog as FileDialog

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.figure import Figure
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from pylab import *
from itertools import cycle
import pandas as pd

#from six.moves import tkinter_filedialog as FileDialog

#from scipy.interpolate import interp1d, UnivariateSpline, InterpolatedUnivariateSpline

#import VerticalScrolledFrame as VSF
#import Tooltip
#import ViewLog as VL

#import ParamsDict as PrDt
#import MakefileDict as MkDt
#import TooltipDict as TtDt

import subprocess, shlex

import ttk
import tkFileDialog, Tkconstants
#import tkinter.filedialog, tkinter.messagebox
#
if sys.version_info[0] < 3:
    from Tkinter import *
    from Tkinter import _setit
else:
    from tkinter import *
    from tkinter import _setit


#from mttkinter import *


import cv2
from PIL import Image, ImageTk
#from CosmoSuite import *

from astropy.cosmology import wCDM
from astropy.io import fits

from emailling import *

from textdictENG import text_dict

#----------------------------------
def destroy(e): sys.exit()

#--- CAM Read ------
width, height = 480, 480; dpi = 100
cap = cv2.VideoCapture(0)
cap.set(cv2.CAP_PROP_FRAME_WIDTH, width)
cap.set(cv2.CAP_PROP_FRAME_HEIGHT, height)

#---- Alignment-----
# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height

sigmaval = 1.; truncateval = 3.

#--- CWD -----
CWD = os.getcwd()#; print CWD
#CWD = '/Users/mahmoud/Desktop/BBFpipeline_gui' #os.path.dirname(os.path.abspath(__file__))
#-----------------------------
def readimage(fileimage):
    image = Image.open(fileimage)
    xsize, ysize = image.size
    return image, xsize, ysize

def buildlens(filelens):
    image_file1 = filelens + "_alpha1.fits"
    hdu_list1 = fits.open(image_file1)
    image_data1 = fits.getdata(image_file1)
    image_file2 = filelens + "_alpha2.fits"
    image_data2 = fits.getdata(image_file2)
    return image_data1, image_data2


def deflect(image_arr,image_data1,image_data2,xsize,ysize, scalefac):
#    
#    for ll in range(0,len(zz)):\n",
#        zs=zz[ll]\n",
#        dl = np.array(cosmo.angular_diameter_distance(zl)*cosmo.H(0.)/100.)
#        ds = np.array(cosmo.angular_diameter_distance(zs)*cosmo.H(0.)/100.)
#        dls = (np.array(cosmo.angular_diameter_distance(zs)*cosmo.H(0.)/100.)*(1+zs) - np.array(cosmo.angular_diameter_distance(zl)*cosmo.H(0.)/100.)*(1+zl))/(1+zs)
    ds = 1156.34008206061
    dl = 564.558513269509
    dls = 803.491011267171 + scalefac
    
    f = ds/dl/dls/xsize*2.5e8
    
    lensed_data = copy(image_arr)
    #print xsize, ysize
    for i in range(0,xsize):
        for j in range(0,ysize):
            ii = i - image_data2[i][j]*f + 0.5
            if int(ii) >= xsize:
                ii = int(ii - xsize)
            if int(ii) < -0.5:
                ii = int(ii + xsize)
            jj = j - image_data1[i][j]*f + 0.5
            if int(jj) >= ysize:
                jj = int(jj - ysize)
            if int(jj) < -0.5:
                jj = int(jj + ysize)
            lensed_data[i][j] = image_arr[int(ii),int(jj)]

    return lensed_data

def center(toplevel):
    toplevel.update_idletasks()
    w = toplevel.winfo_screenwidth()
    h = toplevel.winfo_screenheight()
    size = tuple(int(_) for _ in toplevel.geometry().split('+')[0].split('x'))
    x = w/2 - size[0]/2
    y = h/2 - size[1]/2
    toplevel.geometry("%dx%d+%d+%d" % (size + (x, y)))

class AnchoredHScaleBar(matplotlib.offsetbox.AnchoredOffsetbox):
    """ size: length of bar in data units
        extent : height of bar ends in axes units """
    def __init__(self, size=1, extent = 0.02, label="", loc=2, ax=None, pad=0.4, borderpad=0.5, ppad = 0, sep=2, prop=None, frameon=True, **kwargs):
        if not ax:
            ax = plt.gca()
        trans = ax.get_xaxis_transform()
        size_bar = matplotlib.offsetbox.AuxTransformBox(trans)
        line = Line2D([0,size],[0,0], **kwargs)
        vline1 = Line2D([0,0],[-extent/2.,extent/2.], **kwargs)
        vline2 = Line2D([size,size],[-extent/2.,extent/2.], **kwargs)
        size_bar.add_artist(line)
        size_bar.add_artist(vline1)
        size_bar.add_artist(vline2)
        txt = matplotlib.offsetbox.TextArea(label, minimumdescent=False, textprops=dict(color="white", size=10))
        self.vpac = matplotlib.offsetbox.VPacker(children=[size_bar,txt], align="right", pad=ppad, sep=sep)
        matplotlib.offsetbox.AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad, child=self.vpac, prop=prop, frameon=frameon)

#----------------------------------
class Application(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.grid()
        self.master.title(text_dict['t1'])#("Big Bang Factory")

        #---------------------------------------
        for r in range(7):
            self.master.rowconfigure(r, weight=2)
#        for c in range(5):
#            self.master.columnconfigure(c, weight=2)
        #---------------------------------------
        self.master.columnconfigure(3, weight=2)
        self.master.columnconfigure(4, weight=2)
        self.master.columnconfigure(5, weight=2)
        #self.master.rowconfigure(0, weight=0)
        #self.master.columnconfigure(0, weight=0)
        self.welcome_screen(self.master)
            
    def goback(self):
        self.status.grid_forget(); self.Frame_2.grid_forget(); self.Frame_3.grid_forget(); self.welcome_screen(self.master)
    
    def welcome_screen(self, frame):
        # Frame 0 :------------------------------------------------------
        self.Frame_0 = Frame(frame, bg="white smoke")
        self.Frame_0.grid(row = 0, column = 0, rowspan = 7, columnspan = 6, sticky = W+E+N+S)
        
        self.original = Image.open("BBF-Logo.gif"); self.Background_photo = ImageTk.PhotoImage(self.original)
        self.Welcome_Frame = Canvas(self.Frame_0) #, image = Background_photo
        self.Welcome_Frame.create_image(0, 0, image=self.Background_photo, anchor=NW, tags="IMG")
        self.Welcome_Frame.pack(side=TOP, fill=BOTH, expand=1)
        self.Frame_0.bind("<Configure>", self.resize)
        
        Welcome_buttom = Button(self.Frame_0, text=text_dict['t3'], command=self.Sim_Create)
        Welcome_buttom.pack()
        
        #---------------------------------------
    def Sim_Create(self):
        self.Frame_0.destroy()
        self.initialize()
        self.Menubar()
        #self.Toolbar(self.Frame_1);
        self.fig, self.ax, self.canvas = self.PlotPan(self.Frame_3)
        self.MainFrame(self.Frame_2)

    def initialize(self):

        #--------------------------------------
        style = ttk.Style()
        style.layout('TNotebook.Tab', [])

        # Status Bar :-----
        self.statusVariable = StringVar()
        self.status = Label(self.master, textvariable = self.statusVariable, anchor = "w", fg = "yellow", bg = "blue")
        self.status.grid(column = 0, row = 7,columnspan = 7, sticky = 'EWS')
        self.statusVariable.set(text_dict['t2'])#(u"Welcome to Big Bang Factory")

#        # Frame 1 :------------------------------------------------------
#        self.Frame_1 = PanedWindow(self.master, bg="white smoke")
#        self.Frame_1.grid(row = 0, column = 0, rowspan = 1, columnspan = 5, sticky = W+E+N+S)

        # Frame 2 :-----------------------------------------------------
        self.Frame_2 = Frame(self.master, bg="white", bd= 1, relief= RIDGE)
        self.Frame_2.grid(row = 0, column = 0, rowspan = 7, columnspan = 2, sticky = W+E+N+S)

        # Frame 3 :-----------------------------------------------------
        self.Frame_3 = Frame(self.master, bg="white", bd= 3, relief= GROOVE)
        self.Frame_3.grid(row = 0, column = 3, rowspan = 7, columnspan = 3, sticky = W+E+N+S)
#
#        # Frame 4 :------------------------------------------------------
#        self.Frame_4 = Frame(self.master, bg="white smoke", bd= 5, relief= RIDGE)
#        self.Frame_4.grid(row = 4, column = 2, rowspan = 2, columnspan = 3, sticky = W+E+N+S)
    def Menubar(self):
        self.menubar = Menu(self.master)
        
        menu1 = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label=text_dict['t4'], menu=menu1)
        menu1.add_command(label=text_dict['t5'], command=self.run_reset)
        menu1.add_command(label=text_dict['t6'], command=self.save_movie)
        menu1.add_command(label=text_dict['t7'], command=self.send_movie)
        #menu1.add_command(label=text_dict['t8'], command=self.goback)
        menu1.add_separator()
        menu1.add_command(label=text_dict['t9'], command=root.quit)
        
        menu2 = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label=text_dict['t10'], menu=menu2)
        menu2.add_command(label=text_dict['t11'], command=self.ModelDirectory)
        
        menu3 = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label=text_dict['t12'], menu=menu3)
        menu3.add_command(label=text_dict['t13'], command=self.myDialog)
        
        self.master.config(menu=self.menubar)

    def Toolbar(self, frame):
        Quit_photo = PhotoImage(file="quit.gif")
        Quit = Button(frame, text=u"Quit", image = Quit_photo, command=sys.exit)
        Quit.photo = Quit_photo
        Quit.grid(column=6, row=0, sticky= W+E+N+S, pady = 5)
        Quit.pack(side="right")

        run_photo = PhotoImage(file="run.gif")
        self.run = Button(frame, text = u"Home", fg='blue', image = run_photo, command = self.callback_RunSpicf)
        self.run.photo = run_photo
        self.run.grid(column = 1, row = 0, sticky = W+E+N+S, pady = 5)
        self.run.pack(side="left")

    def PlotPan(self, frame):
        fig = plt.Figure()
        ax = fig.add_subplot(111); ax.axis('off'); ax.get_xaxis().set_visible(False); ax.get_yaxis().set_visible(False)
        fig.set_tight_layout(True) # fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master = frame)
        #self.toolbar = NavigationToolbar2TkAgg(self.canvas, frame)
        #canvas.get_tk_widget().grid(column = 0, row = 0, pady = 5)
        #canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        return fig, ax, canvas
    
    def MainFrame(self, frame):

        #------------------------
        self.main_nb = ttk.Notebook(frame, padding = -5)
        self.Main_page = ttk.Frame(self.main_nb)#; self.Anals_page = ttk.Frame(self.main_nb)

        #------------------------
        self.main_nb.add(self.Main_page, text='Main')
        #self.main_nb.add(self.Anals_page, text='Analysis & Results')
        #self.main_nb.add(self.SysPerf_page, text='System Preferences')
        self.main_nb.pack(side="top", fill="both", expand=True)

        self.main_nb.select(self.Main_page)
        #------------------------
        self.MainPage(self.Main_page)#; self.AnalysResultsPage(self.Anals_page)

    def MainPage(self, page):

#        for r in range(5):
#            page.rowconfigure(r, weight=2)
#
#        for x in range(2):
#            page.columnconfigure(x, weight=2)

        #------------------------
        self.Entry_Frame = Frame(page, bg="white", bd= 1, relief= RIDGE)
        self.Entry_Frame.grid(row = 0, column = 0, rowspan = 4, sticky = W+E+N+S)

        self.NewRun_nb = ttk.Notebook(self.Entry_Frame,  padding = -5)
        self.NewRun_page = ttk.Frame(self.NewRun_nb)#; self.Codes_Compile_page = ttk.Frame(self.NewRun_nb); self.Run_Specif_page = ttk.Frame(self.NewRun_nb)

        self.NewRun_nb.add(self.NewRun_page, text='User Details')
        #self.NewRun_nb.add(self.Codes_Compile_page, text='Codes Compilation')
        #self.NewRun_nb.add(self.Run_Specif_page, text='Run Specifics')

        self.NewRun_nb.pack(side="top", fill="both", expand=True)

        #--------------------------------------
        for x in range(13):
            Grid.columnconfigure(self.NewRun_page, x, weight=2)

        for y in range(6):
            Grid.rowconfigure(self.NewRun_page, y, weight=2)

        self.Name_Dict_group = LabelFrame(self.NewRun_page, text = text_dict['t14'])
        self.Name_Dict_group.grid(row = 0, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 0, weight=0)

        for x in range(2):
            Grid.columnconfigure(self.Name_Dict_group, x, weight=2)
        for y in range(2):
            Grid.rowconfigure(self.Name_Dict_group, y, weight=2)

        Label(self.Name_Dict_group, text=text_dict['t15']).grid(row=0, column=0, sticky= W)
        self.Name_Var = StringVar()
        self.User_Name = Entry(self.Name_Dict_group, textvariable=self.Name_Var)
        self.User_Name.grid(row=0, column=1, sticky= W+E+N+S,  columnspan = 6)
        #self.name_Var.set(None)

        Label(self.Name_Dict_group, text=text_dict['t16']).grid(row=1, column=0, sticky= W)
        self.Email_Var = StringVar()
        self.User_Email = Entry(self.Name_Dict_group, textvariable=self.Email_Var)
        self.User_Email.grid(row=1, column=1, sticky= W+E+N+S,  columnspan = 6)
        #self.Email_Var.set(None)

        self.Image_button = Button(self.Name_Dict_group, text=text_dict['t17'], command = self.open_cam)
        self.Image_button.grid(row = 2, column = 2, columnspan = 1 , sticky = W+E+N+S, pady = 5)

        self.click_button = Button(self.Name_Dict_group, text=text_dict['t18'], command = self.saveImage)
        self.click_button.grid(row = 2, column = 3, columnspan = 1 , sticky = W+E+N+S, pady = 5)

        self.Cosmo_Parms_group = LabelFrame(self.NewRun_page, text = text_dict['t19'])
        self.Cosmo_Parms_group.grid(row = 1, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.Cosmo_Parms_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.Cosmo_Parms_group, y, weight=2)

        Omega_m_photo = PhotoImage(file="Omega_m.gif")
        Omega_Lambda_photo = PhotoImage(file="Omega_lambda.gif")


        Omega_Lambda_label = Label(self.Cosmo_Parms_group, text="Omega_Lambda", image = Omega_Lambda_photo)
        Omega_Lambda_label.photo = Omega_Lambda_photo
        Omega_Lambda_label.grid(row=0, column=0, sticky= W+E+N+S, pady = 5)

        self.Omega_l_Var= DoubleVar()
        #self.Omega_l_Var.trace("w", self.callback_NBodyTrace)
        self.Omega_l=Scale(self.Cosmo_Parms_group, from_=0.0, to=1.0, resolution=0.25, orient=HORIZONTAL, width=15, length=400, variable = self.Omega_l_Var, digits=3)
        #Entry(self.Cosmo_Parms_group, textvariable=self.Omega_l_Var, width=10)
        self.Omega_l.grid(row=0, column=1, sticky= W+E+N+S, pady = 5, columnspan = 3)


        Omega_m_label = Label(self.Cosmo_Parms_group, text="Omega_m", image = Omega_m_photo)
        Omega_m_label.photo = Omega_m_photo
        Omega_m_label.grid(row=1, column=0, sticky= W+E+N+S, pady = 5)

        self.Omega_m_Var= DoubleVar()
        #self.Omega_m_Var.trace("w", self.callback_NBodyTrace)
        self.Omega_m= Scale(self.Cosmo_Parms_group, from_=0.0, to=1.0, resolution=0.25, orient=HORIZONTAL, width=15, length=400,variable = self.Omega_m_Var, digits=3)
        #Entry(self.Cosmo_Parms_group, textvariable=self.Omega_m_Var, width=10)
        self.Omega_m.grid(row=1, column=1, sticky= W+E+N+S, pady = 5, columnspan = 3)


        self.DarkEnergy_group = LabelFrame(self.NewRun_page, text = text_dict['t20'])
        self.DarkEnergy_group.grid(row = 2, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.DarkEnergy_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.DarkEnergy_group, y, weight=2)


        self.Lambda_Var = StringVar()
        self.Lambda_Var.trace('w', self.models_refresh)
        self.Lambda_RadBtt = Radiobutton(self.DarkEnergy_group, text = text_dict['t21'], variable = self.Lambda_Var, value = 'Lambda_')
        self.Lambda_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.Lambda_RadBtt = Radiobutton(self.DarkEnergy_group, text = text_dict['t22'], variable = self.Lambda_Var, value = 'Quint_')
        self.Lambda_RadBtt.grid(row = 0, column = 1, sticky = W)

        self.Lambda_RadBtt = Radiobutton(self.DarkEnergy_group, text = text_dict['t23'], variable = self.Lambda_Var, value = 'Phantom_')
        self.Lambda_RadBtt.grid(row = 0, column = 2, sticky = W)
        self.Lambda_Var.set('Lambda_')

        self.DarkMatter_group = LabelFrame(self.NewRun_page, text = text_dict['t24'])
        self.DarkMatter_group.grid(row = 3, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.DarkMatter_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.DarkMatter_group, y, weight=2)

        self.CDM_Var = StringVar()
        self.CDM_Var.trace('w', self.models_refresh1)
        self.CDM_RadBtt = Radiobutton(self.DarkMatter_group, text = text_dict['t25'], variable = self.CDM_Var, value = 'Lambda_')
        self.CDM_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.CDM_RadBtt = Radiobutton(self.DarkMatter_group, text = text_dict['t26'], variable = self.CDM_Var, value = 'wDM_0.5-')
        self.CDM_RadBtt.grid(row = 0, column = 1, sticky = W)
        self.CDM_Var.set('Lambda_')

        self.IniMatter_group = LabelFrame(self.NewRun_page, text = text_dict['t27'])
        self.IniMatter_group.grid(row = 4, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.IniMatter_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.IniMatter_group, y, weight=2)

        self.IniM_Var = StringVar()
        self.IniM_Var.trace('w', self.models_refresh2)
        self.IniM_RadBtt = Radiobutton(self.IniMatter_group, text = text_dict['t28'], variable = self.IniM_Var, value = 'Lambda_')
        self.IniM_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.IniM_RadBtt = Radiobutton(self.IniMatter_group, text = text_dict['t29'], variable = self.IniM_Var, value = 'LocalPNG_1000-')
        self.IniM_RadBtt.grid(row = 0, column = 1, sticky = W)

        self.IniM_RadBtt = Radiobutton(self.IniMatter_group, text = text_dict['t30'], variable = self.IniM_Var, value = 'LocalPNG_-1000-')
        self.IniM_RadBtt.grid(row = 0, column = 2, sticky = W)
        self.IniM_Var.set('Lambda_')

        self.MG_group = LabelFrame(self.NewRun_page, text = text_dict['t31'])
        self.MG_group.grid(row = 5, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.MG_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.MG_group, y, weight=2)

        self.MG_Var = StringVar()
        self.MG_Var.trace('w', self.models_refresh3)
        self.MG_RadBtt = Radiobutton(self.MG_group, text = text_dict['t32'], variable = self.MG_Var, value = 'Lambda_')
        self.MG_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.MG_RadBtt = Radiobutton(self.MG_group, text = text_dict['t33'], variable = self.MG_Var, value = 'MGfR_-1e-04-')
        self.MG_RadBtt.grid(row = 0, column = 1, sticky = W)
        self.MG_Var.set('Lambda_')


        #------------------------
        self.Control_Frame = Frame(page, bg="white smoke", bd= 1, relief= RIDGE)
        self.Control_Frame.grid(row = 4, column = 0, rowspan = 1, sticky = W+E+N+S)

        for x in range(4):
            Grid.columnconfigure(self.Control_Frame, x, weight=2)

        for y in range(2):
            Grid.rowconfigure(self.Control_Frame, y, weight=2)
        
        self.Simulation_group = LabelFrame(self.Control_Frame, text = text_dict['t34'])
        self.Simulation_group.grid(row = 0, column = 0, columnspan = 4, sticky = W+E+N+S)

        #Grid.rowconfigure(self.NewRun_page, 0, weight=0)

        for x in range(6):
            Grid.columnconfigure(self.Simulation_group, x, weight=2)
        for y in range(2):
            Grid.rowconfigure(self.Simulation_group, y, weight=2)

        Grid.rowconfigure(self.Simulation_group, 0, weight=0)
        self.Simulation_Run = Button(self.Simulation_group, text = text_dict['t35'], foreground = 'red', command = self.start)
        self.Simulation_Run.grid(column = 0, row = 0, pady = 5, sticky= W+E+N+S)

        self.progress_var = DoubleVar()
        self.progress = ttk.Progressbar(self.Simulation_group, variable=self.progress_var,  orient="horizontal", length=300, mode="determinate")
        self.progress.grid(column = 0, row = 1, pady = 5, sticky= W+E+N+S, columnspan = 6)


        self.Lensing_group = LabelFrame(self.Control_Frame, text = text_dict['t36'])
        self.Lensing_group.grid(row = 1, column = 0, columnspan = 4, sticky = W+E+N+S)
        
        #Grid.rowconfigure(self.NewRun_page, 0, weight=0)
        
        for x in range(4):
            Grid.columnconfigure(self.Lensing_group, x, weight=2)
        for y in range(5):
            Grid.rowconfigure(self.Lensing_group, y, weight=2)

        self.Sow_Lensing_Map = Button(self.Lensing_group, text=text_dict['t37'], command = self.showlensMap, width = 15)
        self.Sow_Lensing_Map.grid(row=0, column=0, columnspan = 2, sticky= W)
            
        self.Sow_Lensing_User = Button(self.Lensing_group, text=text_dict['t38'], command = self.MapLensedImage, width = 15)
        self.Sow_Lensing_User.grid(row=0, column=2, columnspan = 2, sticky= W)

        self.Halo_Lensing_Map = Button(self.Lensing_group, text=text_dict['t39'], command = self.showlenscluster, width = 15)
        self.Halo_Lensing_Map.grid(row=1, column=0, columnspan = 2, sticky= W)
        
        self.Halo_Lensing_User = Button(self.Lensing_group, text=text_dict['t40'], command = self.HaloLensedImage, width = 15)
        self.Halo_Lensing_User.grid(row=1, column=2, columnspan = 2, sticky= W)
        
        Label(self.Lensing_group, text=text_dict['t41'], justify=LEFT, anchor=W).grid(row=2, column=0, sticky= W+E+N+S, pady = 5)
            
        self.ComvDist_Var= DoubleVar()
        #self.ComvDist_Var.trace("w", self.callback_NBodyTrace)
        self.ComvDist=Scale(self.Lensing_group, from_=0.0, to=100.0, resolution=10.0, orient=HORIZONTAL, width=15, length=300, variable = self.ComvDist_Var, digits=2)
        self.ComvDist.grid(row=3, column=0, sticky= W+E+N+S, pady = 5, columnspan = 4)


        self.Logo_group = LabelFrame(self.Control_Frame)
        self.Logo_group.grid(row = 3, column = 0, columnspan = 4, sticky = W+E+N+S)

        logoimg = ImageTk.PhotoImage(file="SIMCODE.png")
        LogoButton = Button(self.Logo_group, command = self.goback)
        LogoButton.config(image=logoimg)
        LogoButton.image = logoimg
        LogoButton.pack(fill=BOTH, expand=1)

    def MapLensedImage(self):
        self.ax.clear(); self.ax.axis('off')
        fileimage = CWD + "/tmp/" + self.img_filename + "_Photo.jpg"
        
        Simu_Dir = self.model_select()+"/Lens-Maps/"
        filelens = self.simdir + "/" + Simu_Dir + self.model_name +'kappaBApp_2.fits'
        
        image, xsize, ysize = readimage(fileimage); image_arr = np.array(image)
        
        alpha1, alpha2 = buildlens(filelens)
        self.maplensedimage = deflect(image_arr, alpha1, alpha2, xsize, ysize, self.ComvDist_Var.get()); self.ax.imshow(self.maplensedimage)
        self.ax.axis('off'); self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False); self.canvas.show()
        imsave(self.savedir + "/" + self.img_filename + "_LensedMap_Photo.jpg", self.maplensedimage)
    
    def HaloLensedImage(self):
        self.ax.clear(); self.ax.axis('off')
        fileimage = CWD + "/tmp/" + self.img_filename + "_Photo.jpg"
        
        Simu_Dir = self.model_select()+"/Lens-Maps/"
        filelens = self.simdir + "/" + Simu_Dir + self.model_name +'_Halo.fits'
        
        image, xsize, ysize = readimage(fileimage); image_arr = np.array(image)
        
        alpha1, alpha2 = buildlens(filelens)
        self.halolensedimage = deflect(image_arr, alpha1, alpha2, xsize, ysize, self.ComvDist_Var.get()); self.ax.imshow(self.halolensedimage)
        self.ax.axis('off'); self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False); self.canvas.show()
        imsave(self.savedir + "/" + self.img_filename + "_LensedHalo_Photo.jpg", self.halolensedimage)


    def showlensMap(self):
        self.ax.clear(); self.ax.axis('off')
        Simu_Dir = self.model_select()+"/Lens-Maps/"
        filename = self.simdir + "/" + Simu_Dir + self.model_name +'kappaBApp_2.fits'; self.Lens_map = fits.getdata(filename, ext=0)
        LenImg = self.ax.imshow(self.Lens_map + 1, cmap=matplotlib.cm.magma, norm=matplotlib.colors.LogNorm(), interpolation="bicubic") #vmin=1., vmax=1800., clip = True
        self.ax.axis('off'); self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False); self.canvas.show()
        imsave(self.savedir + "/" + self.img_filename + "_LensedMap.jpg", self.Lens_map)
    
    def showlenscluster(self):
        self.ax.clear(); self.ax.axis('off')
        Simu_Dir = self.model_select()+"/Lens-Maps/"
        filename = self.simdir + "/" + Simu_Dir + self.model_name +'_Halo.fits'; self.Halo_map = fits.getdata(filename, ext=0)
        HaloImg = self.ax.imshow(self.Halo_map + 1, cmap=matplotlib.cm.magma, norm=matplotlib.colors.LogNorm(), interpolation="bicubic") #vmin=1., vmax=1800., clip = True
        self.ax.axis('off'); self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False); self.canvas.show()
        imsave(self.savedir + "/" + self.img_filename + "_LensedHalo.jpg", self.Halo_map)
    
    def run_reset(self):
        self.progress_var.set(0); self.ax.clear(); self.ax.axis('off'); self.canvas.show()
        self.Name_Var.set(''); self.Email_Var.set(''); self.Omega_m_Var.set(0.0); self.Omega_l_Var.set(0.0)
        self.Lambda_Var.set('Lambda_'); self.CDM_Var.set('Lambda_'); self.IniM_Var.set('Lambda_'); self.MG_Var.set('Lambda_')
        shutil.rmtree(CWD + "/tmp/")
        
        #if self._job1 is not None: self.after_cancel(self._job1); self._job1 = None
    

    def model_select(self):
        if self.Lambda_Var.get() != 'Lambda_':
            run_type = self.Lambda_Var.get()
            if run_type == "Quint_":
                self.wx = -0.9
            else:
                self.wx = -1.1

        elif self.CDM_Var.get()  != 'Lambda_':
            run_type = self.CDM_Var.get(); self.wx = -1.0
        elif self.IniM_Var.get() != 'Lambda_':
            run_type = self.IniM_Var.get(); self.wx = -1.0
        elif self.MG_Var.get() != 'Lambda_':
            run_type = self.MG_Var.get(); self.wx = -1.0
        else:
            run_type = 'Lambda_'; self.wx = -1.0

        if self.Omega_m_Var.get() == 0.0:
            Omega_m = 0.1; Omega_m = str(Omega_m)
        else:
            Omega_m = str(self.Omega_m_Var.get())

        model  = "BBF_" + run_type + Omega_m + "-" + str(self.Omega_l_Var.get())
        self.model_name = run_type + Omega_m + "-" + str(self.Omega_l_Var.get())
        print model

        return model

    def save_movie(self):
        self.SaveDirectory()
        writer = animation.writers['ffmpeg'](fps=15)
        self.ani.save(self.savedir + "/" + self.img_filename + "_movie.mp4", writer=writer, dpi=dpi)
        video_file = self.savedir + "/" + self.img_filename + "_movie.mp4"
        muxvideo_file = self.savedir + "/" + self.img_filename + "mux_movie.mp4"
        audio_file = "ChillingMusic.wav"
        cmd = 'ffmpeg -i '+ video_file + ' -i ' + audio_file + ' -shortest ' + muxvideo_file
        subprocess.call(cmd, shell=True); print('Saving and Muxing Done')
        
    def send_movie(self):
        emailling(self.From, self.Email_Var.get(), self.PWD, self.savedir, self.Name_Var.get().split()[-1] + "(" + self.Email_Var.get()  + ")_movie.mp4")

    def Main_destory(self):
        self.Frame_1.destroy(); self.Frame_2.destroy(); self.Frame_3.destroy();

    def Main_recreate(self, event):
        self.initialize(); self.Sim_Create()

    def start(self):
        self.progress_var.set(0); self.frames = 0; self.maxframes = 0
        
        
        Simu_Dir = self.model_select()+"/Dens-Maps/"
        cosmo = wCDM(70.3, self.Omega_m_Var.get(), self.Omega_l_Var.get(), w0=self.wx)
        filenames=sorted(glob.glob(self.simdir + "/" + Simu_Dir +'*.npy')); lga = linspace(log(0.05), log(1.0), 300); a = exp(lga); z = 1./a - 1.0; lktime = cosmo.lookback_time(z).value

        def animate(filename):
            image = np.load(filename); indx = filenames.index(filename)#; image=ndimage.gaussian_filter(image, sigma= sigmaval, truncate=truncateval, mode='wrap')
            im.set_data(image + 1)#; im.set_clim(image.min()+1.,image.max()+1.)
            self.time.set_text('%s %s' %(round(lktime[indx], 4), text_dict['t49']))
            return im

        dens_map = load(filenames[0])#;  dens_map=ndimage.gaussian_filter(dens_map, sigma= sigmaval, truncate=truncateval, mode='wrap') #; dens_map0 = load(filenames[-1]); #print dens_map0.min()+1, dens_map0.max()+1.
        im = self.ax.imshow(dens_map + 1, cmap=matplotlib.cm.magma, norm=matplotlib.colors.LogNorm(vmin=1., vmax=1800., clip = True), interpolation="bicubic")#, clim = (1, 1800.+1.))

        self.ax.annotate(text_dict['t42'] + self.Name_Var.get(), xy=(0.25, 0.45), fontsize='12', fontstyle = 'oblique', color='white', xycoords='data', xytext=(10., 40.), textcoords='data')
        self.time = self.ax.text(0.1, 0.05 , text_dict['t43'] + ' %s Gyr' %round(lktime[0], 4), horizontalalignment='left', verticalalignment='top',color='white', transform = self.ax.transAxes, fontsize=10)


        arr_hand = mpimg.imread(CWD + "/tmp/" + self.img_filename + "_Photo.jpg")
        imagebox = OffsetImage(arr_hand, zoom=.04); xy = [0.30, 0.45] # coordinates to position this image

        ab = AnnotationBbox(imagebox, xy, xybox=(40., -60.), xycoords='data', boxcoords="offset points", pad=0.1)
        self.ax.add_artist(ab)
        
        
        arr_hand1 = mpimg.imread("SIMCODE.png")
        imagebox1 = OffsetImage(arr_hand1, zoom=.1); xy = [950.0, 85.0] # coordinates to position this image
        ab1 = AnnotationBbox(imagebox1, xy, xybox=(0., 0.), xycoords='data', boxcoords="offset points", pad=0.1)
        self.ax.add_artist(ab1)

        #iMpc = lambda x: x*1024/125  #x in Mpc, return in Pixel *3.085e19
        ob = AnchoredHScaleBar(size=0.1, label="10Mpc", loc=4, frameon=False, pad=0.6, sep=2, color="white", linewidth=0.8)
        self.ax.add_artist(ob)
        
        sim_details_text = 'Universe Type: %s\nDark Energy Type: %s\nDark Matter Type: %s\nPrimordial Universe: %s\nGravity Type: %s' %('Flat', 'Constant', 'Cold', 'Natural', 'Einstein')
        
        self.ax.text(0.7, 0.7, sim_details_text, color='white', bbox=dict(facecolor='none', edgecolor='white', boxstyle='round,pad=1', alpha=0.5), transform = self.ax.transAxes, alpha = 0.5)
        
        #self.canvas.mpl_connect('button_press_event', self.onClick)
        self.ani = animation.FuncAnimation(self.fig, animate, filenames, repeat=False, interval=25, blit=False)
        self.ax.axis('off'); self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False); self.canvas.show()

        self.progress["value"] = 0; self.maxframes = 300; self.progress["maximum"] = 300
        self.read_frames()

    def read_frames(self):
        self.frames += 1
        #self.progress["value"] = self.frames
        self.progress_var.set(self.frames)
        if self.frames < self.maxframes:
            self._job1 = self.after(25, self.read_frames)

    def open_cam(self):
        # if self._job is not None:
        #     return None
        if not self.Name_Var.get():
            self.statusVariable.set(text_dict['t44'])
            return None
        self.show_frame()

    def show_frame(self):
        _, frame = cap.read()
        frame = cv2.flip(frame, 1)
        self.cv2image = cv2.cvtColor(frame, cv2.COLOR_BGR2RGBA)

        self.img = Image.fromarray(self.cv2image)

        self.lmain = self.ax.imshow(self.img); self.canvas.show(); self.ax.clear(); self.ax.axis('off')
        self._job = self.after(25, self.show_frame)
        #print self._job

    def saveImage(self):
        if not self.Name_Var.get():
            self.statusVariable.set(text_dict['t44'])
            return None

        try:
            os.stat(CWD + "/tmp/")
        except:
            os.mkdir(CWD + "/tmp/")

        self.img = self.img.resize((1024, 1024), Image.ANTIALIAS)
        self.img_filename = self.Name_Var.get(); self.img_filename = ''.join(e for e in self.img_filename if e.isalnum())
        
        self.img.save(CWD + "/tmp/" + self.img_filename + "_Photo.jpg")
        
        #----: CAM stop
        if self._job is not None:
            self.after_cancel(self._job)
            self._job = None
        #self.cap.release

    def callback_RunSpicf(self):
        self.main_nb.select(self.Main_page); self.NewRun_nb.select(self.NewRun_page)


    def models_refresh(self, *args):
        #return None
        if self.Lambda_Var.get() != 'Lambda_':
            self.CDM_Var.set('Lambda_'); self.IniM_Var.set('Lambda_'); self.MG_Var.set('Lambda_')
    def models_refresh1(self, *args):
        if self.CDM_Var.get() != 'Lambda_':
            self.Lambda_Var.set('Lambda_'); self.IniM_Var.set('Lambda_'); self.MG_Var.set('Lambda_')
    def models_refresh2(self, *args):
        if self.IniM_Var.get() != 'Lambda_':
            self.CDM_Var.set('Lambda_'); self.Lambda_Var.set('Lambda_'); self.MG_Var.set('Lambda_')
    def models_refresh3(self, *args):
        if self.MG_Var.get() != 'Lambda_':
            self.CDM_Var.set('Lambda_'); self.IniM_Var.set('Lambda_'); self.Lambda_Var.set('Lambda_')


    def ModelDirectory(self):
        self.modeldirname = tkFileDialog.askdirectory(parent=root, initialdir='/Users/mahmoud/')
        self.simdir = self.modeldirname

    def SaveDirectory(self):
        self.savedirname = tkFileDialog.askdirectory(parent=root, initialdir='/Users/mahmoud/')
        self.savedir = self.savedirname


    def resize(self, event):
        size = (event.width, event.height)
        self.Background_photo = ImageTk.PhotoImage(self.original.resize(size, Image.ANTIALIAS))
        self.Welcome_Frame.delete("IMG")
        self.Welcome_Frame.create_image(0, 0, image = self.Background_photo, anchor=NW, tags="IMG")



    def myDialog(self):
        Entry_dialog = Toplevel(self.master); Entry_dialog.title(text_dict['t45'])#; center(Entry_dialog)
        
        for x in range(2):
            Grid.columnconfigure(Entry_dialog, x, weight=2)
        
        for y in range(3):
            Grid.rowconfigure(Entry_dialog, y, weight=2)
        
        Label(Entry_dialog, text=text_dict['t46']).grid(row=0, column=0, sticky= W)
        self.UserName_Var = StringVar()
        self.UserName = Entry(Entry_dialog, textvariable=self.UserName_Var)
        self.UserName.grid(row=0, column=1, sticky= W+E+N+S)
        
        Label(Entry_dialog, text=text_dict['t47']).grid(row=1, column=0, sticky= W)
        self.PassWord_Var = StringVar()
        self.PassWord = Entry(Entry_dialog, textvariable=self.PassWord_Var, show="*")
        self.PassWord.grid(row=1, column=1, sticky= W+E+N+S)

        def getDate():
            self.From = self.UserName_Var.get(); self.PWD =self.PassWord_Var.get()
            Entry_dialog.destroy() # close the window
        
        submit = Button(Entry_dialog, text =text_dict['t48'], command = getDate)
        submit.grid(row=3, column=0,columnspan=2, sticky= W+E+N+S)

#--------- RUN ----------------------------
if __name__ == "__main__":
    root = Tk()
    app = Application(master=root)
    app.mainloop()
