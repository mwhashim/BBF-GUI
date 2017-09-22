from __future__ import division
import os, sys
sys.setrecursionlimit(50000) # to solve maximum recursion depth exceeded error !!

import logging
import glob
from threading import Thread
#from KThread import *

import time

from numpy import *
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
#import tkFileDialog, Tkconstants
#import tkinter.filedialog, tkinter.messagebox

if sys.version_info[0] < 3:
    from Tkinter import *
    from Tkinter import _setit
else:
    from tkinter import *
    from tkinter import _setit


import cv2
from PIL import Image, ImageTk
#from CosmoSuite import *

from astropy.cosmology import wCDM

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

#--- CWD -----
#CWD = os.getcwd(); print CWD
CWD = os.path.dirname(os.path.abspath(__file__))

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
        self.master.title("Big Bang Factory")

        #---------------------------------------
        for r in range(6):
            self.master.rowconfigure(r, weight=2)
        for c in range(5):
            self.master.columnconfigure(c, weight=2)
        #---------------------------------------
        self.master.rowconfigure(0, weight=0)
        #self.master.columnconfigure(0, weight=0)

        # Frame 0 :------------------------------------------------------
        self.Frame_0 = Frame(self.master, bg="white smoke")
        self.Frame_0.grid(row = 0, column = 0, rowspan = 5, columnspan = 5, sticky = W+E+N+S)
        
        Background_photo = PhotoImage(file="cosmicenviro.gif")
        Welcome_Frame = Label(self.Frame_0, image = Background_photo)
        Welcome_Frame.photo = Background_photo
        Welcome_Frame.pack() 
        
        Welcome_buttom = Button(self.Frame_0, text=u"Click to Create Your Universe", command=self.Sim_Create)
        Welcome_buttom.pack()
        
        #---------------------------------------
    def Sim_Create(self):
        self.Frame_0.destroy()
        self.initialize()
        self.Toolbar(self.Frame_1);
        self.PlotPan(self.Frame_3)
        self.MainFrame(self.Frame_2)

    def initialize(self):

        #--------------------------------------
        style = ttk.Style()
        style.layout('TNotebook.Tab', [])

        # Status Bar :-----
        self.statusVariable = StringVar()
        self.status = Label(self.master, textvariable = self.statusVariable, anchor = "w", fg = "yellow", bg = "blue")
        self.status.grid(column = 0, row = 7,columnspan = 5, sticky = 'EWS')
        self.statusVariable.set(u"Welcome to Big Bang Factory")

        # Frame 1 :------------------------------------------------------
        self.Frame_1 = PanedWindow(self.master, bg="white smoke")
        self.Frame_1.grid(row = 0, column = 0, rowspan = 1, columnspan = 5, sticky = W+E+N+S)

        # Frame 2 :-----------------------------------------------------
        self.Frame_2 = Frame(self.master, bg="white", bd= 1, relief= RIDGE)
        self.Frame_2.grid(row = 1, column = 0, rowspan = 5, columnspan = 2, sticky = W+E+N+S)

        # Frame 3 :-----------------------------------------------------
        self.Frame_3 = Frame(self.master, bg="white", bd= 3, relief= GROOVE)
        self.Frame_3.grid(row = 1, column = 2, rowspan = 5, columnspan = 3, sticky = W+E+N+S)
#
#        # Frame 4 :------------------------------------------------------
#        self.Frame_4 = Frame(self.master, bg="white smoke", bd= 5, relief= RIDGE)
#        self.Frame_4.grid(row = 4, column = 2, rowspan = 2, columnspan = 3, sticky = W+E+N+S)

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
#        self.f, (ax1, ax2) = subplots(2, sharex=True)
#        self.gs = gridspec.GridSpec(2,1, height_ratios=[3,1])
#        self.ax1 = subplot(self.gs[0]); #self.ax2 = subplot(gs[1])
#        #self.ax3 = axes([.16, .35, .3, .3])
#        self.f.subplots_adjust(hspace=0)#; self.f.tight_layout()

        self.f = plt.Figure()
        self.ax1 = self.f.add_subplot(111); self.ax1.axis('off')

        self.canvas = FigureCanvasTkAgg(self.f, master = frame)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, frame)
        self.canvas.get_tk_widget().grid(column = 0, row = 0, pady = 5)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

    def MainFrame(self, frame):

        #------------------------
        self.main_nb = ttk.Notebook(frame, padding = -5)
        self.Main_page = ttk.Frame(self.main_nb); self.Anals_page = ttk.Frame(self.main_nb)

        #------------------------
        self.main_nb.add(self.Main_page, text='Main')
        self.main_nb.add(self.Anals_page, text='Analysis & Results')
        #self.main_nb.add(self.SysPerf_page, text='System Preferences')
        self.main_nb.pack(side="top", fill="both", expand=True)

        self.main_nb.select(self.Main_page)
        #------------------------
        self.MainPage(self.Main_page)#; self.AnalysResultsPage(self.Anals_page)

    def MainPage(self, page):

        for r in range(5):
            page.rowconfigure(r, weight=2)

        for x in range(2):
            page.columnconfigure(x, weight=2)

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

        self.Name_Dict_group = LabelFrame(self.NewRun_page, text = "User Specifications")
        self.Name_Dict_group.grid(row = 0, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 0, weight=0)

        for x in range(2):
            Grid.columnconfigure(self.Name_Dict_group, x, weight=2)
        for y in range(2):
            Grid.rowconfigure(self.Name_Dict_group, y, weight=2)

        Label(self.Name_Dict_group, text="Name").grid(row=0, column=0, sticky= W)
        self.Name_Var = StringVar()
        self.User_Name = Entry(self.Name_Dict_group, textvariable=self.Name_Var)
        self.User_Name.grid(row=0, column=1, sticky= W+E+N+S,  columnspan = 6)
        #self.name_Var.set(None)

        Label(self.Name_Dict_group, text="Email").grid(row=1, column=0, sticky= W)
        self.Email_Var = StringVar()
        self.User_Email = Entry(self.Name_Dict_group, textvariable=self.Email_Var)
        self.User_Email.grid(row=1, column=1, sticky= W+E+N+S,  columnspan = 6)
        #self.Email_Var.set(None)

        self.Image_button = Button(self.Name_Dict_group, text="Show CAM", command = self.open_cam)
        self.Image_button.grid(row = 2, column = 2, columnspan = 1 , sticky = W+E+N+S, pady = 5)

        self.click_button = Button(self.Name_Dict_group, text="Take Photo", command = self.saveImage)
        self.click_button.grid(row = 2, column = 3, columnspan = 1 , sticky = W+E+N+S, pady = 5)

        self.Cosmo_Parms_group = LabelFrame(self.NewRun_page, text = "Cosmological Parameters")
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


        self.DarkEnergy_group = LabelFrame(self.NewRun_page, text = "Dark Energy Type")
        self.DarkEnergy_group.grid(row = 2, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.DarkEnergy_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.DarkEnergy_group, y, weight=2)


        self.Lambda_Var = StringVar()
        self.Lambda_Var.trace('w', self.models_refresh)
        self.Lambda_RadBtt = Radiobutton(self.DarkEnergy_group, text = 'Constant', variable = self.Lambda_Var, value = 'Lambda_')
        self.Lambda_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.Lambda_RadBtt = Radiobutton(self.DarkEnergy_group, text = 'Quintessence', variable = self.Lambda_Var, value = 'Quint_')
        self.Lambda_RadBtt.grid(row = 0, column = 1, sticky = W)

        self.Lambda_RadBtt = Radiobutton(self.DarkEnergy_group, text = 'Phantom', variable = self.Lambda_Var, value = 'Phantom_')
        self.Lambda_RadBtt.grid(row = 0, column = 2, sticky = W)
        self.Lambda_Var.set('Lambda_')

        self.DarkMatter_group = LabelFrame(self.NewRun_page, text = "Dark Matter Type")
        self.DarkMatter_group.grid(row = 3, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.DarkMatter_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.DarkMatter_group, y, weight=2)

        self.CDM_Var = StringVar()
        self.CDM_Var.trace('w', self.models_refresh)
        self.CDM_RadBtt = Radiobutton(self.DarkMatter_group, text = 'Cold', variable = self.CDM_Var, value = 'Lambda_')
        self.CDM_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.CDM_RadBtt = Radiobutton(self.DarkMatter_group, text = 'Warm', variable = self.CDM_Var, value = 'wDM_0.5-')
        self.CDM_RadBtt.grid(row = 0, column = 1, sticky = W)
        self.CDM_Var.set('Lambda_')

        self.IniMatter_group = LabelFrame(self.NewRun_page, text = "Initial Matter Distribution")
        self.IniMatter_group.grid(row = 4, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.IniMatter_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.IniMatter_group, y, weight=2)

        self.IniM_Var = StringVar()
        self.IniM_Var.trace('w', self.models_refresh)
        self.IniM_RadBtt = Radiobutton(self.IniMatter_group, text = 'Gaussian', variable = self.IniM_Var, value = 'Lambda_')
        self.IniM_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.IniM_RadBtt = Radiobutton(self.IniMatter_group, text = 'Positive non-Gaussian', variable = self.IniM_Var, value = 'LocalPNG_1000-')
        self.IniM_RadBtt.grid(row = 0, column = 1, sticky = W)

        self.IniM_RadBtt = Radiobutton(self.IniMatter_group, text = 'Negative non-Gaussian', variable = self.IniM_Var, value = 'LocalPNG_-1000-')
        self.IniM_RadBtt.grid(row = 0, column = 2, sticky = W)
        self.IniM_Var.set('Lambda_')

        self.MG_group = LabelFrame(self.NewRun_page, text = "Gravity Type")
        self.MG_group.grid(row = 5, column = 0, columnspan = 6, sticky = W+E+N+S)
        Grid.rowconfigure(self.NewRun_page, 1, weight=0)

        for x in range(13):
            Grid.columnconfigure(self.MG_group, x, weight=2)

        for y in range(4):
            Grid.rowconfigure(self.MG_group, y, weight=2)

        self.MG_Var = StringVar()
        self.MG_Var.trace('w', self.models_refresh)
        self.MG_RadBtt = Radiobutton(self.MG_group, text = 'Einstein', variable = self.MG_Var, value = 'Lambda_')
        self.MG_RadBtt.grid(row = 0, column = 0, sticky = W)

        self.MG_RadBtt = Radiobutton(self.MG_group, text = 'Modified Garvity', variable = self.MG_Var, value = 'MGfR_-1e-04-')
        self.MG_RadBtt.grid(row = 0, column = 1, sticky = W)
        self.MG_Var.set('Lambda_')


        #------------------------
        self.Control_Frame = Frame(page, bg="white smoke", bd= 1, relief= RIDGE)
        self.Control_Frame.grid(row = 4, column = 0, rowspan = 1, sticky = W+E+N+S)

        for x in range(4):
            Grid.columnconfigure(self.Control_Frame, x, weight=2)

        for y in range(2):
            Grid.rowconfigure(self.Control_Frame, y, weight=2)

        Grid.rowconfigure(self.Control_Frame, 0, weight=0)
        self.Simulation_Run = Button(self.Control_Frame, text = u"Simulation Run", foreground = 'red', command = self.start)
        self.Simulation_Run.grid(column = 0, row = 0, pady = 5, sticky= W+E+N+S)

        self.Save_Run = Button(self.Control_Frame, text = u"Save Movie", foreground = 'red', command = self.save_movie)
        self.Save_Run.grid(column = 1, row = 0, pady = 5, sticky= W+E+N+S)

        self.Reset_Run = Button(self.Control_Frame, text = u"Reset", foreground = 'red', command = self.run_reset)
        self.Reset_Run.grid(column = 4, row = 0, pady = 5, sticky= W+E+N+S)

        self.progress_var = DoubleVar()
        self.progress = ttk.Progressbar(self.Control_Frame, variable=self.progress_var,  orient="horizontal", length=300, mode="determinate")
        self.progress.grid(column = 0, row = 1, pady = 5, sticky= W+E+N+S, columnspan = 6)


    def run_reset(self):
        self.progress_var.set(0); self.ax1.clear(); self.ax1.axis('off'); self.canvas.show()

    def onClick(self, event):
        global pause
        pause == True

    def model_select(self):
        if self.Lambda_Var.get() != 'Lambda_':
            run_type = self.Lambda_Var.get()
            if run_type == "Quint_":
                self.wx = -0.9
            else:
                self.wx = -1.1

        elif self.CDM_Var.get()  != 'Lambda_':
            run_type = self.CDM_Var.get()
        elif self.IniM_Var.get() != 'Lambda_':
            run_type = self.IniM_Var.get()
        elif self.MG_Var.get() != 'Lambda_':
            run_type = self.MG_Var.get()
        else:
            run_type = 'Lambda_'; self.wx = -1.0

        if self.Omega_m_Var.get() == 0.0:
            Omega_m = 0.1; Omega_m = str(Omega_m)
        else:
            Omega_m = str(self.Omega_m_Var.get())

        model  = "BBF_" + run_type + Omega_m + "-" + str(self.Omega_l_Var.get())
        print model

        return model

    def save_movie(self):
        writer = animation.writers['ffmpeg'](fps=15)
        self.ani.save(CWD + "/users_photo/" + self.Name_Var.get().split()[-1] + "'s_movie.mp4", writer=writer, dpi=dpi)

    def start(self):
        self.frames = 0; self.maxframes = 0
        Simu_Dir = "/BBF/"+self.model_select()+"/Dens-Maps/"
        cosmo = wCDM(70.3, self.Omega_m_Var.get(), self.Omega_l_Var.get(), w0=self.wx)
        filenames=sorted(glob.glob(CWD + Simu_Dir +'*.npy')); lga = linspace(log(0.05), log(1.0), 300); a = exp(lga); z = 1./a - 1.0; lktime = cosmo.lookback_time(z).value

        def animate(filename):
            image = np.load(filename); indx = filenames.index(filename) #; image=ndim.gaussian_filter(image,sigma=sys.argv[2],mode='wrap')
            im.set_data(image + 1)#; im.set_clim(image.min()+1.,image.max()+1.)
            self.time.set_text('Age of the Universe: %s Gyr' %round(lktime[indx],3))
            return im

        dens_map = load(filenames[0]); dens_map0 = load(filenames[-1]); print dens_map0.min()+1, dens_map0.max()+1.
        im = self.ax1.imshow(dens_map + 1, cmap=matplotlib.cm.magma, norm=matplotlib.colors.LogNorm(vmin=1., vmax=1800., clip = True), interpolation="bicubic")#, clim = (1, 1800.+1.))

        self.ax1.annotate("This is the Universe by " + self.Name_Var.get(), xy=(0.25, 0.45), fontsize=10, color='white', xycoords='data', xytext=(10., 40.), textcoords='data')
        self.time = self.ax1.text(0.15, 0.1 , 'Age of the Universe: %s Gyr' %round(lktime[0],3), horizontalalignment='left', verticalalignment='top',color='white', transform = self.ax1.transAxes, fontsize=10)


        arr_hand = mpimg.imread(CWD + "/users_photo/" + self.Name_Var.get().split()[-1] + "'s_Photo.jpg")
        imagebox = OffsetImage(arr_hand, zoom=.04); xy = [0.30, 0.45] # coordinates to position this image

        ab = AnnotationBbox(imagebox, xy, xybox=(30., -40.), xycoords='data', boxcoords="offset points", pad=0.1)
        self.ax1.add_artist(ab)

        #iMpc = lambda x: x*1024/125  #x in Mpc, return in Pixel *3.085e19
        ob = AnchoredHScaleBar(size=0.1, label="10Mpc", loc=4, frameon=False, pad=0.6, sep=2, color="white", linewidth=0.8)
        self.ax1.add_artist(ob)

        self.canvas.mpl_connect('button_press_event', self.onClick)
        self.ani = animation.FuncAnimation(self.f, animate, filenames, repeat=False, interval=25, blit=False)
        self.ax1.axis('off'); self.canvas.show()

        self.progress["value"] = 0
        self.maxframes = 300
        self.progress["maximum"] = 300
        self.read_frames()

    def read_frames(self):
        self.frames += 1
        #self.progress["value"] = self.frames
        self.progress_var.set(self.frames)
        if self.frames < self.maxframes:
            self.after(25, self.read_frames)

    def open_cam(self):
        # if self._job is not None:
        #     return None
        if not self.Name_Var.get():
            self.statusVariable.set(u"Please Enter Your Name !!")
            return None
        self.show_frame()

    def show_frame(self):
        _, frame = cap.read()
        frame = cv2.flip(frame, 1)
        self.cv2image = cv2.cvtColor(frame, cv2.COLOR_BGR2RGBA)

        self.img = Image.fromarray(self.cv2image)

        self.lmain = self.ax1.imshow(self.img); self.canvas.show(); self.ax1.clear(); self.ax1.axis('off')
        self._job = self.after(25, self.show_frame)
        #print self._job

    def saveImage(self):
        if not self.Name_Var.get():
            self.statusVariable.set(u"Please Enter Your Name !!")
            return None

        try:
            os.stat(CWD + "/users_photo/")
        except:
            os.mkdir(CWD + "/users_photo/")

        self.img = self.img.resize((1024, 1024), Image.ANTIALIAS)
        self.img.save(CWD + "/users_photo/" + self.Name_Var.get().split()[-1] + "'s_Photo.jpg")
        #----: CAM stop
        if self._job is not None:
            self.after_cancel(self._job)
            self._job = None
        #self.cap.release

    def callback_RunSpicf(self):
        self.main_nb.select(self.Main_page); self.NewRun_nb.select(self.NewRun_page)


    def models_refresh(self, *args):
        #return None
        if self.Lambda_Var.get() != 'Lambda_':#  or self.MG_Var.get() != 'Lambda_' or self.IniM_Var.get() != 'Lambda_' or self.CDM_Var.get() != 'Lambda_':
            self.CDM_Var.set('Lambda_'); self.IniM_Var.set('Lambda_'); self.MG_Var.set('Lambda_')
        elif self.CDM_Var.get() != 'Lambda_':# or self.Lambda_Var.get() != 'Lambda_' or self.MG_Var.get() != 'Lambda_' or self.IniM_Var.get() != 'Lambda_':
            self.Lambda_Var.set('Lambda_'); self.IniM_Var.set('Lambda_'); self.MG_Var.set('Lambda_')
        elif self.IniM_Var.get() != 'Lambda_':# or self.MG_Var.get() != 'Lambda_' or self.Lambda_Var.get() != 'Lambda_' or self.CDM_Var.get() != 'Lambda_' :
            self.CDM_Var.set('Lambda_'); self.Lambda_Var.set('Lambda_'); self.MG_Var.set('Lambda_')
        elif self.MG_Var.get() != 'Lambda_':# or self.IniM_Var.get() != 'Lambda_' or self.Lambda_Var.get() != 'Lambda_' or self.CDM_Var.get() != 'Lambda_':
            self.CDM_Var.set('Lambda_'); self.IniM_Var.set('Lambda_'); self.Lambda_Var.set('Lambda_')
#        else:
#            self.Lambda_Var.set('Lambda_')

#--------- RUN ----------------------------
if __name__ == "__main__":
    root = Tk()
    app = Application(master=root)
    app.mainloop()
