#################################################################
### Optimize and plot the gain/VSWR/Pattern of a 3-element Yagi antenna
###
### Adapted from: Examples of using pyNEC (necpp) python libraries
### 
### Ref: Timothy C.A. Molteno, 
###      ''NEC2++: An NEC-2 compatible Numerical Electromagnetics Code'', 
### 	    Electronics Technical Reports No. 2014-3, ISSN 1172-496X, October 2014.
### 
### See: http://tmolteno.github.io/necpp/libnecpp_8h.html
###      http://astroelec.blogspot.com/2015/05/modeling-antennas-in-python-with-nec2.html
###    	 https://github.com/tmolteno/necpp
###      https://github.com/tmolteno/python-necpp
###      https://github.com/tmolteno/python-necpp/blob/master/PyNEC/example/logperiodic_opt.py
### 
### Install using: pip install --user PyNEC==1.7.3.6
#################################################################
# TODO: Add ground capability
# TODO: Add image of antenna structure with dimensions
# TODO: Add handler to respond appropriately when there is  geometry problem (e.g. element diameter too large)
# TODO: Investigate why results differ from 4NEC2 and EZNEC for vswr, gain, etc.
# TODO: Determine why pattern is sometimes mirror imaged (compared to 4NEC2, etc.)



import sys
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QApplication, QMainWindow, QSpinBox, QDoubleSpinBox, QMessageBox, QRadioButton, QFileDialog

#################################################################
#### import matplotlib classes to enable plots in PyQt5
#################################################################
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import mplcursors  # used to display info when hovering over plot lines

import logging

import numpy as np
import scipy.optimize
import math

from PyNEC import *
from antenna_util import *
from context_clean import *


logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARNING)
version = "1.0"
date = "Nov 27, 2024"

# Global iteration counter to show optimization progress
iteration_count = 0

# Global variables used in the optimization mapping plot
sampledResults = []
sampledL1 = []
sampledL2 = []
sampledL3 = []
sampledD1 = []
sampledD2 = []
        



#################################################################
class MplFigure(object):
    def __init__(self, parent):
        self.figure = plt.figure(facecolor='white')
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, parent) 
#################################################################

        
#################################################################
class uiMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        uic.loadUi('Yagi Optimizer UI.ui', self)  #### load Qt Designer UI design
        self.setWindowTitle("3-Element Yagi Optimizer v" + version)
        self.status_label.setText("Status: Initial Unoptimized Design")

        #### Populate default values into spinner & combo boxes ####
        self.system_impedance = 50 # Hard coded for now

        self.foldedDipoleSpacing = 0.05 # distance between driven element and folded element in meters 

        self.design_freq_mhz = 144.1 # The frequency of interest (in MHz)

        self.start_frq = 130.0 # min frequency for freq range plots
        self.stop_frq  = 160.0 # max frequency for freq range plots (inclusive)
        self.step_frq = 1.0 # frequency step for plots

        self.start_frq_opt = 144.0 # start of optimization freq range
        self.stop_frq_opt = 146.0  # end of optimization freq range (inclusive)
        self.step_frq_opt = 1.0 # frequncy steps for optimization calculation

        self.elementDiameter = 0.125 # 1/8-inch diameter (3 mm)
        self.elementMaterial = "Aluminum"

        vswrWeight = 30    # default weights for initialization
        fwdGainWeight = 10
        fbRatioWeight = 1

        self.DesignFrqDoubleSpinBox.setValue(self.design_freq_mhz)

        self.OptFrqStrtDoubleSpinBox.setValue(self.start_frq_opt)
        self.OptFrqStopDoubleSpinBox.setValue(self.stop_frq_opt)
        self.OptFrqStepDoubleSpinBox.setValue(self.step_frq_opt)

        self.PlotFrqMinDoubleSpinBox.setValue(self.start_frq)
        self.PlotFrqMaxDoubleSpinBox.setValue(self.stop_frq)
        self.PlotFrqStepDoubleSpinBox.setValue(self.step_frq)

        self.ElementDiameterDoubleSpinBox.setValue(self.elementDiameter)

        self.FoldedDipoleSpcDoubleSpinBox.setValue(self.foldedDipoleSpacing / 0.0254) # meters to inches

        self.VswrWeightDoubleSpinBox.setValue(vswrWeight)
        self.FwdGainWeightDoubleSpinBox.setValue(fwdGainWeight)
        self.FBRatioWeightDoubleSpinBox.setValue(fbRatioWeight)


        #### Populate values into combo boxes  ####
        self.ElementMaterialComboBox.addItems(["Aluminum", "Brass", "Stainless Steel", "Copper"])
        self.ElementMaterialComboBox.setCurrentIndex(0) # default "Aluminum"
        
        self.OptAlgorithmComboBox.addItems(["Diff Evolution", "Basin Hopping", 
                                            "Dual Annealing", "Gradient Descent",
                                            "Dividing Rectangles", "Simplicial Homology", 
                                            "Brute Force"])
        self.OptAlgorithmComboBox.setCurrentIndex(0) # default "Diff Evolution"


       
        #### Connect signals from menu item selections ####
        self.actionOpen_Yagi_File.triggered.connect(lambda: self.clicked("Open Yagi File Was Clicked"))
        self.actionSave_Yagi_File.triggered.connect(lambda: self.clicked("Save Yagi File Was Clicked"))
        self.actionAbout.triggered.connect(lambda: self.clicked("About Was Clicked"))
        self.actionExit.triggered.connect(sys.exit)
        
        #### Connect signals from entry box changes ####
        self.DesignFrqDoubleSpinBox.valueChanged.connect(self.onValueChanged)

        self.OptFrqStrtDoubleSpinBox.valueChanged.connect(self.onValueChanged)
        self.OptFrqStopDoubleSpinBox.valueChanged.connect(self.onValueChanged)
        self.OptFrqStepDoubleSpinBox.valueChanged.connect(self.onValueChanged)

        self.PlotFrqMinDoubleSpinBox.valueChanged.connect(self.onValueChanged)
        self.PlotFrqMaxDoubleSpinBox.valueChanged.connect(self.onValueChanged)
        self.PlotFrqStepDoubleSpinBox.valueChanged.connect(self.onValueChanged)

        self.ElementDiameterDoubleSpinBox.valueChanged.connect(self.onValueChanged)
        self.ElementMaterialComboBox.currentIndexChanged.connect(self.onValueChanged)

        self.FoldedDipoleSpcDoubleSpinBox.valueChanged.connect(self.onValueChanged)

        self.OptAlgorithmComboBox.currentIndexChanged.connect(self.onValueChanged)
        self.VswrWeightDoubleSpinBox.valueChanged.connect(self.onValueChanged)
        self.FwdGainWeightDoubleSpinBox.valueChanged.connect(self.onValueChanged)
        self.FBRatioWeightDoubleSpinBox.valueChanged.connect(self.onValueChanged)

        self.useFoldedDipoleRadioButton.clicked.connect(self.onValueChanged)
        self.optimizePushButton.clicked.connect(self.updatePlot)

    
        #### Setup and initialize plotting area in the verticalLayout box ####
        self.main_figure = MplFigure(self)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.addWidget(self.main_figure.toolbar)
        self.verticalLayout.addWidget(self.main_figure.canvas)

        self.initMplWidget()

        self.show()



    #################################################################
    #### initialize MatPlotLib plot and update with default values
    #################################################################
    def initMplWidget(self):
        self.ax = self.main_figure.figure.add_subplot(111)
        self.onValueChanged() # call event handler that reads inputs
        self.initializePlot() # add plots of initial guess at dimensions 
    #################################################################


    #################################################################
    #### Method to handle clicking menu items and display
    #### text that is passed using the triggered/clicked connects
    #################################################################
    def clicked(self, text):
        self.setStatusTip(text)

        if text == 'Open Yagi File Was Clicked':
            fileName, _ = QFileDialog.getOpenFileName(self, "Open Yagi File", "",
                                                     "Yagi Optimizer Files (*.ygi);;All Files (*)")
            if fileName:
                try:
                    with open(fileName, "r") as yagiFile:
                        self.system_impedance = float(yagiFile.readline())
                        self.DesignFrqDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.PlotFrqMinDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.PlotFrqMaxDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.PlotFrqStepDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.OptFrqStrtDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.OptFrqStopDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.OptFrqStepDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.ElementDiameterDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.ElementMaterialComboBox.setCurrentIndex(int(yagiFile.readline()))
                        self.OptAlgorithmComboBox.setCurrentIndex(int(yagiFile.readline()))
                        self.VswrWeightDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.FwdGainWeightDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.FBRatioWeightDoubleSpinBox.setValue(float(yagiFile.readline()))
                        self.FoldedDipoleSpcDoubleSpinBox.setValue(float(yagiFile.readline()))

                    self.setStatusTip("File: " + fileName + " loaded successfully")

                except Exception:
                    QMessageBox.warning(self,"Yagi Optimizer - File Open Warning",
                        "Open: " + fileName + "\n\nFile not compatible with Yagi Optimizer v" + version)
                    self.setStatusTip("File: " + fileName + " not compatible")
        
            self.onValueChanged() # update to reflect read parameters


        if text == 'Save Yagi File Was Clicked':
            fileName, _ = QFileDialog.getSaveFileName(self,"Save Yagi File","*.ygi",
                                                      "Yagi Optimizer Files (*.ygi);;All Files (*)")
            if fileName:
                with open(fileName, "w") as yagiFile:
                    yagiFile.write(str(self.system_impedance)+"\n")
                    yagiFile.write(str(self.DesignFrqDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.PlotFrqMinDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.PlotFrqMaxDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.PlotFrqStepDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.OptFrqStrtDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.OptFrqStopDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.OptFrqStepDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.ElementDiameterDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.ElementMaterialComboBox.currentIndex()) + "\n")
                    yagiFile.write(str(self.OptAlgorithmComboBox.currentIndex()) + "\n")
                    yagiFile.write(str(self.VswrWeightDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.FwdGainWeightDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.FBRatioWeightDoubleSpinBox.value())+"\n")
                    yagiFile.write(str(self.FoldedDipoleSpcDoubleSpinBox.value()))

                self.setStatusTip("File: " + fileName + " saved successfully")


        if text == 'About Was Clicked':
            QMessageBox.about(self, "Yagi Optimizer", 
                                    "Yagi Optimizer \n" +
                                    "Version: " + version + "\n" +
                                    "Author: Steve Harbin \n" +
                                    "Date: " + date + "\n" + 
                                    "----------------------- \n" +
                                    "3-element yagi array optimizer \n" +
                                    "Based on necpp references " )

        # TODO: Add handlers for other menu actions if added                         
    #################################################################



    #################################################################
    #### This method handles changes in spinner boxes and combo boxes
    #### Enables Optimize button when changes occur
    #################################################################
    #@pyqtSlot()   # Do I need this??
    def onValueChanged(self):
        
        self.designFrq = self.DesignFrqDoubleSpinBox.value()

        self.optFrqStrt = self.OptFrqStrtDoubleSpinBox.value()
        self.optFrqStop = self.OptFrqStopDoubleSpinBox.value()
        self.optFrqStep = self.OptFrqStepDoubleSpinBox.value()

        self.plotFrqMin = self.PlotFrqMinDoubleSpinBox.value()
        self.plotFrqMax = self.PlotFrqMaxDoubleSpinBox.value()
        self.plotFrqStep = self.PlotFrqStepDoubleSpinBox.value()

        self.elementDiameter = self.ElementDiameterDoubleSpinBox.value()
        self.elementMaterial = self.ElementMaterialComboBox.currentText()

        self.foldedDipoleSpacing = self.FoldedDipoleSpcDoubleSpinBox.value() * 0.0254 # convert inches to meters
        self.foldedDipoleCheck()

        self.optAlgorithm = self.OptAlgorithmComboBox.currentText()         
        self.vswrWeight = self.VswrWeightDoubleSpinBox.value()
        self.fwdGainWeight = self.FwdGainWeightDoubleSpinBox.value()
        self.fbRatioWeight = self.FBRatioWeightDoubleSpinBox.value()

        self.optimizePushButton.setText('Optimize')
        self.optimizePushButton.setEnabled(True)

        #self.updatePlot() 
    #################################################################


    #################################################################
    ####  This method handles the folded dipole check box
    ####  Changes the spacing to meters (or zero if not used)
    #################################################################
    def foldedDipoleCheck(self):
        if self.useFoldedDipoleRadioButton.isChecked():
            self.foldedDipoleSpacing = self.FoldedDipoleSpcDoubleSpinBox.value() * 0.0254 # convert inches to meters
            self.FoldedDipoleSpcDoubleSpinBox.setEnabled(True)
        else: 
            self.foldedDipoleSpacing = 0
            self.FoldedDipoleSpcDoubleSpinBox.setEnabled(False)
        return 
    #################################################################




 
    #####################################################################
    #### Calculate and show results when initializing the User Interface
    #### Uses "rule of thumb" values for element lengths (not optimized)
    #####################################################################
    def initializePlot(self):
        self.design_freq_mhz = self.designFrq # The frequency of interest (in MHz)
        wavelength = 299792e3/(self.design_freq_mhz*1000000)

        self.start_frq = self.plotFrqMin # min frequency for freq range plots
        self.stop_frq  = self.plotFrqMax # max frequency for freq range plots (inclusive)
        self.step_frq  = self.plotFrqStep # frequency step for plots

        self.start_frq_opt = self.optFrqStrt # start of optimization freq range
        self.stop_frq_opt = self.optFrqStop  # end of optimization freq range (inclusive)
        self.step_frq_opt = self.optFrqStep # frequncy steps for optimization calculation
        '''
        initial_l1  = 1.05 * self.wavelength / 2 # reflector length initial guess .989
        initial_l2  = 0.98 * self.wavelength / 2 # driven element length initial guess .946
        initial_l3  = 0.88 * self.wavelength / 2 # director length initial guess .88
        initial_d1  = 0.60 * self.wavelength / 4 # distance reflector to driven element initial guess
        initial_d2  = 0.50 * self.wavelength / 4 # distance driven element to director initial guess
        '''
        initial_l1  = 0.487 * wavelength # reflector length initial guess 
        initial_l2  = 0.437 * wavelength # driven element length initial guess 
        initial_l3  = 0.450 * wavelength # director length initial guess 
        initial_d1  = 0.111 * wavelength # distance reflector to driven element initial guess
        initial_d2  = 0.184 * wavelength # distance driven element to director initial guess

        print("Wavelength is %0.4fm, initial geometry is L1=%0.4fm, L2=%0.4fm, L3=%0.4fm, D1=%0.4fm, D2=%0.4fm, Spc=%0.4fm" % 
              (wavelength, initial_l1, initial_l2, initial_l3, initial_d1, initial_d2, self.foldedDipoleSpacing))  
        
        print("\nUnoptimized antenna parameters...")
        self.show_report(initial_l1, initial_l2, initial_l3, initial_d1, initial_d2) # display results for initial guess
    #####################################################################################


    #####################################################################################
    #### This method is called when the Optimize button is clicked.  It Kicks off the 
    #### optimization process based on the algorithm selected.  Then recalculates
    ####  and calls the function to redraw plots                               
    #####################################################################################
    def updatePlot(self):

        self.design_freq_mhz = self.designFrq # The frequency of interest (in MHz)
        wavelength = 299792e3/(self.design_freq_mhz*1000000)

        self.start_frq = self.plotFrqMin # min frequency for freq range plots
        self.stop_frq  = self.plotFrqMax # max frequency for freq range plots (inclusive)
        self.step_frq  = self.plotFrqStep # frequency step for plots

        self.start_frq_opt = self.optFrqStrt # start of optimization freq range
        self.stop_frq_opt = self.optFrqStop  # end of optimization freq range (inclusive)
        self.step_frq_opt = self.optFrqStep # frequncy steps for optimization calculation
        
        # Starting "rule of thumb" values for element lengths pre-optimization
        '''
        initial_l1  = 1.05 * self.wavelength / 2 # reflector length initial guess
        initial_l2  = 0.98 * self.wavelength / 2 # driven element length initial guess
        initial_l3  = 0.88 * self.wavelength / 2 # director length initial guess
        initial_d1  = 0.60 * self.wavelength / 4 # distance reflector to driven element initial guess
        initial_d2  = 0.50 * self.wavelength / 4 # distance driven element to director initial guess
        '''
        initial_l1  = 0.487 * wavelength # reflector length initial guess 
        initial_l2  = 0.437 * wavelength # driven element length initial guess 
        initial_l3  = 0.450 * wavelength # director length initial guess 
        initial_d1  = 0.111 * wavelength # distance reflector to driven element initial guess
        initial_d2  = 0.184 * wavelength # distance driven element to director initial guess

        print("\nOptimizing antenna...")
        self.optimizePushButton.setText('Optimizing...')
        self.optimizePushButton.setEnabled(False)

        target = self.create_optimization_target()

        # Assume that the initial guess is within +/-25% of the optimum dimensions and use as optimizatin bounds
        wire_radius = (self.elementDiameter)/2 * 0.0254 # diameter in inches, radius in meters
        foldedDipoleSpacing = self.foldedDipoleSpacing  # distance between driven element and folded element in meters
        if 1.5*initial_d2 <= foldedDipoleSpacing + 2*wire_radius: # max d2 interferes with folded element
            logging.critical("ERROR: Frequency too high for folded dipole spacing of %6.3f in." % (foldedDipoleSpacing/0.0254))
            logging.critical("       Results cannot be fully optimized")
            QMessageBox.critical(self,"ERROR!", 
                                "Frequency too high for folded dipole spacing of %6.3f in. \nResults cannot be fully optimized.  Exiting application." % (foldedDipoleSpacing/0.0254))
            sys.exit() # not possible to set director element far enough away. So exit.
        elif 0.5*initial_d2 <= foldedDipoleSpacing + 2*wire_radius: # min d2 interferes with folded element, but max is ok
            logging.warning("WARNING: Frequency too high for folded dipole spacing of %6.3f in." % (foldedDipoleSpacing/0.0254))
            logging.warning("         Results may not be fully optimized.")
            QMessageBox.warning(self,"WARNING!", 
                                "Frequency too high for folded dipole spacing of %6.3f in. \nResults may not be fully optimized." % (foldedDipoleSpacing/0.0254))

        # Set minimum dimensions so there are no mechanical interference problems
        minL1 = max(0.75*initial_l1, 2*wire_radius)
        minL2 = max(0.75*initial_l2, 2*wire_radius)
        minL3 = max(0.75*initial_l3, 2*wire_radius)
        minD1 = max(0.5*initial_d1, 2*wire_radius+0.0001)
        minD2 = max(0.5*initial_d2, foldedDipoleSpacing+2*wire_radius+0.0001)

        maxL1 = max(1.25*initial_l1, 2*wire_radius)
        maxL2 = max(1.25*initial_l2, 2*wire_radius)
        maxL3 = max(1.25*initial_l3, 2*wire_radius)
        maxD1 = max(1.5*initial_d1, 2*wire_radius+0.0002)
        maxD2 = max(1.5*initial_d2, foldedDipoleSpacing+2*wire_radius+0.0002)

        bounds = [ (minL1, maxL1), (minL2, maxL2), (minL3, maxL3), (minD1, maxD1), (minD2, maxD2) ]
        
        if self.optAlgorithm == "Gradient Descent":
            # Use gradient descent: (finds local minimum only; Fast )
            optimized_result = scipy.optimize.minimize(target, np.array([initial_l1, initial_l2, initial_l3, initial_d1, initial_d2]), method='Nelder-Mead', bounds=bounds)

        elif self.optAlgorithm == "Diff Evolution":
            # Use differential evolution: (Good results @ 70626 iterations)
            optimized_result = scipy.optimize.differential_evolution(target, bounds, seed=42, disp=False, popsize=20)

        elif self.optAlgorithm == "Basin Hopping":
            # Use basin hopping: (works good. May need to adjust parameters T, step size etc.)
            minimizer_kwargs = dict(method='Nelder-Mead')
            optimized_result = scipy.optimize.basinhopping(target, np.array([initial_l1, initial_l2, initial_l3, initial_d1, initial_d2]), minimizer_kwargs=minimizer_kwargs, niter=10, stepsize=0.015, T=2.0, seed=42, disp=True)
            
        elif self.optAlgorithm == "Dual Annealing":
            # Use Dual Annealing: 
            minimizer_kwargs = dict(method='Nelder-Mead')
            optimized_result = scipy.optimize.dual_annealing(target, bounds, seed=42, minimizer_kwargs=minimizer_kwargs)

        elif self.optAlgorithm == "Dividing Rectangles":
            # Use Dividing Rectangles: (Not very good with these parameters)
            optimized_result = scipy.optimize.direct(target, bounds)

        elif self.optAlgorithm == "Simplicial Homology":
            # Use simplicial homology global optimization
            minimizer_kwargs = dict(method='Nelder-Mead')
            optimized_result = scipy.optimize.shgo(target, bounds, minimizer_kwargs=minimizer_kwargs)

        elif self.optAlgorithm == "Brute Force":
            # Use brute force optimization: (With Ns=10: ~100k iterations. With Ns=14: >500k iterations [Ns^5]. Takes hours)
            ranges = (bounds[0], bounds[1],bounds[2], bounds[3], bounds[4])
            optimized_result_brute = scipy.optimize.brute(target, ranges, full_output=True, Ns=10)
            optimized_result = optimized_result_brute[0]
            optimized_l1, optimized_l2, optimized_l3, optimized_d1, optimized_d2 =  optimized_result[0], optimized_result[1], optimized_result[2], optimized_result[3], optimized_result[4]

        print("\nOptimized antenna parameters...")
        if self.optAlgorithm != "Brute Force": # Brute Force algorithm returns a different structure than the others.
            optimized_l1, optimized_l2, optimized_l3, optimized_d1, optimized_d2 =  optimized_result.x[0], optimized_result.x[1], optimized_result.x[2], optimized_result.x[3], optimized_result.x[4]
        self.show_report(optimized_l1, optimized_l2, optimized_l3, optimized_d1, optimized_d2)

        self.status_label.setText("Status: Optimized Design Complete")
        #self.optimizePushButton.setStyleSheet("background-color: lightgreen") 
        self.optimizePushButton.setText('Optimize')  # Reset button label when complete
    #####################################################################################




    #####################################################################################
    # Geometry for a generic 3-element yagi parallel to the y-axis and pointing in the +x direction
    # Yagi dimensions are defined as:
    #     l1: length of reflector element in meters
    #     l2: length of driven element in meters
    #     l3: length of director element in meters
    #     d1: distance between reflector and driven elements in meters
    #     d1: distance between director and driven elements in meters
    #
    # The driven dipole is at the origin. It is a folded dipole with 5-cm spacing 
    # Dipoles are parallel to the y axis; the driven element is centered at (0, 0, 0)
    #####################################################################################
    def geometry_yagi(self, l1, l2, l3, d1, d2):
        wire_radius = (self.elementDiameter)/2 * 0.0254 # diameter in inches, radius in meters
        foldedDipoleSpacing = self.foldedDipoleSpacing  # distance between driven element and folded element in meters
        wavelength = 299792e3/(self.designFrq*1000000)

        if self.elementMaterial == "Brass":
            conductivity = 15600000 # mhos/m
        elif self.elementMaterial == "Aluminum":
            conductivity = 25000000 # mhos/m (6061-T6)
        elif self.elementMaterial == "Stainless Steel":
            conductivity = 1450000 # mhos/m
        elif self.elementMaterial == "Copper":
            conductivity = 57471264 # mhos/m
        else:
            return # should never get here (leave placeholder for additional materials)

        extTWKernel = True # use extended thin wire kernel in the simulation

        # calculate minimum segment length for good accuracy
        if extTWKernel:
            #length_segments_min = 2 * wire_radius # at non-junctions L/R>2
            length_segments_min = 6 * wire_radius # at junctions L/R>6
        else:
            length_segments_min = 8 * wire_radius

        # calculate maximum segment length for good accuracy
        length_segments_max = wavelength/18 # max length per Cebik. Others use wl/10

        # NOTE: for length_segments_min < length_segments_max,
        #       f_mhz < 2.77777 / radius_wire

        
        # set segment length to value within NEC2 accuracy bounds 
        if (length_segments_max > length_segments_min):
            length_segments = math.sqrt(length_segments_max * length_segments_min) # set to geometric mean
            #length_segments = length_segments_min # make as small as practical (slows things down at low freqs)
        else: # This happens if freq > 1,749.8 MHz for 1/8" diameter wire, 218.7 MHz for 1" diameter (see NOTE above)
            length_segments = length_segments_max # make as small as practical (max is smaller than min here)
            logging.warning("GEOMETRY WARNING: Element diameter too large at this frequency. length_segments_min=%6.3f >= length_segments_max=%6.3f" %(length_segments_min, length_segments_max))
            #QMessageBox.critical(self,"GEOMETRY ERROR!", 
            #                    " Element diameter too large at this frequency. \n\nCannot proceed.")
            #sys.exit()

        if length_segments < wavelength/1000:  
            length_segments = wavelength/1000 # set to minimum allowed
            logging.critical("GEOMETRY ERROR: Wire segment length too short (less than wl/1000). Results may be invalid.")
            QMessageBox.critical(self,"GEOMETRY ERROR!", 
                                 "Wire segment length too short (less than wl/1000). Results may be invalid.")
            # TODO: Maybe set length to wl/1000+.0000001 and don't exit
            sys.exit()


        nec = context_clean(nec_context())
        nec.set_extended_thin_wire_kernel(extTWKernel)
        geo = geometry_clean(nec.get_geometry())


        # Set mnimimum dimensions to ensure geometry is valid and no mechanical interference occurs
        # d2 must exceed (Wire diameter + foldedDipoleSpacing) to prevent mechanical interference
        # d1 must exceed the wire diameter
        #print ("number reflector segs = ", l1 / length_segments, l1, length_segments)
        if (l1 < 2*length_segments) or (math.isnan(l1)): l1=2*length_segments
        if (l2 < 2*length_segments) or (math.isnan(l2)): l2=2*length_segments
        if (l3 < 2*length_segments) or (math.isnan(l3)): l3=2*length_segments
        if d1 < wire_radius*2 + length_segments/1000: 
            logging.info("geometry_yagi: d1 too small - d1=%6.3f  element dia=%6.3f" % (d1, 2*wire_radius))
            d1 = wire_radius*2 + length_segments/1000
        if d2 < foldedDipoleSpacing + (wire_radius*2) + length_segments/1000: 
            logging.info("geometry_yagi: d2 too small - d2=%6.3f  element dia=%6.3f" % (d2, 2*wire_radius))
            d2 = foldedDipoleSpacing + (wire_radius*2) + length_segments/1000



        # define the elements' coordinates and parameters
        # Note that nec tags start at 1 (not zero-indexed)

        #### Reflector Element
        nr_segments = int(l1 / length_segments)
        if nr_segments % 2 == 0: nr_segments = nr_segments - 1 # Set to an odd number of segments (reduce segs)
        geo.wire(tag_id=1, nr_segments=nr_segments, src=[-d1, -l1/2, 0], dst=[-d1, l1/2, 0], radius=wire_radius)

        #### Driven element 
        nr_segments = int(l2 / length_segments)  # <--- MUST make this an odd number so that source can be in center
        if nr_segments % 2 == 0: nr_segments = nr_segments - 1 # Set to an odd number of segments (reduce segs)
        geo.wire(tag_id=2, nr_segments=nr_segments, src=[0, -l2/2, 0], dst=[0, l2/2, 0], radius=wire_radius)
        driver_center_seg = int(nr_segments/2) + 1
        
        #### Folded dipole element and connection segments
        if self.useFoldedDipoleRadioButton.isChecked():
            nrEndSegments = max(1, int(foldedDipoleSpacing/length_segments)) # ensure there are at minimum 1 segments
            # Check junctions (l_seg_big/l_seg_small<5 and wire_radius1/wire_radius2<5 to avoid errors)
            lengthEndSegments = foldedDipoleSpacing/nrEndSegments
            if max(lengthEndSegments, length_segments)/min(lengthEndSegments, length_segments) > 5: #(Error if Lbig/Lsmall>5 at junction)
                logging.critical("GEOMETRY ERROR: Folded dipole segment length too short at junctions. Results may be invalid.")
                logging.critical("ratio=%6.3f  segment 1: %6.3f, segment 2: %6.3f" % 
                    (max(lengthEndSegments, length_segments)/min(lengthEndSegments, length_segments), 
                    length_segments, lengthEndSegments))
                logging.critical("This error often occurs if the element radius is too large.")
                QMessageBox.critical(self,"GEOMETRY ERROR!", 
                                    "Folded dipole segment length too short at junctions. Results may be invalid.\n" + 
                                    "This error often occurs if the element radius is too large.")
                sys.exit()
            
            geo.wire(tag_id=4, nr_segments=nr_segments, src=[foldedDipoleSpacing, -l2/2, 0], dst=[foldedDipoleSpacing, l2/2, 0], radius=wire_radius)
            geo.wire(tag_id=5, nr_segments=nrEndSegments, src=[0, -l2/2, 0], dst=[foldedDipoleSpacing, -l2/2, 0], radius=wire_radius)
            geo.wire(tag_id=6, nr_segments=nrEndSegments, src=[0, l2/2, 0], dst=[foldedDipoleSpacing, l2/2, 0], radius=wire_radius)
        
        #### Director Element
        nr_segments = int(l3 / length_segments)
        if nr_segments % 2 == 0: nr_segments = nr_segments - 1 # Set to an odd number of segments (reduce segs)
        geo.wire(tag_id=3, nr_segments=nr_segments, src=[d2, -l3/2, 0], dst=[d2, l3/2, 0], radius=wire_radius)


        nec.set_wire_conductivity(conductivity)
        nec.geometry_complete(ground_plane=False)


        nec.voltage_excitation(wire_tag=2, segment_nr=driver_center_seg, voltage=1.0)

        return nec
    #####################################################################################



    #####################################################################################
    #### Perform an NEC2 analysis on a 3-element yagi passed to this method over a range
    #### of frequencies provided.  Returns lists of Freq, Fwd Gain, VSWR, and Rev Gain
    #### Yagi dimensions are defined as:
    ####     l1: length of reflector element in meters
    ####     l2: length of driven element in meters
    ####     l3: length of director element in meters
    ####     d1: distance between reflector and driven elements in meters
    ####     d1: distance between director and driven elements in meters
    #####################################################################################
    def get_gain_swr_range(self, l1, l2, l3, d1, d2, start=156.0, stop=158.0, step=1.0):
        fwd_gains = []
        frequencies = []
        vswrs = []
        rev_gains = []

        def frange(start, stop, step): # Using generator function to define range prevents memory problems
            i = start
            while i <= stop:
                yield i
                i += step
    
        for freq in frange(start, stop, step):
        #for freq in range(start, stop + 1, step):
            nec = self.geometry_yagi(l1, l2, l3, d1, d2) # must get new context to change freq
            nec.set_frequency(freq) 
            nec.radiation_pattern(thetas=Range(90, 90, count=1), phis=Range(0,360,count=2))
            # confusingly, phi range 0..360 with count=2 produces points at 0 and 180 (step=(360-0)/2 )

            rp = nec.context.get_radiation_pattern(0) 
            ipt = nec.get_input_parameters(0)
            z = ipt.get_impedance()
            
            fwd_gains.append(rp.get_gain()[0][0]) # Gain at phi = 0 degrees
            rev_gains.append(rp.get_gain()[0][1]) # Gain at phi = 180 degrees
            vswrs.append(vswr(z, self.system_impedance))
            frequencies.append(ipt.get_frequency())
            
        return frequencies, fwd_gains, vswrs, rev_gains
    #####################################################################################


    #####################################################################################
    #### This method contains the objective function (called target).  It is used by
    #### the optimization algorithm to seek the lowest score given the weights and
    #### geometries that are iterated over.
    #####################################################################################
    def create_optimization_target(self):
        global iteration_count 
        iteration_count = 0

        global sampledResults, sampledL1, sampledL2, sampledL3, sampledD1, sampledD2
        sampledResults = []
        sampledL1 = []
        sampledL2 = []
        sampledL3 = []
        sampledD1 = []
        sampledD2 = []


        # Objective function to calculate a score to be minimized by the optimization algorithms.
        # This procedure is called iteratively by the sciPy minimization procedures to search for 
        # a global minimum of the objective function calculated for variables: Fwd gain, VSWR, Rev Gain
        def target(args):
            l1, l2, l3, d1, d2 = args

            elementDiameterMeters = self.elementDiameter * 0.0254
            foldedDipoleSpacing = self.foldedDipoleSpacing  # distance between driven element and folded element in meters

            # Check that dimensions are valid and no mechanical interference problems
            if l1<elementDiameterMeters or l2<elementDiameterMeters or l3<elementDiameterMeters or d1<=elementDiameterMeters or d2<=foldedDipoleSpacing+elementDiameterMeters:
                logging.warning("TARGET WARNING: Geometry dimensions out of range l1=%6.3f l2=%6.3f, l3=%6.3f, d1=%6.3f d2=%6.3f" % (l1, l2, l3, d1, d2))
                return float('inf')            
            if math.isnan(l1) or math.isnan(l2) or math.isnan(l3) or math.isnan(d1) or math.isnan(d2):
                logging.warning("TARGET WARNING: Geometry dimensions invalid (NaN) l1=%6.3f l2=%6.3f, l3=%6.3f, d1=%6.3f d2=%6.3f" % (l1, l2, l3, d1, d2))
                return float('inf')
            
            try: # simulate yagi and return the parameters
                freqs, gains, vswrs, rev_gains = self.get_gain_swr_range(l1, l2, l3, d1, d2, 
                                                                         start=self.start_frq_opt, stop=self.stop_frq_opt, 
                                                                         step=self.step_frq_opt)
            except RuntimeError as err: # Usually a geometry error occured if we get here
                logging.warning("TARGET WARNING: Exception while calculating optimization score: %s" % (err))
                logging.warning("l1=%6.3f l2=%6.3f, l3=%6.3f, d1=%6.3f d2=%6.3f diameter=%6.3f" % (l1, l2, l3, d1, d2, elementDiameterMeters))
                return float('inf')
            

            # Calculate a score for VSWR, Fwd Gain and Rev Gain over the optimization freq range
            result = 0
            vswr_score = 0
            gains_score = 0
            rev_gains_score = 0
            for gain in gains:
                gains_score += gain     # gain values are in dBi
            for vswr in vswrs:
                #if vswr >= 1.8:         # invoke a severe penalty when VSWR 
                #    vswr = np.exp(vswr) # exceeds 1.8 within the freq ranges :)
                vswr_score += vswr
            for rev_gain in rev_gains:
                if rev_gain > -30:   # shoot for f/b of around 30 (no improvement to score if less)
                    rev_gains_score += rev_gain
                #rev_gains_score += rev_gain

            # Calculate an objective function to be minimized based on the scores
            #    * VSWR should minimal, gains maximal, front-to-back minimal
            #    * "result" variable decreases as these criteria are met
            #    * Note that each contributer is weighted so that each carries appr. equal contribution
            result = self.vswrWeight*vswr_score - self.fwdGainWeight*gains_score + self.fbRatioWeight*(rev_gains_score - gains_score)

            # experiment with alternative objective function scoring system that goes very negative 
            # when all 3 objectives near target simultaneously
            #  --- Gives worse result, and with 3X the iterations (26,126)
            #num_frqs = (stop_frq_opt-start_frq_opt)/step_frq_opt + 1
            #result = -1/((vswr_score-num_frqs) + (-gains_score + 9.5*num_frqs) + (rev_gains_score + 30*num_frqs))


            # These variables are used to create an optimization path map.
            # Since it is only 3 dimensions, can only select two parameters to map 
            global sampledResults, sampledL1, sampledL2, sampledL3, sampledD1, sampledD2
            sampledL1.append(l1)
            sampledL2.append(l2)
            sampledL3.append(l3)
            sampledD1.append(d1)
            sampledD2.append(d2)
            sampledResults.append(result)

            # Bump the iteration counter and display feedback to user on optimization progress
            global iteration_count 
            iteration_count += 1            
            if vswr>1.8: # print correct vswr if severe penalty applied
                #statusText = ("Iteration:%5d   Fwd Gain=%6.2f  Rev Gain=%6.2f VSWR=%6.2f  Score=%6.1f     args=%5.3f, %5.3f, %5.3f, %5.3f, %5.3f" % (iteration_count, gain, rev_gain, np.log(vswr), min(result, 999.9), l1, l2, l3, d1, d2))
                statusText = ("Iteration:%5d  Fwd Gain=%6.2f  Rev Gain=%6.2f VSWR=%6.2f  Score=%6.1f   args=%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f" % (iteration_count, gain, rev_gain, vswr, min(result, 999.9), l1, l2, l3, d1, d2, foldedDipoleSpacing))            
            else:
                statusText = ("Iteration:%5d  Fwd Gain=%6.2f  Rev Gain=%6.2f VSWR=%6.2f  Score=%6.1f   args=%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f" % (iteration_count, gain, rev_gain, vswr, min(result, 999.9), l1, l2, l3, d1, d2, foldedDipoleSpacing))            
            print(statusText, end="\r")
            self.status_label.setText(statusText)
            QApplication.processEvents()  ########## Give UI opportunity to update -- can use threading instead) ######
            
            return result # return result of objective function to calling procedure
        
        return target # returns the name of the scoring procedure to be used
    #####################################################################################



    #####################################################################################
    #### Perform simulation of the geometry passed in the nec parameter and at the 
    #### frequency passed in the freq_mhz parameter.  Returns the impedance.
    #####################################################################################
    def simulate_and_get_impedance(self, nec, freq_mhz):
        nec.set_frequency(freq_mhz) 
        nec.xq_card(0)
        index = 0
        return nec.get_input_parameters(index).get_impedance()
    #####################################################################################




    #####################################################################################
    #### Draw vertical lines at the start and stop optimization frequencies
    #####################################################################################
    def draw_frequency_ranges(self, ax):
        ax.axvline(x=self.start_frq_opt, color='red', linewidth=1, label="Min Frq")
        ax.axvline(x=self.stop_frq_opt, color='red', linewidth=1, label="Max Frq")
    #####################################################################################




    #####################################################################################
    #### Simulate the antenna and plot the results 
    #####################################################################################
    def show_report(self, l1, l2, l3, d1, d2):   
        Omega = "\u2126"
        ohm = ' %s' % Omega

        #### Plot 2-D elevation/azimuth patterns at the theta/phi cuts that include the max gain 
        def plot_pattern_2d(gains_db, thetas, phis, subplot_num, line = "solid", label = None): 
            max_idx = np.unravel_index(np.argmax(gains_db, axis=None), gains_db.shape) # find indices of the max gain_db element
            max_gain = gains_db[max_idx]
            max_theta = thetas[max_idx[0]]
            max_phi = phis[max_idx[1]]

            ax = plt.subplot(subplot_num, polar=True)
            ax.plot(thetas, gains_db[ :, max_idx[1]], color='b', linewidth=1, linestyle=line, label="Elevation Pattern at "+label)  # el pattern at az of max gain
            ax.plot(phis, gains_db[max_idx[0], : ], color='r', linewidth=1, linestyle=line, label="Azimuth Pattern at "+label) # az pattern at el of max gain

            ax.set_thetalim(-np.pi, np.pi)
            ax.set_thetagrids(np.linspace(-180,  180, 24, endpoint=False), fontsize=6)
            ax.set_theta_zero_location("N") # theta is zero due north
            ax.set_theta_direction(-1) # theta increases going clockwise
            perimeter = math.ceil(max_gain/5) * 5 # set the perimeter to the 5 dB value >= the max gain
            #perimeter = max_gain
            center = perimeter - 40
            ax.set_rlim(center, perimeter) # Set outer limit to max gain with 40 dB range
            #ax.set_rgrids(np.linspace(center, perimeter, 5, endpoint=False), angle=0, fontsize=6) 
            #ax.set_rgrids([perimeter-50, perimeter-40, perimeter-30, perimeter-20, perimeter-15, perimeter-10, perimeter-6, perimeter-3, perimeter], angle=0, fontsize=6) 
            ax.set_rgrids([perimeter-40, perimeter-30, perimeter-20, perimeter-10, perimeter], angle=0, fontsize=6) 
            ax.grid(True)
            ax.set_title("Gain patterns at az/el of max gain", va='bottom', fontsize=10)

            return max_gain, max_theta, max_phi
        

        #### Plot Plot 3-D representation of the antenna pattern
        # TODO: Figure out why i get a pick support warning. I think it is mplcursors related. 
        #       Need to turn it off when hovering over 3D plot?)
        def plot_pattern_3d(gains_db, thetas, phis, subplot_num, scaling = "ARRL"): 

            ### Function used by 3-D pattern plotting code to center surface plot on axes
            ### Adapted from plot-antenna.py authored by Ralf Schlatterbeck
            def scene_ranges (matrix = None, add_ground = False):
                """ Create cubic bounding box to force equal aspect ratio """
                x, y, z = matrix
                min_x = x.min()
                max_x = x.max()
                min_y = y.min()
                max_y = y.max()
                min_z = z.min()
                max_z = z.max()
                if add_ground and min_z > 0:
                    min_z = 0
                max_range = np.array([max_x - min_x, max_y - min_y, max_z - min_z ]).max() / 2.0
                mid_x = (max_x + min_x) / 2
                mid_y = (max_y + min_y) / 2
                mid_z = (max_z + min_z) / 2
                xr = np.array([mid_x - max_range, mid_x + max_range])
                yr = np.array([mid_y - max_range, mid_y + max_range])
                zr = np.array([mid_z - max_range, mid_z + max_range])
                # Avoid that something isn't shown due to rounding errors
                if xr[0] > min_x:
                    xr[0] = min_x
                if xr[1] < max_x:
                    xr[1] = max_x
                if yr[0] > min_y:
                    yr[0] = min_y
                if yr[1] < max_y:
                    yr[1] = max_y
                if zr[0] > min_z:
                    zr[0] = min_z
                if zr[1] < max_z:
                    zr[1] = max_z
                if add_ground and min_z == 0:
                    zr = np.array([min_z, min_z + 2 * max_range])
                return np.array([xr, yr, zr])
            # end def scene_ranges 


            # Get parameters associated with the max gain coordinates
            max_idx = np.unravel_index(np.argmax(gains_db, axis=None), gains_db.shape) # find indices of the max gain_db element
            max_gain = gains_db[max_idx]
            max_theta = thetas[max_idx[0]]
            max_phi = phis[max_idx[1]]            

            if scaling == "Linear":
                gains = 10.0**((gains_db - gains_db.max()) / 10.0) # Linear Scaling
            if scaling == "ARRL":
                gains = (1 / 0.89) ** ((gains_db - gains_db.max()) / 2) # ARRL scaling

            n_phis, n_thetas = np.meshgrid(phis, thetas)

            # convert from polar to cartesian coords.
            X = gains * np.sin(n_thetas) * np.cos(n_phis)
            Y = gains * np.sin(n_thetas) * np.sin(n_phis)
            Z = gains * np.cos(n_thetas)
            V = np.sqrt(X**2 + Y**2 + Z**2)

            # Assign colors corresponding to vector length V
            norm = mpl.colors.Normalize(vmin=V.min().min(), vmax=V.max().max())
            mycolors = mpl.cm.rainbow(norm(V)) 

            # Plot the surface
            ax = plt.subplot(subplot_num, projection='3d')
            rc, cc = gains.shape
            ax.plot_surface(X, Y, Z, rcount=rc, ccount=cc, facecolors=mycolors, shade=False)    

            xr, yr, zr = scene_ranges((X, Y, Z)) # center surface plot on cartesian axes
            ax.set_xlim(xr)
            ax.set_ylim(yr)
            ax.set_zlim(zr)
            ax.set_xlabel('X', fontsize=8)
            ax.set_ylabel('Y', fontsize=8)
            ax.set_zlabel('Z', fontsize=8)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            #ax.format_coord = lambda x, y: '' # suppress display of az/el of plot
            ax.set_title("3-D Pattern - 3-element Yagi",  va='bottom', fontsize=10)  

            return max_gain, max_theta, max_phi




        # Display dimensions of the antenna elements (meters)
        print("L1= %5.3f  L2= %5.3f  L3= %5.3f  D1= %5.3f  D2= %5.3f  Folded Dipole Spc=%5.3f" % (l1, l2, l3, d1, d2, self.foldedDipoleSpacing))
        self.dimensionsLabel.setText("Dimensions (m): L1= %5.3f  L2= %5.3f  L3= %5.3f  D1= %5.3f  D2= %5.3f  Folded Dipole Spc=%5.3f" % 
                                     (l1, l2, l3, d1, d2, self.foldedDipoleSpacing))

        #### Run simulation and get impedance at the design frequency 
        nec = self.geometry_yagi(l1, l2, l3, d1, d2)
        freq_mhz = self.design_freq_mhz
        z = self.simulate_and_get_impedance(nec, freq_mhz) # gets impedance results at freq of interest
        print("Impedance: (%0.2f,%+0.2fj)%s @ %0.2f MHz" % (z.real[0], z.imag[0], Omega, freq_mhz))
        print("VSWR @ %0.1f%s is %5.2f @ %0.2f MHz" % (self.system_impedance,Omega, vswr(z, self.system_impedance), freq_mhz))

        #### Run simulation and get vswrs & gains over the start/stop range of freqs
        nec = self.geometry_yagi(l1, l2, l3, d1, d2)
        freqs, fwd_gains, vswrs, rev_gains = self.get_gain_swr_range(l1, l2, l3, d1, d2, start=self.start_frq, stop=self.stop_frq, step=self.step_frq)
        freqs = np.array(freqs) / 1000000 # In MHz


        # Prepare plotting area
        self.main_figure.figure.clear()  # Clears plots, titles and gridlines, etc. 
        plt.subplots_adjust(wspace=0.5, hspace=0.85)  

        #### Plot gains vs frequency 
        ax = plt.subplot(221)
        ax.plot(freqs, fwd_gains, label="Forward Gain (dBi)")
        ax.plot(freqs, rev_gains, linestyle='dotted', label="Reverse Gain (dBi)")
        self.draw_frequency_ranges(ax)

        ax.set_title("Gain - 3-element Yagi", fontsize=10)
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Boresight Gain (dBi)")

        yMaxLimit = math.ceil(max(fwd_gains)/5) * 5
        ax.set_ylim(yMaxLimit-10, yMaxLimit) # Set upper limit to the 5 dB value >= the max gain with 10 dB range

        majorLocator = mpl.ticker.MultipleLocator(1)
        majorFormatter = mpl.ticker.FormatStrFormatter('%d')
        minorLocator = mpl.ticker.MultipleLocator(0.5)
        minorFormatter = mpl.ticker.FormatStrFormatter('')
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_formatter(minorFormatter)
        ax.yaxis.grid(True, which='minor', color='0.6', linestyle='dotted')
        ax.yaxis.grid(True, which='major', color='0.6', linestyle='solid')
        ax.grid(True)
        

        #### Plot vswrs vs frequency 
        ax = plt.subplot(222)
        ax.plot(freqs, vswrs, label="VSWR")
        self.draw_frequency_ranges(ax)

        ax.set_yscale("log")
        ax.set_title("VSWR - 3-element Yagi", fontsize=10)
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("VSWR")

        ax.set_ylim(1, 6)

        majorLocator = mpl.ticker.MultipleLocator(1)
        majorFormatter = mpl.ticker.FormatStrFormatter('%d')
        minorLocator = mpl.ticker.MultipleLocator(0.5)
        minorFormatter = mpl.ticker.FormatStrFormatter('')
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_minor_formatter(minorFormatter)    
        ax.yaxis.grid(True, which='minor', color='0.6', linestyle='dotted')
        ax.yaxis.grid(True, which='major', color='0.6', linestyle='solid')
        ax.grid(True)

        
        #### Optionally plot optimization map in lieu of F/B in slot (224)
        #### Optimization map is useful to test whether opt bounds are capturing the best point
        optMap = True # "True" plots opimizaiton map.  "False" plots F/B ratio vs frequency
        
        global sampledResults, sampledL3, sampledD2
        if optMap and sampledResults:
            '''
            i=0
            for result in sampledResults:
                if result > -150: sampledResults[i]=100000
                i=i+1
            sampledResults = np.log10(-min(sampledResults) + np.array(sampledResults)+.00001)
            norm = mpl.colors.Normalize(vmin=-5, vmax=5, clip=True)
            #sampledResults[::-1] # sort descending
            print(min(sampledResults), max(sampledResults))
            '''
            i=0
            for result in sampledResults:
                if result > -150: sampledResults[i]=-150
                i=i+1
            #norm = mpl.colors.Normalize(vmin=min(sampledResults), vmax=max(sampledResults))
            norm = mpl.colors.Normalize(vmin=-200, vmax=-150, clip=True)

            cmap = mpl.cm.cool
            ax = plt.subplot(224)
            ax.scatter(sampledL3, sampledD2, c=sampledResults, 
                cmap=cmap, norm=norm, edgecolors='face', alpha=0.25)
            ax.scatter(sampledL3, sampledD1, c=sampledResults, 
                cmap=cmap, norm=norm, edgecolors='face', alpha=0.25)
            ax.scatter(sampledL3, sampledL2, c=sampledResults, 
                cmap=cmap, norm=norm, edgecolors='face', alpha=0.25)
            ax.scatter(sampledL1, sampledL2, c=sampledResults, 
                cmap=cmap, norm=norm, edgecolors='face', alpha=0.25)
            ax.set_title("Optimization path", fontsize=10)
            ax.set_xlabel("Director Length (m)")
            ax.set_ylabel("Director Distance (m)")
            self.main_figure.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            ax.grid(True)
        else:
            #### Plot F/B vs frequency 
            f_to_b = np.array(fwd_gains) - np.array(rev_gains)
            ax = plt.subplot(224)
            ax.plot(freqs, f_to_b, label="F/B Ratio (dB)")
            self.draw_frequency_ranges(ax)

            ax.set_title("F/B - 3-element Yagi", fontsize=10)
            ax.set_xlabel("Frequency (MHz)")
            ax.set_ylabel("F/B Ratio (dB)")

            #yMaxLimit = math.ceil(min(40, math.ceil(max(f_to_b)))/5) * 5 # Set upper limit to the 5 dB value >= the max F/B
            #ax.set_ylim(yMaxLimit-40, yMaxLimit)                         #  but cap at 40 dB with 40 dB range
            ax.set_ylim(0, 40)                    

            majorLocator = mpl.ticker.MultipleLocator(10)
            majorFormatter = mpl.ticker.FormatStrFormatter('%d')
            minorLocator = mpl.ticker.MultipleLocator(5)
            minorFormatter = mpl.ticker.FormatStrFormatter('')
            ax.yaxis.set_major_locator(majorLocator)
            ax.yaxis.set_major_formatter(majorFormatter)
            ax.yaxis.set_minor_locator(minorLocator)
            ax.yaxis.set_minor_formatter(minorFormatter)
            ax.yaxis.grid(True, which='minor', color='0.6', linestyle='dotted')
            ax.yaxis.grid(True, which='major', color='0.6', linestyle='solid')
            ax.grid(True)
        ##############################################



        ##############################################
        # repeat simulation and get full spherical radiation pattern at design frequency
        plot_3d = True # "True" selects 3-D plot.  "False" selects 2-D plot for subplot (223)        

        geo_opt = self.geometry_yagi(l1, l2, l3, d1, d2)
        geo_opt.set_frequency(self.design_freq_mhz)
        if plot_3d:
            geo_opt.radiation_pattern(thetas=Range(0, 184, count=46), phis=Range(0,364,count=91), average=0) # python crashes if both counts>1 and average<>0 
        else:
            geo_opt.radiation_pattern(thetas=Range(-180, 180, count=180), phis=Range(0,360,count=180), average=0) # python crashes if both counts>1 and average<>0 
        rp = geo_opt.context.get_radiation_pattern(0) #get the radiation_pattern

        gains_db = rp.get_gain() # Is an array of theta,phi -> gain (in dBi).
        thetas = rp.get_theta_angles()
        phis = rp.get_phi_angles()

        # convert angles from degrees to radians
        thetas = np.deg2rad(thetas)
        phis = np.deg2rad(phis)
     
        if plot_3d:
            #### Plot 3-D representation of the antenna pattern 
            # TODO: Figure out why i get a pick support warning. I think it is mplcursors related. 
            #       Need to turn it off when hovering over 3D plot?)
            max_gain, max_theta, max_phi = plot_pattern_3d(gains_db, thetas, phis, subplot_num=223, scaling = "ARRL")

            print ("Maximal gain is %0.2f dBi" % (max_gain))
            print ("   at an elevation angle of %0.1f degrees" % (max_theta * 180.0 / np.pi))
            print ("   at an azimuthal angle of %0.1f degrees" % (max_phi * 180.0 / np.pi))

            self.performanceLabel.setText("Max Gain: %0.2f dBi   VSWR: %5.2f  Impedance: (%0.2f,%+0.2fj)%s @ %0.2f MHz" % 
                                        (max_gain, vswr(z, self.system_impedance), z.real[0], z.imag[0], Omega, freq_mhz))


        else:
            #### Plot 2-D representation of the antenna pattern 
            max_gain, max_theta, max_phi = plot_pattern_2d(gains_db, thetas, phis, 
                                                           subplot_num=223, label="Design Frq")

            print ("Maximal gain is %0.2f dBi" % (max_gain))
            print ("   at an elevation angle of %0.1f degrees" % (max_theta * 180.0 / np.pi))
            print ("   at an azimuthal angle of %0.1f degrees" % (max_phi * 180.0 / np.pi))

            self.performanceLabel.setText("Max Gain: %0.2f dBi   VSWR: %5.2f  Impedance: (%0.2f,%+0.2fj)%s @ %0.2f MHz" % 
                                        (max_gain, vswr(z, self.system_impedance), z.real[0], z.imag[0], Omega, freq_mhz))


            ''' These additional patterns are interesting but clutter the plot. So make them optional
            '''
            plot_2d_more_freqs = False
            if plot_2d_more_freqs == True:

                # Analyze and plot another pattern at the start optimization freq
                geo_opt = self.geometry_yagi(l1, l2, l3, d1, d2)
                geo_opt.set_frequency(self.start_frq_opt)
                geo_opt.radiation_pattern(thetas=Range(-180, 180, count=180), phis=Range(0,360,count=180), average=0) # python crashes if both counts>1 and average<>0 
                rp = geo_opt.context.get_radiation_pattern(0) #get the radiation_pattern

                gains_db = rp.get_gain() # Is an array of theta,phi -> gain (in dBi).
                thetas = rp.get_theta_angles() * np.pi / 180.0
                phis = rp.get_phi_angles() * np.pi / 180.0

                plot_pattern_2d(gains_db, thetas, phis, subplot_num=223, 
                                line = "dotted", label = "Min Frq")
                

                # Analyze and plot another pattern at the stop optimization freq
                geo_opt = self.geometry_yagi(l1, l2, l3, d1, d2)
                geo_opt.set_frequency(self.stop_frq_opt)
                geo_opt.radiation_pattern(thetas=Range(-180, 180, count=180), phis=Range(0,360,count=180), average=0) # python crashes if both counts>1 and average<>0 
                rp = geo_opt.context.get_radiation_pattern(0) #get the radiation_pattern

                gains_db = rp.get_gain() # Is an array of theta,phi -> gain (in dBi).
                thetas = rp.get_theta_angles() * np.pi / 180.0
                phis = rp.get_phi_angles() * np.pi / 180.0

                plot_pattern_2d(gains_db, thetas, phis, subplot_num=223, 
                                line = "dotted", label = "Max Frq")
                

                # Analyze and plot another pattern at the geometric mean of the optimization freqs
                geo_opt = self.geometry_yagi(l1, l2, l3, d1, d2)
                geo_opt.set_frequency(math.sqrt(self.start_frq_opt*self.stop_frq_opt))
                geo_opt.radiation_pattern(thetas=Range(-180, 180, count=180), phis=Range(0,360,count=180), average=0) # python crashes if both counts>1 and average<>0 
                rp = geo_opt.context.get_radiation_pattern(0) #get the radiation_pattern

                gains_db = rp.get_gain() # Is an array of theta,phi -> gain (in dBi).
                thetas = rp.get_theta_angles() * np.pi / 180.0
                phis = rp.get_phi_angles() * np.pi / 180.0
                
                plot_pattern_2d(gains_db, thetas, phis, subplot_num=223,  
                                line = "dotted", label = "Midband")

        

        # Interactively display the label when cursor hovers over a line
        mplcursors.cursor(hover=mplcursors.HoverMode.Transient).connect(
           "add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))

        ax.figure.canvas.draw()   
    #####################################################################################




#######################################
#### Start application and UI Loop ####
#######################################
if __name__ == "__main__":
    app = QApplication(sys.argv)        
    ui = uiMainWindow()
    ui.show()
    sys.exit(app.exec_())
######### End of application ##########

