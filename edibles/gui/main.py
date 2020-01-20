#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Requirements: Python3, pandas, matplotlib, PyQt5

import sys
import numpy as np
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, qApp, QFileDialog
from edibles.edibles.gui.gui import Ui_MainWindow
from edibles.edibles.functions.edibles_spectrum import EdiblesSpectrum as edspec
from edibles.edibles.gui.models import PandasModel, SelectionModel
from edibles.edibles import EDIBLES_PYTHONDIR, DATADIR, DATARELEASE


class MainWindow(QMainWindow, Ui_MainWindow):
    """
    Main GUI window: used to interactively view EDIBLES spectra.

    Command line usage: $ python main.py

    Alternatively, add the following to .bashrc/bashprofile for command line
    usage anywhere:

    $ alias edigui="python /FULL/PATH/TO/edibles/gui/main.py"

    Alternate command line usage: $ edigui

    Current functionality:
    - Use filter tab to select EDIBLES spectra
    - Use plot tab to display matplotlib plot of highlighted spectra
    - Basic I/O functions for subselection
    - Basic filtering by object and wavelength (value or range)

    TODO:
    - Add FITS header/stellar/etc info into new panel below MPL canvas?
    - More filter functions for overview table?
    - Expand overview file with additional parameters? (exptime, S/N, etc)
    - Integrate fitting/science functions into GUI
      (i.e. interactive lambda selection for profile fitting, etc)
    - add more TODO points
    """
    def __init__(self, parent=None):
        # Initialise the Main window
        QMainWindow.__init__(self, parent=parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.selected_data = []

        # Load overview file from /catalog/DR3_obslist.ext.txt
        self.load_overview()
        # Initialise Filter/Data overview
        self.filtertab(self.overview)
        # Initialise MPL figure and toolbar
        self.add_mpl()

        # Connect Main window buttons to relevant functions
        self.ui.FilterAddButton.clicked.connect(lambda: self.filter_add())
        self.ui.PlotButton.clicked.connect(lambda: self.update_plot())
        self.ui.SearchPushButton.clicked.connect(lambda: self.Objectfilter())
        self.ui.Selectedpushbutton.clicked.connect(lambda: self.Remove_selected())
        self.ui.menuExit.triggered.connect(qApp.quit)
        # TODO Probably make these names shorter...
        self.ui.menuSelectedDataImport.triggered.connect(self.ImportSelData)
        self.ui.menuSelectedDataExport.triggered.connect(self.ExportSelData)

    def ImportSelData(self):
        # retrieve file to load in
        fname = QFileDialog.getOpenFileName(self, 'Open file', './')
        if fname[0]:
            with open(fname[0], 'r') as f:
                # for filename in file
                for filename in f.read().splitlines():
                    if 'ascii' in filename:
                        filename = filename.replace('ascii', 'fits')
                    if filename in self.selected_data:
                        # if filename already in selected data, skip
                        continue
                    # add filename to selected data
                    self.selected_data.append(filename)
                # refresh selected data table
                if len(self.selected_data):
                    self.selectionmodel = SelectionModel(self.selected_data)
                    self.ui.SelectedDataTable.setModel(self.selectionmodel)

    def ExportSelData(self):
        # create file to export selected data to
        fname = QFileDialog.getSaveFileName(self, 'Save file', './')
        if fname[0]:
            with open(fname[0], 'w') as f:
                for row in self.selected_data:
                    if self.ui.actionEnable_ASCII_format.isChecked():
                        f.write(row[:-4] + 'ascii\n')
                    else:
                        f.write(row + '\n')

    def keyPressEvent(self, event):
        # Delete key used to remove highlighted spectra from selected sidebar
        if event.key() == QtCore.Qt.Key_Delete:
            self.Remove_selected()

    def Remove_selected(self):
        # Get highlighted rows
        idx = self.ui.SelectedDataTable.selectionModel().selectedRows()

        # Retrieve filename corresponding to highlighted rows
        # ( TODO probably a nicer way to do this?)
        filenames = []
        for idxxx in idx:
            filenames.append(self.selectionmodel.data(
                             self.selectionmodel.index(idxxx.row(), 0)))

        # Iterate over highlighted files, remove rom selection list
        for filename in filenames:
            if len(self.selected_data):
                if filename in self.selected_data:
                    self.selected_data.remove(filename)
                    continue

        # Regenerate selection model and tableview
        self.selectionmodel = SelectionModel(self.selected_data)
        self.ui.SelectedDataTable.setModel(self.selectionmodel)

    def load_overview(self):
        # Load obslist overview into pandas frame
        self.overview = pd.read_csv(EDIBLES_PYTHONDIR + '/edibles/data/' + DATARELEASE + '_ObsLog.csv')

    def add_mpl(self):
        # Initial MPL figure and toolbar
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self.ui.matplotlib,
                                         coordinates=True)
        # Add MPL canvas and toolbar to Main window
        self.ui.mpl_layout.addWidget(self.toolbar)
        self.ui.mpl_layout.addWidget(self.canvas)
        self.canvas.draw()

    def update_plot(self):
        try:
            # Retrieve highlighted row from tableview
            idx = self.ui.SelectedDataTable.selectionModel().selectedRows()
            self.fig.clf()
            self.ax = self.fig.add_subplot(111)
            for idxxx in idx:
                filename = self.selectionmodel.data(
                           self.selectionmodel.index(idxxx.row(), 0))
                # Load wav and flux of corresponding spectrum
                if self.ui.actionEnable_ASCII_format.isChecked():
                    wav, flux = np.loadtxt(DATADIR+filename[:-4]+'ascii',unpack=True)
                else:
                    df = edspec(filename).getSpectrum()
                    wav = df['wave']
                    flux = df['flux']

                # Refresh and replot figure
                self.ax.plot(wav, flux)
                self.ax.set_xlabel(r'Wavelength ($\AA$)')
                self.ax.set_ylabel(r'Flux')

            self.canvas.draw()
        except(AttributeError, IndexError):
            pass
        except(FileNotFoundError):
            print('FileNotFoundError: \"' + filename + '\" not found...')

    def filtertab(self, input):
        # Initialise overview table from pandas dataframe, and show the dialog
        self.model = PandasModel(input)
        self.ui.FiltertableView.setModel(self.model)

    def filter_add(self):
            idx = self.ui.FiltertableView.selectionModel().selectedRows()
            skiptotal = 0
            for idxxx in idx:
                filename = self.model.data(self.model.index(idxxx.row(), 7))

                if len(self.selected_data):
                    if filename in self.selected_data:
                        skiptotal += 1
                        continue
                self.selected_data.append(filename)
            if len(self.selected_data):
                self.selectionmodel = SelectionModel(self.selected_data)
                self.ui.SelectedDataTable.setModel(self.selectionmodel)

            if skiptotal:
                if skiptotal == len(idx):
                    statustext = str(skiptotal) + ' spectra already added...'
                else:
                    statustext = (str(skiptotal) + ' spectra already added. ' +
                                  str(len(idx)-skiptotal) +
                                  ' new spectra added...')
            elif len(idx):
                statustext = str(len(idx))+' spectra addded...'
            else:
                statustext = ''
            self.ui.statusBar.showMessage(statustext, 3000)

    def Objectfilter(self):
        newdata = self.overview
        obj_in = self.ui.ObjLineEdit.text()
        wav_in = self.ui.WavLineEdit.text()

        if len(obj_in):
            newdata = newdata[newdata['Object'].str.contains(obj_in, na=False,
                                                             case=False)]

        if len(wav_in):
            if '-' in wav_in:
                wav_min, wav_max = wav_in.split('-')
                try:
                    newdata = newdata[((newdata['WaveMin'] < float(wav_min)) &
                                      (newdata['WaveMax'] > float(wav_min))) |
                                      ((newdata['WaveMin'] < float(wav_max)) &
                                      (newdata['WaveMax'] > float(wav_max))) |
                                      ((newdata['WaveMin'] < float(wav_max)) &
                                      (newdata['WaveMax'] > float(wav_min)))]
                except(ValueError):
                    pass
            else:
                try:
                    newdata = newdata[(newdata['WaveMin'] < float(wav_in)) &
                                      (newdata['WaveMax'] > float(wav_in))]
                except(ValueError):
                    pass
        self.filtertab(newdata)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())
