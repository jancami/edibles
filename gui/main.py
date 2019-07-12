#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Requirements: Python3, pandas, matplotlib, PyQt5

import sys
# import numpy as np
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt5.QtCore import QAbstractTableModel, Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog
from gui import Ui_MainWindow
from filter import Ui_FilterDialog
from edibles_spectrum import edibles_spectrum as edspec


# Directory spectra are contained in (TODO: sensible way of initialising this)
path2spectra = '../../spectra/'


class PandasModel(QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """
    def __init__(self, data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._data.columns[col]
        return None


class mainwindow_exec():
    """
    Main GUI window: used to interactively view EDIBLES spectra.

    Command line usage: $ python main.py

    Current functionality:
    - Use filter button to select EDIBLES spectrum (from overview file)
    - Use plot button to display matplotlib plot of highlighted spectrum

    TODO:
    - Implement sub-selection sidebar in MainWindow to allow easy viewing
      of multiple/subselection of spectra (i.e. all spectra of a given object)
      (I/O functions for subselection? overplotted spec?, etc)
    - Add FITS header/stellar/etc info into new panel below MPL canvas
    - Implement filter functions for overview table (by starname, etc)
    - Expand overview file with additional parameters? (exptime, S/N, etc)
    - Integrate fitting/science functions into GUI
      (i.e. interactive lambda selection for profile fitting, etc)
    - add more TODO points
    """
    def __init__(self):
        # Initialise the Main window, and filter dialog
        self.window = QMainWindow()
        self.filter = QDialog(parent=self.window)

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.window)
        self.filterui = Ui_FilterDialog()
        self.filterui.setupUi(self.filter)

        # Load overview file from /catalog/DR3_obslist.ext.txt
        self.load_overview()
        # Initialise MPL figure and toolbar
        self.add_mpl()

        self.window.show()

        # Connect Main window buttons to relevant functions
        self.ui.Button_filter.clicked.connect(lambda: self.filter_dialog())
        self.ui.Button_plot.clicked.connect(lambda: self.update_plot())

    def load_overview(self):
        # Load obslist overview into pandas frame
        self.overview = pd.read_csv('../catalog/DR3_obslist_ext.txt',
                                    delim_whitespace=True)
        cols = list(self.overview)
        # move the Filename column to end of list using index, pop and insert
        cols.insert(len(cols), cols.pop(cols.index('Filename')))
        # use ix to reorder
        self.overview = self.overview.ix[:, cols]

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
            idx = self.filterui.FiltertableView.selectionModel().selectedRows()
            filename = self.model.data(self.model.index(idx[0].row(), 4))
            starname = self.model.data(self.model.index(idx[0].row(), 0))

            # Load wav and flux of corresponding spectrum
            wav, flux = edspec(path2spectra+filename).GetSpectrum()

            # Refresh and replot figure
            self.fig.clf()
            self.ax = self.fig.add_subplot(111)
            self.ax.plot(wav, flux)
            self.ax.set_title(starname)
            self.canvas.draw()
        except(AttributeError, IndexError):
            pass
        except(FileNotFoundError):
            print('FileNotFoundError: \"' + filename + '\" not found...')

    def filter_dialog(self):
        # Initialise overview table from pandas dataframe, and show the dialog
        self.model = PandasModel(self.overview)
        self.filterui.FiltertableView.setModel(self.model)
        self.filter.show()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainwindow_exec()
    sys.exit(app.exec_())
