# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'filter.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_FilterDialog(object):
    def setupUi(self, FilterDialog):
        FilterDialog.setObjectName("FilterDialog")
        FilterDialog.resize(900, 757)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(FilterDialog.sizePolicy().hasHeightForWidth())
        FilterDialog.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(FilterDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.FiltertableView = QtWidgets.QTableView(FilterDialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.FiltertableView.sizePolicy().hasHeightForWidth())
        self.FiltertableView.setSizePolicy(sizePolicy)
        self.FiltertableView.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        self.FiltertableView.setEditTriggers(QtWidgets.QAbstractItemView.AnyKeyPressed|QtWidgets.QAbstractItemView.DoubleClicked|QtWidgets.QAbstractItemView.EditKeyPressed|QtWidgets.QAbstractItemView.SelectedClicked)
        self.FiltertableView.setAlternatingRowColors(False)
        self.FiltertableView.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.FiltertableView.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.FiltertableView.setObjectName("FiltertableView")
        self.FiltertableView.horizontalHeader().setSortIndicatorShown(False)
        self.FiltertableView.horizontalHeader().setStretchLastSection(True)
        self.verticalLayout.addWidget(self.FiltertableView)

        self.retranslateUi(FilterDialog)
        QtCore.QMetaObject.connectSlotsByName(FilterDialog)

    def retranslateUi(self, FilterDialog):
        _translate = QtCore.QCoreApplication.translate
        FilterDialog.setWindowTitle(_translate("FilterDialog", "Dialog"))

