# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'src/foxsisim_gui/ui/mainwindow.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(549, 540)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 9, 531, 501))
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.doubleSpinBox = QtGui.QDoubleSpinBox(self.tab)
        self.doubleSpinBox.setGeometry(QtCore.QRect(140, 280, 121, 27))
        self.doubleSpinBox.setDecimals(4)
        self.doubleSpinBox.setMinimum(0.0001)
        self.doubleSpinBox.setMaximum(999999.9999)
        self.doubleSpinBox.setSingleStep(0.1)
        self.doubleSpinBox.setObjectName(_fromUtf8("doubleSpinBox"))
        self.doubleSpinBox_2 = QtGui.QDoubleSpinBox(self.tab)
        self.doubleSpinBox_2.setGeometry(QtCore.QRect(140, 310, 121, 27))
        self.doubleSpinBox_2.setDecimals(4)
        self.doubleSpinBox_2.setMinimum(0.0001)
        self.doubleSpinBox_2.setMaximum(999999.9999)
        self.doubleSpinBox_2.setSingleStep(0.1)
        self.doubleSpinBox_2.setObjectName(_fromUtf8("doubleSpinBox_2"))
        self.label = QtGui.QLabel(self.tab)
        self.label.setGeometry(QtCore.QRect(20, 310, 131, 31))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(self.tab)
        self.label_2.setGeometry(QtCore.QRect(20, 280, 121, 31))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.checkBox = QtGui.QCheckBox(self.tab)
        self.checkBox.setGeometry(QtCore.QRect(300, 250, 211, 31))
        self.checkBox.setChecked(True)
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.doubleSpinBox_3 = QtGui.QDoubleSpinBox(self.tab)
        self.doubleSpinBox_3.setGeometry(QtCore.QRect(140, 340, 121, 27))
        self.doubleSpinBox_3.setDecimals(4)
        self.doubleSpinBox_3.setMinimum(-999999.9999)
        self.doubleSpinBox_3.setMaximum(999999.9999)
        self.doubleSpinBox_3.setSingleStep(0.1)
        self.doubleSpinBox_3.setObjectName(_fromUtf8("doubleSpinBox_3"))
        self.label_3 = QtGui.QLabel(self.tab)
        self.label_3.setGeometry(QtCore.QRect(20, 340, 121, 31))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(self.tab)
        self.label_4.setGeometry(QtCore.QRect(20, 370, 121, 31))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_5 = QtGui.QLabel(self.tab)
        self.label_5.setGeometry(QtCore.QRect(20, 400, 121, 31))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.doubleSpinBox_4 = QtGui.QDoubleSpinBox(self.tab)
        self.doubleSpinBox_4.setGeometry(QtCore.QRect(140, 370, 121, 27))
        self.doubleSpinBox_4.setDecimals(4)
        self.doubleSpinBox_4.setMinimum(0.0001)
        self.doubleSpinBox_4.setMaximum(999999.9999)
        self.doubleSpinBox_4.setSingleStep(0.1)
        self.doubleSpinBox_4.setObjectName(_fromUtf8("doubleSpinBox_4"))
        self.doubleSpinBox_5 = QtGui.QDoubleSpinBox(self.tab)
        self.doubleSpinBox_5.setGeometry(QtCore.QRect(140, 400, 121, 27))
        self.doubleSpinBox_5.setDecimals(4)
        self.doubleSpinBox_5.setMinimum(0.0001)
        self.doubleSpinBox_5.setMaximum(999999.9999)
        self.doubleSpinBox_5.setSingleStep(0.1)
        self.doubleSpinBox_5.setObjectName(_fromUtf8("doubleSpinBox_5"))
        self.tableWidget = QtGui.QTableWidget(self.tab)
        self.tableWidget.setGeometry(QtCore.QRect(10, 10, 502, 231))
        self.tableWidget.setAlternatingRowColors(True)
        self.tableWidget.setObjectName(_fromUtf8("tableWidget"))
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        self.tableWidget.horizontalHeader().setDefaultSectionSize(250)
        self.tableWidget.horizontalHeader().setMinimumSectionSize(250)
        self.tableWidget.verticalHeader().setVisible(False)
        self.tableWidget.verticalHeader().setDefaultSectionSize(23)
        self.tableWidget.verticalHeader().setMinimumSectionSize(23)
        self.toolButton_3 = QtGui.QToolButton(self.tab)
        self.toolButton_3.setGeometry(QtCore.QRect(20, 250, 24, 25))
        self.toolButton_3.setObjectName(_fromUtf8("toolButton_3"))
        self.toolButton_4 = QtGui.QToolButton(self.tab)
        self.toolButton_4.setGeometry(QtCore.QRect(50, 250, 24, 25))
        self.toolButton_4.setObjectName(_fromUtf8("toolButton_4"))
        self.label_9 = QtGui.QLabel(self.tab)
        self.label_9.setGeometry(QtCore.QRect(20, 430, 141, 31))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.spinBox_2 = QtGui.QSpinBox(self.tab)
        self.spinBox_2.setGeometry(QtCore.QRect(140, 430, 59, 27))
        self.spinBox_2.setMinimum(1)
        self.spinBox_2.setMaximum(9999)
        self.spinBox_2.setObjectName(_fromUtf8("spinBox_2"))
        self.spinBox_3 = QtGui.QSpinBox(self.tab)
        self.spinBox_3.setGeometry(QtCore.QRect(200, 430, 59, 27))
        self.spinBox_3.setMinimum(1)
        self.spinBox_3.setMaximum(9999)
        self.spinBox_3.setObjectName(_fromUtf8("spinBox_3"))
        self.toolButton_5 = QtGui.QToolButton(self.tab)
        self.toolButton_5.setGeometry(QtCore.QRect(80, 250, 24, 25))
        self.toolButton_5.setArrowType(QtCore.Qt.UpArrow)
        self.toolButton_5.setObjectName(_fromUtf8("toolButton_5"))
        self.toolButton_6 = QtGui.QToolButton(self.tab)
        self.toolButton_6.setGeometry(QtCore.QRect(110, 250, 24, 25))
        self.toolButton_6.setArrowType(QtCore.Qt.DownArrow)
        self.toolButton_6.setObjectName(_fromUtf8("toolButton_6"))
        self.pushButton_5 = QtGui.QPushButton(self.tab)
        self.pushButton_5.setGeometry(QtCore.QRect(300, 280, 201, 27))
        self.pushButton_5.setObjectName(_fromUtf8("pushButton_5"))
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setToolTip(_fromUtf8(""))
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.tableWidget_2 = QtGui.QTableWidget(self.tab_2)
        self.tableWidget_2.setGeometry(QtCore.QRect(10, 10, 500, 171))
        self.tableWidget_2.setAlternatingRowColors(True)
        self.tableWidget_2.setObjectName(_fromUtf8("tableWidget_2"))
        self.tableWidget_2.setColumnCount(6)
        self.tableWidget_2.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        self.tableWidget_2.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget_2.setHorizontalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget_2.setHorizontalHeaderItem(2, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget_2.setHorizontalHeaderItem(3, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget_2.setHorizontalHeaderItem(4, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget_2.setHorizontalHeaderItem(5, item)
        self.tableWidget_2.horizontalHeader().setDefaultSectionSize(83)
        self.tableWidget_2.horizontalHeader().setMinimumSectionSize(83)
        self.tableWidget_2.verticalHeader().setVisible(False)
        self.tableWidget_2.verticalHeader().setDefaultSectionSize(23)
        self.tableWidget_2.verticalHeader().setMinimumSectionSize(23)
        self.toolButton = QtGui.QToolButton(self.tab_2)
        self.toolButton.setGeometry(QtCore.QRect(20, 190, 24, 25))
        self.toolButton.setObjectName(_fromUtf8("toolButton"))
        self.toolButton_2 = QtGui.QToolButton(self.tab_2)
        self.toolButton_2.setGeometry(QtCore.QRect(50, 190, 24, 25))
        self.toolButton_2.setObjectName(_fromUtf8("toolButton_2"))
        self.checkBox_2 = QtGui.QCheckBox(self.tab_2)
        self.checkBox_2.setGeometry(QtCore.QRect(300, 190, 211, 31))
        self.checkBox_2.setChecked(True)
        self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        self.toolButton_7 = QtGui.QToolButton(self.tab_2)
        self.toolButton_7.setGeometry(QtCore.QRect(80, 190, 24, 25))
        self.toolButton_7.setArrowType(QtCore.Qt.UpArrow)
        self.toolButton_7.setObjectName(_fromUtf8("toolButton_7"))
        self.toolButton_8 = QtGui.QToolButton(self.tab_2)
        self.toolButton_8.setGeometry(QtCore.QRect(110, 190, 24, 25))
        self.toolButton_8.setArrowType(QtCore.Qt.DownArrow)
        self.toolButton_8.setObjectName(_fromUtf8("toolButton_8"))
        self.tabWidget.addTab(self.tab_2, _fromUtf8(""))
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.groupBox = QtGui.QGroupBox(self.tab_3)
        self.groupBox.setGeometry(QtCore.QRect(10, 10, 281, 251))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.pushButton_2 = QtGui.QPushButton(self.groupBox)
        self.pushButton_2.setGeometry(QtCore.QRect(10, 210, 261, 31))
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setGeometry(QtCore.QRect(10, 30, 121, 31))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.spinBox = QtGui.QSpinBox(self.groupBox)
        self.spinBox.setGeometry(QtCore.QRect(140, 30, 131, 27))
        self.spinBox.setMaximum(99999999)
        self.spinBox.setObjectName(_fromUtf8("spinBox"))
        self.label_7 = QtGui.QLabel(self.groupBox)
        self.label_7.setGeometry(QtCore.QRect(10, 60, 121, 31))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.label_10 = QtGui.QLabel(self.groupBox)
        self.label_10.setGeometry(QtCore.QRect(10, 90, 121, 31))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.label_8 = QtGui.QLabel(self.groupBox)
        self.label_8.setGeometry(QtCore.QRect(140, 60, 131, 31))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.label_11 = QtGui.QLabel(self.groupBox)
        self.label_11.setGeometry(QtCore.QRect(140, 90, 131, 31))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.pushButton = QtGui.QPushButton(self.groupBox)
        self.pushButton.setGeometry(QtCore.QRect(10, 150, 261, 31))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.pushButton_6 = QtGui.QPushButton(self.groupBox)
        self.pushButton_6.setGeometry(QtCore.QRect(10, 180, 261, 31))
        self.pushButton_6.setObjectName(_fromUtf8("pushButton_6"))
        self.progressBar = QtGui.QProgressBar(self.groupBox)
        self.progressBar.setGeometry(QtCore.QRect(10, 120, 261, 21))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.groupBox_2 = QtGui.QGroupBox(self.tab_3)
        self.groupBox_2.setGeometry(QtCore.QRect(320, 10, 191, 71))
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.pushButton_3 = QtGui.QPushButton(self.groupBox_2)
        self.pushButton_3.setGeometry(QtCore.QRect(10, 30, 171, 31))
        self.pushButton_3.setObjectName(_fromUtf8("pushButton_3"))
        self.groupBox_3 = QtGui.QGroupBox(self.tab_3)
        self.groupBox_3.setGeometry(QtCore.QRect(320, 110, 191, 211))
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.checkBox_3 = QtGui.QCheckBox(self.groupBox_3)
        self.checkBox_3.setGeometry(QtCore.QRect(10, 140, 171, 31))
        self.checkBox_3.setChecked(True)
        self.checkBox_3.setObjectName(_fromUtf8("checkBox_3"))
        self.pushButton_4 = QtGui.QPushButton(self.groupBox_3)
        self.pushButton_4.setGeometry(QtCore.QRect(10, 170, 171, 31))
        self.pushButton_4.setObjectName(_fromUtf8("pushButton_4"))
        self.listWidget = QtGui.QListWidget(self.groupBox_3)
        self.listWidget.setGeometry(QtCore.QRect(10, 30, 171, 111))
        self.listWidget.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.listWidget.setObjectName(_fromUtf8("listWidget"))
        self.tabWidget.addTab(self.tab_3, _fromUtf8(""))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 549, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionAbout_FOXSISIM = QtGui.QAction(MainWindow)
        self.actionAbout_FOXSISIM.setObjectName(_fromUtf8("actionAbout_FOXSISIM"))
        self.actionAbout_FOXSISIM_2 = QtGui.QAction(MainWindow)
        self.actionAbout_FOXSISIM_2.setObjectName(_fromUtf8("actionAbout_FOXSISIM_2"))
        self.menuHelp.addAction(self.actionAbout_FOXSISIM_2)
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "FOXSISIM", None))
        self.label.setToolTip(_translate("MainWindow", "Length of each segment", None))
        self.label.setText(_translate("MainWindow", "Segment Length", None))
        self.label_2.setToolTip(_translate("MainWindow", "Focal length, measured in centimeters from the center of the module", None))
        self.label_2.setText(_translate("MainWindow", "Focal Length", None))
        self.checkBox.setToolTip(_translate("MainWindow", "Autocalculate the angle given a radius (and vice versa)", None))
        self.checkBox.setText(_translate("MainWindow", "Autocalculate Dimensions", None))
        self.label_3.setToolTip(_translate("MainWindow", "Offset of the detector from the focal point. A positive value means the detector is farther from the module.", None))
        self.label_3.setText(_translate("MainWindow", "Detector Offset", None))
        self.label_4.setToolTip(_translate("MainWindow", "Width of the detector", None))
        self.label_4.setText(_translate("MainWindow", "Detector Width", None))
        self.label_5.setToolTip(_translate("MainWindow", "Height of the detector", None))
        self.label_5.setText(_translate("MainWindow", "Detector Height", None))
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Shell Radius", None))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Shell Angle", None))
        self.toolButton_3.setToolTip(_translate("MainWindow", "Add a row to the radii/angles table", None))
        self.toolButton_3.setText(_translate("MainWindow", "+", None))
        self.toolButton_4.setToolTip(_translate("MainWindow", "Remove a row from the radii/angles table", None))
        self.toolButton_4.setText(_translate("MainWindow", "-", None))
        self.label_9.setToolTip(_translate("MainWindow", "Pixel resolution of detector", None))
        self.label_9.setText(_translate("MainWindow", "Detector Reso", None))
        self.toolButton_5.setToolTip(_translate("MainWindow", "Move current row up", None))
        self.toolButton_5.setText(_translate("MainWindow", "...", None))
        self.toolButton_6.setToolTip(_translate("MainWindow", "Move current row down", None))
        self.toolButton_6.setText(_translate("MainWindow", "...", None))
        self.pushButton_5.setToolTip(_translate("MainWindow", "Plot a cross section of the module", None))
        self.pushButton_5.setText(_translate("MainWindow", "Plot Cross Section", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Module", None))
        item = self.tableWidget_2.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Type", None))
        item = self.tableWidget_2.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Center", None))
        item = self.tableWidget_2.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "Normal", None))
        item = self.tableWidget_2.horizontalHeaderItem(3)
        item.setText(_translate("MainWindow", "Width", None))
        item = self.tableWidget_2.horizontalHeaderItem(4)
        item.setText(_translate("MainWindow", "Height", None))
        item = self.tableWidget_2.horizontalHeaderItem(5)
        item.setText(_translate("MainWindow", "Color", None))
        self.toolButton.setToolTip(_translate("MainWindow", "Add row to the sources table", None))
        self.toolButton.setText(_translate("MainWindow", "+", None))
        self.toolButton_2.setToolTip(_translate("MainWindow", "Remove row from sources table", None))
        self.toolButton_2.setText(_translate("MainWindow", "-", None))
        self.checkBox_2.setToolTip(_translate("MainWindow", "Autocalculate as many source settings as is possible", None))
        self.checkBox_2.setText(_translate("MainWindow", "Autocalculate Dimensions", None))
        self.toolButton_7.setToolTip(_translate("MainWindow", "Move current row up", None))
        self.toolButton_7.setText(_translate("MainWindow", "...", None))
        self.toolButton_8.setToolTip(_translate("MainWindow", "Move current row down", None))
        self.toolButton_8.setText(_translate("MainWindow", "...", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Sources", None))
        self.groupBox.setTitle(_translate("MainWindow", "Simulate", None))
        self.pushButton_2.setToolTip(_translate("MainWindow", "Reset session, deleting all simulation variables", None))
        self.pushButton_2.setText(_translate("MainWindow", "Reset", None))
        self.label_6.setToolTip(_translate("MainWindow", "Number of rays to simulate for each source", None))
        self.label_6.setText(_translate("MainWindow", "Rays Per Source", None))
        self.label_7.setToolTip(_translate("MainWindow", "Total rays to be simulated", None))
        self.label_7.setText(_translate("MainWindow", "Rays to Simulate", None))
        self.label_10.setToolTip(_translate("MainWindow", "Total rays simulated thus far", None))
        self.label_10.setText(_translate("MainWindow", "Total Simulated", None))
        self.label_8.setText(_translate("MainWindow", "0", None))
        self.label_11.setText(_translate("MainWindow", "0", None))
        self.pushButton.setToolTip(_translate("MainWindow", "Start simulation. Once pressed, module/detector/source settings cannot be modified until session is reset.", None))
        self.pushButton.setText(_translate("MainWindow", "Simulate", None))
        self.pushButton_6.setToolTip(_translate("MainWindow", "Stop current simulation", None))
        self.pushButton_6.setText(_translate("MainWindow", "Stop", None))
        self.groupBox_2.setTitle(_translate("MainWindow", "Detector Pixel Plot", None))
        self.pushButton_3.setToolTip(_translate("MainWindow", "Plot the image captured by the detector", None))
        self.pushButton_3.setText(_translate("MainWindow", "Plot", None))
        self.groupBox_3.setTitle(_translate("MainWindow", "Scatter Plot", None))
        self.checkBox_3.setToolTip(_translate("MainWindow", "Rays will be colored by how many times they bounce in the module: none (red), once (green), twice (blue), more (black)", None))
        self.checkBox_3.setText(_translate("MainWindow", "Color Bounce Count", None))
        self.pushButton_4.setToolTip(_translate("MainWindow", "Generate a scatter plot of captured rays from one or more sources", None))
        self.pushButton_4.setText(_translate("MainWindow", "Plot", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Simulation", None))
        self.menuHelp.setTitle(_translate("MainWindow", "Help", None))
        self.actionAbout_FOXSISIM.setText(_translate("MainWindow", "About FOXSISIM", None))
        self.actionAbout_FOXSISIM_2.setText(_translate("MainWindow", "About FOXSISIM", None))

