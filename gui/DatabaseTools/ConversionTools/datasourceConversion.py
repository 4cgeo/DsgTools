# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DsgTools
                                 A QGIS plugin
 Brazilian Army Cartographic Production Tools
                             -------------------
        begin                : 2018-09-05
        git sha              : $Format:%H$
        copyright            : (C) 2018 by João P. Esperidião - Cartographic Engineer @ Brazilian Army
        email                : esperidiao.joao@eb.mil.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

from qgis.PyQt import QtWidgets, uic
from qgis.PyQt.QtCore import Qt, pyqtSignal, pyqtSlot
from qgis.PyQt.QtGui import QIcon

from DsgTools.gui.CustomWidgets.BasicInterfaceWidgets.genericDialogLayout import GenericDialogLayout

from functools import partial
import os, json

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'datasourceConversion.ui'))

class DatasourceConversion(QtWidgets.QWizard, FORM_CLASS):
    # enum for column ordering
    COLUMN_COUNT = 8
    InDs, Filter, InEdgv, inCrs, OutDs, OutEdgv, outCrs, ConversionMode = list(range(COLUMN_COUNT))

    def __init__(self, manager, parentMenu, parent=None):
        """
        Class constructor.
        :param manager:
        :param parentMenu: menu to which conversion tool will be listed
        :param parent: (QWidget) widget parent to a datasource conversion GUI widget.
        """
        super(DatasourceConversion, self).__init__()
        self.setupUi(self)
        # self.setTitle(self.tr('Datasource Conversion Wizard'))
        self.manager = manager
        self.parentMenu = parentMenu
        self.parentButton = None
        self.connectToolSignals()
        # initiate widget dicts
        self.outDs = self.getWidgetNameDict(self.datasourceManagementWidgetOut.activeDrivers)
        self.inDs = self.getWidgetNameDict(self.datasourceManagementWidgetIn.activeDrivers)

    def connectToolSignals(self):
        """
        Connects all tool generic signals.
        """
        # if any widget was turned active/inactive
        self.datasourceManagementWidgetIn.activeWidgetsChanged.connect(self.setTableInitialState)
        self.datasourceManagementWidgetOut.activeWidgetsChanged.connect(self.setTableInitialState)
        # if datasource is changed (e.g. user changed his postgis database selection, for instance)
        self.datasourceManagementWidgetIn.datasourceChangedSignal.connect(self.setTableInitialState)
        self.datasourceManagementWidgetOut.datasourceChangedSignal.connect(self.setTableInitialState)
        # conversion mapping start
        self.button(QtWidgets.QWizard.FinishButton).clicked.connect(self.startConversion)

    def disconnectToolSignals(self):
        """
        Connects all tool generic signals.
        """
        # if any widget was turned active/inactive
        self.datasourceManagementWidgetIn.activeWidgetsChanged.disconnect(self.setTableInitialState)
        self.datasourceManagementWidgetOut.activeWidgetsChanged.disconnect(self.setTableInitialState)
        # if datasource is changed (e.g. user changed his postgis database selection, for instance)
        self.datasourceManagementWidgetIn.datasourceChangedSignal.disconnect(self.setTableInitialState)
        self.datasourceManagementWidgetOut.datasourceChangedSignal.disconnect(self.setTableInitialState)
        self.button(QtWidgets.QWizard.FinishButton).clicked.disconnect(self.startConversion)

    def getWidgetNameDict(self, d):
        """
        Gets the name translated into widget dict from a given dict.
        :param d: (dict) dictionary  - { (str)driver : (QWidget)widget } - to have its widgets translated.
        :return: (dict) translated dict - { (str)datasource_name : (QWidget)widget }.
        """
        returnDict = dict()
        for k in d:
            for w in d[k]:
                returnDict[w.groupBox.title()] = w
        return returnDict

    def resetTable(self, enabled=False):
        """
        Resets table to initial state.
        :param enabled: (bool) indicates whether table will be enabled.
        """
        # clear possible content in it
        self.tableWidget.setEnabled(enabled)
        self.tableWidget.setColumnCount(DatasourceConversion.COLUMN_COUNT)
        # if a signal issue comes along, maybe cell shoul've been systematically removed
        # and signals disconnect for each row along. So that's where shit might be re-coded (=
        self.tableWidget.setRowCount(0)
        # set policy to make cell size adjust to content
        self.tableWidget.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        # map header to its enum
        headerDict = {
            DatasourceConversion.InDs : self.tr("Input"),
            DatasourceConversion.Filter : self.tr("Filters"),
            DatasourceConversion.inCrs : self.tr("In: CRS"),
            DatasourceConversion.InEdgv : self.tr("In: EDGV Version"),
            DatasourceConversion.OutDs : self.tr("Output"),
            DatasourceConversion.outCrs : self.tr("Out: CRS"),
            DatasourceConversion.OutEdgv : self.tr("Out: EDGV Version"),
            DatasourceConversion.ConversionMode : self.tr("Conversion Mode")
        }
        # make the order always follow as presented at enum
        self.tableWidget.setHorizontalHeaderLabels(list([headerDict[i] for i in range(DatasourceConversion.COLUMN_COUNT)]))

    def getRowContents(self, row):
        """
        Retrieves all filled row info.
        :param row: (int) row index to have its output columns populated.
        :return: (list-of-object) a list of all information found on a given row (may be a widget or string).
        """
        # input datasource, input and output EDGV versions are always a text
        inDs = self.tableWidget.item(row, DatasourceConversion.InDs).text()
        inEdgv = self.tableWidget.item(row, DatasourceConversion.InEdgv).text()
        outEdgvItem = self.tableWidget.item(row, DatasourceConversion.OutEdgv)
        outEdgv = outEdgvItem.text() if outEdgvItem else ''
        conversionMode = self.tableWidget.cellWidget(row, DatasourceConversion.ConversionMode).currentText()
        # output datasource, input and output SRC and filter status are always a QWidget
        outDs = self.tableWidget.cellWidget(row, DatasourceConversion.OutDs)
        # however, output info might not yet have been filled
        try:
            _filter = self.tableWidget.cellWidget(row, DatasourceConversion.Filter)
        except:
            _filter = None
        try:
            inCrs = self.tableWidget.cellWidget(row, DatasourceConversion.inCrs)
        except:
            inCrs = None
        try:
            outCrs = self.tableWidget.cellWidget(row, DatasourceConversion.outCrs)
        except:
            outCrs = None
        return [inDs, _filter, inEdgv, inCrs, outDs, outEdgv, outCrs, conversionMode]

    def clearOutDsInforRow(self, row):
        """
        Clears output information for a given row.
        :param row: (int) row to have its output info cleared.
        """
        # add an empty item to it
        self.addItemToTable(col=DatasourceConversion.OutEdgv, row=row, isEditable=False)
        outCrsPb = self.getRowContents(row=row)[6]
        if outCrsPb:
            outCrsPb.setEnabled(False)

    def fillOutDsInfoRow(self, row):
        """
        Fills out row with output info for each output column. In here, ouput SRC and EDGV are filled. 
        :param row: (int) row index to have its output columns populated.
        :return: (list-of-object) return a list containing (str) output EDGV version and (QPushButton) output SRC.
        """
        # clear current content, if any
        self.clearOutDsInforRow(row=row)
        # get only outDs widget
        outDs = self.getRowContents(row=row)[DatasourceConversion.OutDs]
        # widget dict keys are defined as group title, which is part of outDs current text
        if outDs:
            groupTitle = outDs.currentText().split(':')[0]
            crsIcon = QIcon(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'icons', 'CRS_qgis.svg'))
        else:
            return []
        if groupTitle in self.outDs:
            widget = self.outDs[groupTitle]
            # only fills line if dictionary is a controlled widget
            # new push button for SRC
            outCrs = QtWidgets.QPushButton()
            outCrs.setIcon(crsIcon)
            # get new text item to add output datasource
            edgvOut = widget.connWidget.getDatasourceEdgvVersion()
            itemEdgvOut = QtWidgets.QTableWidgetItem()
            itemEdgvOut.setText(edgvOut)
            itemEdgvOut.setFlags(Qt.ItemIsEditable) # not editable
            # add both to table
            self.tableWidget.setCellWidget(row, DatasourceConversion.outCrs, outCrs)
            self.tableWidget.setItem(row, DatasourceConversion.OutEdgv, itemEdgvOut)
            return [edgvOut, outCrs]
        else:
            # if is not controlled, clear line
            return []

    def getTableItem(self, text='', isEditable=True):
        """
        Gets an item to be added to the table that may be set to not be editable.
        :param text: (str) name to be exposed on table cell.
        :param isEditable: (bool) boolean indicating whether cell content should be editable.
        :return: (QTableWidgetItem) item to be added as a table cell.
        """
        item = QtWidgets.QTableWidgetItem()
        item.setText(text)
        if not isEditable:
            item.setFlags(Qt.ItemIsEditable) # not editable
        return item

    def addItemToTable(self, col, row, text='', isEditable=True):
        """
        Adds an item to the mapping table into a given column and row.
        :param col: (int) column containing new item.
        :param row: (int) row containing new item.
        :param text: (str) name to be exposed on table cell.
        :param isEditable: (bool) boolean indicating whether cell content should be editable.
        :return: (QTableWidgetItem) item added.
        """
        newItem = self.getTableItem(text=text, isEditable=isEditable)
        self.tableWidget.setItem(row, col, newItem)
        return newItem

    def prepareRowFilterDialog(self, row):
        """
        Prepares filter dialog for current dataset in a given row.
        :param row: (int) row containing dataset information.
        """
        # get row information
        inDs, _filter, inEdgv, inCrs, outDs, outEdgv, outCrs, conversionMode = self.getRowContents(row=row)
        # retrieve input widget
        inWidget = self.inDs[inDs.split(':')[0]]
        # instantiate a new filter dialog
        filterDlg = GenericDialogLayout()
        # prepare its own GUI
        filterDlg.hideButtons() # hide Ok and Cancel
        title = '{0}: {2} ({1})'.format(inWidget.groupBox.title(), inWidget.connWidget.getDatasourcePath(), \
                                     inWidget.connWidget.getDatasourceConnectionName())
        # set dialog title to current datasource path
        filterDlg.setWindowTitle(title)
        # get layers dict
        layers = inWidget.connWidget.getLayersDict()
        # control dict for each new checkbox added
        checkBoxes = dict()
        for layerName, featCount in layers.items():
            # add a new checkbox widget to layout for each layer found
            checkBoxes[layerName] = QtWidgets.QCheckBox()
            msg = self.tr('{0} ({1} features)') if featCount > 1 else self.tr('{0} ({1} feature)')
            checkBoxes[layerName].setText(msg.format(layerName, featCount))
            if not inWidget.filters['layer'] or (inWidget.filters['layer'] and layerName in inWidget.filters['layer']):
                # in case no filters are added or if layer is among the filtered ones, set it checked
                checkBoxes[layerName].setChecked(True)
            # those are only for confirmation, so it should be disabled
            checkBoxes[layerName].setEnabled(False)
            filterDlg.layout.addWidget(checkBoxes[layerName])
        # connect filter pushbutton signal to newly created dialog
        openDialog = lambda : filterDlg.exec_()
        _filter.clicked.connect(openDialog)

    def setTableInitialState(self):
        """
        Sets the mapping table to its initial (populated) state. Each row is composed by input datasource name, widget with
        all output datasource and conversion mode (whether -9999 would be set to all non-null restriction not respected).
        """
        self.resetTable(enabled=True)
        # output ds dict/list
        self.outDs = self.getWidgetNameDict(self.datasourceManagementWidgetOut.activeDrivers)
        # create a 'function' to get datasource exposing name and create the output list
        getNameAlias = lambda widget : '{0}: {1}'.format(widget.groupBox.title(), widget.getDatasourceConnectionName())
        outDsList = [self.tr('Select a datasource')] + list(map(getNameAlias, self.outDs.values()))
        # input ds dict/list
        self.inDs = self.getWidgetNameDict(self.datasourceManagementWidgetIn.activeDrivers)
        inDsList = list(self.inDs.values())
        # set the table rows # the same as the # of input ds
        self.tableWidget.setRowCount(len(inDsList))
        # prepare widgets control dict
        outWidgets = dict()
        filterIcon = QIcon(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'icons', 'filter.png'))
        crsIcon = QIcon(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'icons', 'CRS_qgis.svg'))
        for idx, w in enumerate(inDsList):
            # initiate widgets map for current row
            outWidgets[idx] = dict()
            # create the item containing current loop's input ds
            t = '{0}: {1}'.format(w.groupBox.title(), w.getDatasourceConnectionName())
            self.addItemToTable(col=DatasourceConversion.InDs, row=idx, text=t, isEditable=False)
            # populate edgv versions column
            # input is always text
            t = w.connWidget.getDatasourceEdgvVersion()
            self.addItemToTable(col=DatasourceConversion.InEdgv, row=idx, text=t, isEditable=False)
            # create filter push button
            outWidgets[idx][DatasourceConversion.Filter] = QtWidgets.QPushButton()
            outWidgets[idx][DatasourceConversion.Filter].setIcon(filterIcon)
            # create input SRC push button
            outWidgets[idx][DatasourceConversion.inCrs] = QtWidgets.QPushButton()
            outWidgets[idx][DatasourceConversion.inCrs].setIcon(crsIcon)
            # create the combobox containing all output ds
            outWidgets[idx][DatasourceConversion.OutDs] = QtWidgets.QComboBox()
            outWidgets[idx][DatasourceConversion.OutDs].addItems(outDsList)
            # create combobox containing conversion mode options
            outWidgets[idx][DatasourceConversion.ConversionMode] = QtWidgets.QComboBox()
            outWidgets[idx][DatasourceConversion.ConversionMode].addItems(['Mode 1', 'Mode 2'])
            # set each widget to their column
            self.tableWidget.setCellWidget(idx, DatasourceConversion.OutDs, outWidgets[idx][DatasourceConversion.OutDs])
            self.tableWidget.setCellWidget(idx, DatasourceConversion.Filter, outWidgets[idx][DatasourceConversion.Filter])
            self.tableWidget.setCellWidget(idx, DatasourceConversion.inCrs, outWidgets[idx][DatasourceConversion.inCrs])
            self.tableWidget.setCellWidget(idx, DatasourceConversion.ConversionMode, outWidgets[idx][DatasourceConversion.ConversionMode])
            # start filter widget
            self.prepareRowFilterDialog(row=idx)
            # set table output information population to its widget
            outWidgets[idx][DatasourceConversion.OutDs].currentIndexChanged.connect(partial(self.fillOutDsInfoRow, row=idx))
        # resize to contents
        self.tableWidget.resizeColumnsToContents()

    def createConversionMap(self):
        """
        Creates the conversion map JSON based on table widget contents.
        :return: (dict) the conversion map. (SPECIFY FORMAT!)
        """
        # initiate conversion map dict
        conversionMap = dict()
        # get the row count from table widget and iterate over it
        for row in range(self.tableWidget.rowCount()):
            # get row contents
            inDs, _filter, inEdgv, inCrs, outDs, outEdgv, outCrs, conversionMode = self.getRowContents(row=row)
            # initiate this row's mapping dict and fill it
            rowMapping = dict()
            title = inDs.split(':')[0] # get input group box title (widget's dict key)
            inputDatasourcePath =  self.inDs[title].connWidget.getDatasourcePath()
            inputFilteredLayers = self.inDs[title].filters
            rowMapping['filter'] = inputFilteredLayers # TEMPORARY - CHANGE FOR ACTUAL FILTER RETRIEVING METHOD LATER
            # replace group title to datasource path
            title = outDs.currentText().split(':')[0] # get output group box title (widget's dict key)
            # retrieve widget's datasource path
            outputDsPath = self.outDs[title].connWidget.getDatasourcePath()
            rowMapping['outDs'] = outputDsPath # still to decide what to fill up in here
            # parameter indicating whether it is a new datasource
            rowMapping['createDb'] = str(self.tr('new') in outDs.currentText())
            rowMapping['conversionMode'] = conversionMode
            # it is possible for the same dataset to be chosen for different outputs, in order to prevent instantiating it
            # more than once, map it all to the same dict entry and control layer/feature flux through filter entry
            if inputDatasourcePath not in conversionMap:
                # setting the conversion map key to input datasource path will secure that conversions using the same ds
                # as data origin will make it be read only once
                conversionMap[inputDatasourcePath] = [rowMapping]
            else:
                conversionMap[inputDatasourcePath].append(rowMapping)
        return conversionMap

    def startConversion(self):
        """
        Converts the datasets.
        """
        conversionMap = self.createConversionMap()
        convMap = os.path.join(os.path.dirname(__file__), 'conversion_map.json')
        with open(convMap, 'w') as fp:
            json.dump(conversionMap, fp, indent=4, sort_keys=True)

    def initGui(self):
        """
        Instantiate GUI for user, including button shortcut (if necessary) and tool insertion on DSGTools tab on QGIS. 
        """
        callback = lambda : self.exec_() 
        self.manager.addTool(
            text=self.tr('Convert Databases'),
            callback=callback,
            parentMenu=self.parentMenu,
            icon='install.png',
            parentButton=self.parentButton,
            defaultButton=False
        )

    def unload(self):
        """
        Unloads GUI.
        """
        try:
            # disconnect all tool signals
            self.disconnectToolSignals()
        except:
            pass
        # remove every widget added to interface (in and output) and, consequently, disconnect all signals
        for d in [self.datasourceManagementWidgetIn.activeDrivers, self.datasourceManagementWidgetOut.activeDrivers]:
            for driverName, wList in d.items():
                for w in wList:
                    # removeWidget method disconnects all widget signals
                    self.datasourceManagementWidgetIn.removeWidget(w)
        # for last, removes itself
        self.setParent(None)
        del self
