# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DsgTools
                                 A QGIS plugin
 Brazilian Army Cartographic Production Tools
                              -------------------
        begin                : 2019-11-14
        git sha              : $Format:%H$
        copyright            : (C) 2019 by João P. Esperidião - Cartographic Engineer @ Brazilian Army
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

from functools import partial

from qgis.core import QgsProject, QgsVectorLayer, QgsMapLayerProxyModel
from qgis.gui import QgsMapLayerComboBox, QgsFieldExpressionWidget
from qgis.PyQt.QtCore import QRegExp
from qgis.PyQt.QtGui import QRegExpValidator
from qgis.PyQt.QtWidgets import (QComboBox,
                                 QLineEdit)
from processing.gui.wrappers import (WidgetWrapper,
                                     DIALOG_STANDARD,
                                     DIALOG_MODELER,
                                     DIALOG_BATCH)

from DsgTools.core.GeometricTools.spatialRelationsHandler import SpatialRelationsHandler
from DsgTools.gui.CustomWidgets.OrderedPropertyWidgets.orderedTableWidget import OrderedTableWidget

class EnforceSpatialRuleWrapper(WidgetWrapper):
    __ATTRIBUTE_MAP_VERSION = 0.1
    def __init__(self, *args, **kwargs):
        super(EnforceSpatialRuleWrapper, self).__init__(*args, **kwargs)

    def ruleNameWidget(self):
        """
        Retrieves the widget for reading/setting rule name.
        :return: (QLineEdit)
        """
        le = QLineEdit()
        le.setPlaceholderText(self.tr("Set a name for this spatial rule..."))
        return le

    def mapLayerComboBox(self):
        """
        Retrieves the configured map layer selection combo box.
        :return: (QgsMapLayerComboBox) configured layer selection widget. 
        """
        cb = QgsMapLayerComboBox()
        cb.setFilters(QgsMapLayerProxyModel.VectorLayer)
        return cb

    def mapLayerModelDialog(self):
        """
        Retrieves widget for map layer selection in a model dialog setup.
        :return: (QLineEdit) map layer setter widget for processing dialog
                 mode.
        """
        le = QLineEdit()
        le.setPlaceholderText(self.tr("Type a vector layer's name..."))
        return le

    def filterExpressionWidget(self):
        """
        Retrieves a new widget for filtering expression setting.
        :return: (QgsFieldExpressionWidget) snap mode selection widget.
        """
        filterWidget = QgsFieldExpressionWidget()
        return filterWidget

    def predicateComboBox(self):
        """
        Retrieves widget for spatial predicate selection.
        :return: (QComboBox) a combo box with all available predicates.
        """
        cb = QComboBox()
        cb.addItems(
            list(SpatialRelationsHandler().availablePredicates().values())
        )
        return cb

    def cardinalityWidget(self):
        """
        Retrieves a widget for cardinality setting.
        :return: (QLineEdit) cardinality widget with its content validation
                 applied.
        """
        le = QLineEdit()
        regex = QRegExp("[0-9\*]\.\.[0-9\*]")
        le.setValidator(QRegExpValidator(regex, le))
        le.setPlaceholderText("1..*")
        return le

    def postAddRowStandard(self, row):
        """
        Sets up widgets to work as expected right after they are added to GUI.
        """
        # in standard GUI, the layer selectors are QgsMapLayerComboBox, and its
        # layer changed signal should be connected to the filter expression
        # widget setup
        for col in [1, 4]:
            mapLayerComboBox = self.panel.itemAt(row, col)
            filterWidget = self.panel.itemAt(row, col + 1)
            mapLayerComboBox.layerChanged.connect(filterWidget.setLayer)
            mapLayerComboBox.layerChanged.connect(
                partial(filterWidget.setExpression, "")
            )
            # first setup is manual though
            vl = mapLayerComboBox.currentLayer()
            if vl:
                filterWidget.setLayer(vl)
        def checkCardinalityAvailability(r):
            predicate = self.panel.getValue(r, 3)
            handler = SpatialRelationsHandler()
            noCardinality = predicate in (
                handler.DISJOINT, handler.NOTEQUALS, handler.NOTINTERSECTS,
                handler.NOTTOUCHES, handler.NOTCROSSES, handler.NOTWITHIN,
                handler.NOTOVERLAPS, handler.NOTCONTAINS
            )
            self.panel.itemAt(row, 6).setEnabled(not noCardinality)
            if noCardinality:
                self.panel.setValue(row, 6, "")
        predicateWidget = self.panel.itemAt(row, 3)
        predicateWidget.currentIndexChanged.connect(
            partial(checkCardinalityAvailability, row)
        )
        # also triggers the action for the first time it is open
        checkCardinalityAvailability(row)

    def postAddRowModeler(self, row):
        """
        Sets up widgets to work as expected right after they are added to GUI.
        """
        def checkLayerBeforeConnect(le, filterExp):
            lName = le.text().strip()
            for layer in QgsProject.instance().mapLayersByName(lName):
                if isinstance(layer, QgsVectorLayer) and layer.name() == lName:
                    filterExp.setLayer(layer)
                    return
            filterExp.setLayer(None)
        for col in [1, 4]:
            le = self.panel.itemAt(row, col)
            filterWidget = self.panel.itemAt(row, col + 1)
            le.editingFinished.connect(
                partial(checkLayerBeforeConnect, le, filterWidget)
            )
        def checkCardinalityAvailability(row):
            predicate = self.panel.getValue(row, 3)
            handler = SpatialRelationsHandler()
            noCardinality = predicate in (
                handler.DISJOINT, handler.NOTEQUALS, handler.NOTINTERSECTS,
                handler.NOTTOUCHES, handler.NOTCROSSES, handler.NOTWITHIN,
                handler.NOTOVERLAPS, handler.NOTCONTAINS
            )
            self.panel.itemAt(row, 6).setText("")
            self.panel.itemAt(row, 6).setEnabled(not noCardinality)   
        predicateWidget = self.panel.itemAt(row, 3)
        predicateWidget.currentIndexChanged.connect(
            partial(checkCardinalityAvailability, row)
        )

    def standardPanel(self):
        """
        Returns the table prepared for the standard Processing GUI.
        :return: (OrderedTableWidget) DSGTools customized table widget.
        """
        otw = OrderedTableWidget(headerMap={
            0 : {
                "header" : self.tr("Rule name"),
                "type" : "widget",
                "widget" : self.ruleNameWidget,
                "setter" : "setText",
                "getter" : "text"
            },
            1 : {
                "header" : self.tr("Layer A"),
                "type" : "widget",
                "widget" : self.mapLayerComboBox,
                "setter" : "setCurrentText",
                "getter" : "currentText"
            },
            2 : {
                "header" : self.tr("Filter A"),
                "type" : "widget",
                "widget" : self.filterExpressionWidget,
                "setter" : "setExpression",
                "getter" : "currentText"
            },
            3 : {
                "header" : self.tr("Predicate"),
                "type" : "widget",
                "widget" : self.predicateComboBox,
                "setter" : "setCurrentIndex",
                "getter" : "currentIndex"
            },
            4 : {
                "header" : self.tr("Layer B"),
                "type" : "widget",
                "widget" : self.mapLayerComboBox,
                "setter" : "setCurrentText",
                "getter" : "currentText"
            },
            5 : {
                "header" : self.tr("Filter B"),
                "type" : "widget",
                "widget" : self.filterExpressionWidget,
                "setter" : "setExpression",
                "getter" : "currentText"
            },
            6 : {
                "header" : self.tr("Cardinality"),
                "type" : "widget",
                "widget" : self.cardinalityWidget,
                "setter" : "setText",
                "getter" : "text"
            }
        })
        otw.setHeaderDoubleClickBehaviour("replicate")
        otw.rowAdded.connect(self.postAddRowStandard)
        return otw

    def batchPanel(self):
        """
        Returns the table prepared for the batch Processing GUI.
        :return: (OrderedTableWidget) DSGTools customized table widget.
        """
        return self.standardPanel()

    def modelerPanel(self):
        """
        Returns the table prepared for the modeler Processing GUI.
        :return: (OrderedTableWidget) DSGTools customized table widget.
        """
        otw = OrderedTableWidget(headerMap={
            0 : {
                "header" : self.tr("Rule name"),
                "type" : "widget",
                "widget" : self.ruleNameWidget,
                "setter" : "setText",
                "getter" : "text"
            },
            1 : {
                "header" : self.tr("Layer A"),
                "type" : "widget",
                "widget" : self.mapLayerModelDialog,
                "setter" : "setText",
                "getter" : "text"
            },
            2 : {
                "header" : self.tr("Filter A"),
                "type" : "widget",
                "widget" : self.filterExpressionWidget,
                "setter" : "setExpression",
                "getter" : "currentText"
            },
            3 : {
                "header" : self.tr("Predicate"),
                "type" : "widget",
                "widget" : self.predicateComboBox,
                "setter" : "setCurrentIndex",
                "getter" : "currentIndex"
            },
            4 : {
                "header" : self.tr("Layer B"),
                "type" : "widget",
                "widget" : self.mapLayerModelDialog,
                "setter" : "setText",
                "getter" : "text"
            },
            5 : {
                "header" : self.tr("Filter B"),
                "type" : "widget",
                "widget" : self.filterExpressionWidget,
                "setter" : "setExpression",
                "getter" : "currentText"
            },
            6 : {
                "header" : self.tr("Cardinality"),
                "type" : "widget",
                "widget" : self.cardinalityWidget,
                "setter" : "setText",
                "getter" : "text"
            }
        })
        otw.setHeaderDoubleClickBehaviour("replicate")
        otw.rowAdded.connect(self.postAddRowModeler)
        return otw

    def createPanel(self):
        return {
            DIALOG_MODELER : self.modelerPanel,
            DIALOG_STANDARD : self.standardPanel,
            DIALOG_BATCH : self.batchPanel
        }[self.dialogType]()
    
    def createWidget(self):
        self.panel = self.createPanel()
        self.panel.showSaveLoadButtons(True)
        self.panel.extension = ".rules"
        self.panel.fileType = self.tr("Set of DSGTools Spatial Rules")
        self.panel.setMetadata({
            "version": self.__ATTRIBUTE_MAP_VERSION
        })
        return self.panel
    
    def parentLayerChanged(self, layer=None):
        pass
    
    def setLayer(self, layer):
        pass
    
    def setValue(self, value):
        """
        Sets back parameters to the GUI. Method reimplementation.
        :param value: (str) value to be set to GUI to retrieve its last state.
        """
        if value is None:
            return
        for valueMap in value:
            self.panel.addRow({
                0 : valueMap["name"],
                1 : valueMap["layer_a"],
                2 : valueMap["filter_a"],
                3 : valueMap["predicate"],
                4 : valueMap["layer_b"],
                5 : valueMap["filter_b"],
                6 : valueMap["cardinality"]
            })

    def readStandardPanel(self):
        """
        Reads widget's contents when process' parameters are set from an 
        algorithm call (e.g. Processing toolbox).
        """
        valueMaplist = list()
        for row in range(self.panel.rowCount()):
            values = dict()
            values["name"] = self.panel.getValue(row, 0).strip() or \
                             self.tr("Spatial Rule #{n}".format(n=row + 1))
            values["layer_a"] = self.panel.getValue(row, 1)
            values["filter_a"] = self.panel.getValue(row, 2)
            values["predicate"] = self.panel.getValue(row, 3)
            values["layer_b"] = self.panel.getValue(row, 4)
            values["filter_b"] = self.panel.getValue(row, 5)
            values["cardinality"] = self.panel.getValue(row, 6) or "1..*"
            valueMaplist.append(values)
        return valueMaplist

    def readModelerPanel(self):
        """
        Reads widget's contents when process' parameters are set from a modeler
        instance.
        """
        return self.readStandardPanel()

    def readBatchPanel(self):
        """
        Reads widget's contents when process' parameters are set from a batch
        processing instance.
        """
        return self.readStandardPanel()

    def value(self):
        """
        Retrieves parameters from current widget. Method reimplementation.
        :return: (dict) value currently set to the GUI.
        """
        return {
            DIALOG_STANDARD : self.readStandardPanel,
            DIALOG_MODELER : self.readModelerPanel,
            DIALOG_BATCH : self.readBatchPanel
        }[self.dialogType]()
    
    def postInitialize(self, wrappers):
        pass
