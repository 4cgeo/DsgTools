# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DsgTools
                                 A QGIS plugin
 Brazilian Army Cartographic Production Tools
                              -------------------
        begin                : 2018-02-16
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Philipe Borba - Cartographic Engineer @ Brazilian Army
        email                : borba.philipe@eb.mil.br
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
import os

# Qt imports
from PyQt4 import QtGui, uic, QtCore
from PyQt4.QtCore import pyqtSlot, Qt, pyqtSignal
from PyQt4.QtGui import QMessageBox, QApplication, QCursor, QFileDialog, QMenu, QHeaderView

#DsgTools imports
from DsgTools.CustomWidgets.listSelector import ListSelector
from DsgTools.Utils.utils import Utils
from DsgTools.dsgEnums import DsgEnums

from qgis.core import QgsMessageLog
import json

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'genericCompactPropertyManagerWidget.ui'))

class GenericCompactPropertyManagerWidget(QtGui.QWidget, FORM_CLASS):
    Add, Remove, Import, Export, Update = range(5)
    def __init__(self, genericDbManager = None, parent = None):
        """
        Constructor
        """
        super(GenericCompactPropertyManagerWidget, self).__init__(parent)
        self.setupUi(self)
        self.changeTooltips('')
        self.genericDbManager = genericDbManager
        self.textDict = {'EarthCoverage':self.tr('Earth Coverage'), 
                    'Customization':self.tr('Customization'), 
                    'Style':self.tr('Style'), 
                    'ValidationConfig':self.tr('Validation'), 
                    'FieldToolBoxConfig':self.tr('Field Toolbox Configuration'),
                    'Permission':self.tr('Permissions'),
                    'AttributeRuleConfig':self.tr('Attribute Rule Configuration'),
                    'SpatialRuleConfig':self.tr('Spatial Rule Configuration')}
    
    def changeTooltips(self, propertyName):
        """
        Changes all buttons' tooltips according to propertyName
        """
        self.createPropertyPushButton.setToolTip(self.tr('Add {0}').format(propertyName))
        self.removePropertyPushButton.setToolTip(self.tr('Remove {0}').format(propertyName))
        self.importPropertyPushButton.setToolTip(self.tr('Import {0}').format(propertyName))
        self.exportPropertyPushButton.setToolTip(self.tr('Export {0}').format(propertyName))
        self.updatePropertyPushButton.setToolTip(self.tr('Update {0}').format(propertyName))
    
    def enableButtons(self, enabled):
        """
        Enables or disables all buttons according to boolean enabled
        """
        self.createPropertyPushButton.setEnabled(enabled)
        self.removePropertyPushButton.setEnabled(enabled)
        self.importPropertyPushButton.setEnabled(enabled)
        self.exportPropertyPushButton.setEnabled(enabled)
        self.updatePropertyPushButton.setEnabled(enabled)
    
    def setGenericDbManager(self, genericDbManager):
        """
        Sets generic dbManager and adjusts tooltips
        """
        self.genericDbManager = genericDbManager
        self.changeTooltips(self.textDict[self.genericDbManager.getManagerType()])


    def getWhoAmI(self):
        return str(self.__class__).split('.')[-1].replace('\'>', '').replace('GenericCompactPropertyManagerWidget','')  

    def getParameterDict(self):
        """
        This is used to get ui state
        Returns a dict in the format:
        {
            'selectedConfig':--name of the selected config--, 
        }
        """
        parameterDict = dict()
        parameterDict['selectedConfig'] = self.propertyComboBox.currentText()
        return parameterDict
    
    def setInterface(self, parameterDict):
        """
        Uses parameterDict to populate interface
        Sets the interface with a dict in the format:
        {
            'selectedConfig':--name of the selected config--, 
        }
        """
        if 'selectedConfig' in parameterDict:
            idx = self.propertyComboBox.findText(parameterDict['selectedConfig'], flags = Qt.MatchFlags)
            if idx != -1:
                self.propertyComboBox.setCurrentIndex(idx)
    
    @pyqtSlot(bool)
    def on_createPropertyPushButton_clicked(self):
        """
        Creates property
        """
        pass
    
    @pyqtSlot(bool)
    def on_removePropertyPushButton_clicked(self):
        """
        Removes property
        """
        pass

    @pyqtSlot(bool)
    def on_importPropertyPushButton_clicked(self):
        """
        Imports a property file into dsgtools_admindb
        """
        fd = QFileDialog()
        widgetType = self.getWhoAmI()
        filename = fd.getOpenFileName(caption=self.captionDict[widgetType],filter=self.filterDict[widgetType])
        if filename == '':
            QMessageBox.warning(self, self.tr('Warning!'), self.tr('Warning! Select a file to import!'))
            return
        try:
            QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            self.genericDbManager.importSetting(filename)
            QApplication.restoreOverrideCursor()
            QMessageBox.information(self, self.tr('Success!'), self.widgetName + self.tr(' successfully imported.'))
        except Exception as e:
            QApplication.restoreOverrideCursor()
            QMessageBox.critical(self, self.tr('Error!'), self.tr('Error! Problem importing ') +self.widgetName + ': '  + ':'.join(e.args))
        self.refresh()
    
    @pyqtSlot(bool)
    def on_exportPropertyPushButton_clicked(self):
        """
        Export selected properties.
        """
        exportPropertyList = self.selectConfig()
        if exportPropertyList == []:
            QMessageBox.warning(self, self.tr('Warning!'), self.tr('Warning! Select a profile to export!'))
            return
        fd = QFileDialog()
        folder = fd.getExistingDirectory(caption = self.tr('Select a folder to output'))
        if folder == '':
            QMessageBox.warning(self, self.tr('Warning!'), self.tr('Warning! Select a output!'))
            return
        edgvVersion = self.genericDbManager.edgvVersion
        try:
            QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            for exportProperty in exportPropertyList:
                self.genericDbManager.exportSetting(exportProperty, edgvVersion, folder)
            QApplication.restoreOverrideCursor()
            QMessageBox.information(self, self.tr('Success!'), self.widgetName + self.tr(' successfully exported.'))
        except Exception as e:
            QApplication.restoreOverrideCursor()
            QMessageBox.critical(self, self.tr('Error!'), self.tr('Error! Problem exporting ') + self.widgetName + ': ' + ':'.join(e.args))

    @pyqtSlot(bool)
    def on_updatePropertyPushButton_clicked(self):
        """
        Removes property
        """
        pass