# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DsgTools
                                 A QGIS plugin
 Brazilian Army Cartographic Production Tools
                              -------------------
        begin                : 2016-09-14
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Luiz Andrade - Cartographic Engineer @ Brazilian Army
        email                : luiz.claudio@dsg.eb.mil.br
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
from qgis.core import QgsMessageLog, QgsVectorLayer, QgsMapLayerRegistry, QgsGeometry, QgsVectorDataProvider, QgsFeatureRequest, QgsExpression, QgsFeature
from DsgTools.ValidationTools.ValidationProcesses.validationProcess import ValidationProcess
import processing, binascii

class SnapGeometriesProcess(ValidationProcess):
    def __init__(self, postgisDb, iface):
        """
        Constructor
        """
        super(self.__class__,self).__init__(postgisDb, iface)
        self.processAlias = self.tr('Snap Geometries')
        
        classesWithElem = self.abstractDb.listClassesWithElementsFromDatabase(useComplex = False, primitiveFilter = ['a', 'l'])
        self.parameters = {'Snap': 1.0, 'MinArea':0.001, 'Classes':classesWithElem.keys()}
        
    def runProcessinAlg(self, layer, tempTableName):
        """
        Runs the actual grass process
        """
        alg = 'grass7:v.clean.advanced'
        
        #creating vector layer
        input = QgsVectorLayer(self.abstractDb.getURI(tempTableName, True).uri(), tempTableName, "postgres")
        crs = input.crs()
        epsg = self.abstractDb.findEPSG()
        crs.createFromId(epsg)
        input.setCrs(crs)
        
        #Adding to registry
        QgsMapLayerRegistry.instance().addMapLayer(input)
        
        #setting tools
        tools = 'snap'
        threshold = -1
        minArea = self.parameters['MinArea']
        snap = self.parameters['Snap']

        #getting table extent (bounding box)
        tableSchema, tableName = self.abstractDb.getTableSchema(tempTableName)        
        (xmin, xmax, ymin, ymax) = self.abstractDb.getTableExtent(tableSchema, tableName)
        extent = '{0},{1},{2},{3}'.format(xmin, xmax, ymin, ymax)
        
        ret = processing.runalg(alg, input, tools, threshold, extent, snap, minArea, None, None)

        #updating original layer
        outputLayer = processing.getObject(ret['output'])
        self.updateOriginalLayer(layer, outputLayer)
          
        #getting error flags
        errorLayer = processing.getObject(ret['error'])
        #removing from registry
        QgsMapLayerRegistry.instance().removeMapLayer(input.id())
        return self.getProcessingErrors(errorLayer)

    def execute(self):
        """
        Reimplementation of the execute method from the parent class
        """
        #abstract method. MUST be reimplemented.
        QgsMessageLog.logMessage(self.tr('Starting ')+self.getName()+self.tr(' Process.'), "DSG Tools Plugin", QgsMessageLog.CRITICAL)
        try:
            self.setStatus(self.tr('Running'), 3) #now I'm running!
            self.abstractDb.deleteProcessFlags(self.getName()) #erase previous flags
            classesWithElem = self.parameters['Classes']
            if len(classesWithElem) == 0:
                self.setStatus(self.tr('Empty database.'), 1) #Finished
                QgsMessageLog.logMessage(self.tr('Empty database.'), "DSG Tools Plugin", QgsMessageLog.CRITICAL)
                return 1
            error = False
            for cl in classesWithElem:
                # preparation
                processTableName, lyr = self.prepareExecution(cl)
                #running the process in the temp table
                result = self.runProcessinAlg(lyr, processTableName)
                self.abstractDb.dropTempTable(processTableName)
                if len(result) > 0:
                    error = True
                    recordList = []
                    for tupple in result:
                        recordList.append((cl,tupple[0],self.tr('Snapping error.'),tupple[1]))
                        self.addClassesToBeDisplayedList(cl) 
                    numberOfProblems = self.addFlag(recordList)
                    QgsMessageLog.logMessage(self.tr('{0} feature(s) of class {1} with snapping errors. Check flags.').format(numberOfProblems, cl), "DSG Tools Plugin", QgsMessageLog.CRITICAL)
                else:
                    QgsMessageLog.logMessage(self.tr('There are no snapping errors on {}.').format(cl), "DSG Tools Plugin", QgsMessageLog.CRITICAL)
            if error:
                self.setStatus(self.tr('There are snapping errors. Check log.'), 4) #Finished with errors
            else:
                self.setStatus(self.tr('There are no snapping errors.'), 1) #Finished
            return 1
        except Exception as e:
            QgsMessageLog.logMessage(':'.join(e.args), "DSG Tools Plugin", QgsMessageLog.CRITICAL)
            self.finishedWithError()
            return 0