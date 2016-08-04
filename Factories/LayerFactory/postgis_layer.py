# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DsgTools
                                 A QGIS plugin
 Brazilian Army Cartographic Production Tools
                             -------------------
        begin                : 2015-11-24
        git sha              : $Format:%H$
        copyright            : (C) 2015 by Brazilian Army - Geographic Service Bureau
        email                : suporte.dsgtools@dsg.eb.mil.br
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
from PyQt4.QtCore import pyqtSlot, pyqtSignal
from PyQt4.Qt import QObject

# QGIS imports
from qgis.core import QgsMapLayerRegistry, QgsVectorLayer,QgsDataSourceURI, QgsMessageLog
from qgis.utils import iface

#DsgTools imports
from DsgTools.Factories.LayerFactory.edgv_layer import EDGVLayer

class PostGISLayer(EDGVLayer):
    def __init__(self, iface, abstractDb, codeList, table=None):
        """Constructor."""
        super(PostGISLayer, self).__init__(iface, bstractDb, codeList)
        
        self.provider = 'postgres'
#         self.qmlName = layer.replace('\r','')
        self.setDatabaseConnection()
        self.setDataSource(schema, layer, geomColumn, sql)
        self.geomDict = self.abstractDb.getGeomDict()

    def checkLoaded(self, name, loadedLayers):
        loaded = None
        for ll in loadedLayers:
            if ll.name() == name:
                candidateUri = QgsDataSourceURI(ll.dataProvider().dataSourceUri())
                if self.host == candidateUri.host() and self.database == candidateUri.database() and self.port == int(candidateUri.port()):
                    return ll
        return loaded
    
    def setDatabaseConnection(self):
        self.host = self.abstractDb.db.hostName()
        self.port = self.abstractDb.db.port()
        self.database = self.abstractDb.db.databaseName()
        self.user = self.abstractDb.db.userName()
        self.password = self.abstractDb.db.password()
    
    def buildUri(self):
        self.uri.setConnection(str(self.host),str(self.port), str(self.database), str(self.user), str(self.password))
    
    def setDataSource(self, schema, layer, geomColumn, sql):
        self.uri.setDataSource(schema, layer, geomColumn, sql, 'id')
        self.uri.disableSelectAtId(True)
    
    def getDatabaseGroup(self, groupList):
        dbName = self.abstractDb.getDatabaseName()
        if groupName in groupList:
            return groupList.index(groupName)
        else:
            return self.iface.legendInterface().addGroup(groupName, parent)
    
    def getLyrDict(self, lyrList):
        lyrDict = dict()
        lyrList.sort()
        for lyr in lyrList:
            cat = lyr.split('_')[0]
            if lyr[-1] == 'p':
                if self.tr('Point') not in lyrDict.keys():
                    lyrDict[self.tr('Point')] = dict()
                if cat not in lyrDict[self.tr('Point')].keys():
                    lyrDict[self.tr('Point')][cat] = []
                lyrDict[self.tr('Point')][cat].append(lyr)
            if lyr[-1] == 'l':
                if self.tr('Line') not in lyrDict.keys():
                    lyrDict[self.tr('Line')] = dict()
                if cat not in lyrDict[self.tr('Line')].keys():
                    lyrDict[self.tr('Line')][cat] = []
                lyrDict[self.tr('Line')][cat].append(lyr)
            if lyr[-1] == 'a':
                if self.tr('Area') not in lyrDict.keys():
                    lyrDict[self.tr('Area')] = dict()
                if cat not in lyrDict[self.tr('Area')].keys():
                    lyrDict[self.tr('Area')][cat] = []
                lyrDict[self.tr('Area')][cat].append(lyr)
        return lyrDict
    
    def prepareGroups(self, groupList, parent, lyrDict):
        aux = dict()
        groupDict = dict()
        for geomNode in lyrDict.keys():
            groupDict[geomNode] = dict()
            aux = self.createGroup(groupList, geomNode, parent)
            for catNode in lyrDict[geomNode].keys():
                groupDict[geomNode][catNode] = self.createGroup(groupList, catNode, geomNode)
        return groupDict
    
    def createGroup(self, groupList, groupName, parent):
        if groupName in groupList:
            return groupList.index(groupName) #verificar
        else:
            return self.iface.legendInterface().addGroup(groupName, parent)
    
    def filterLayerList(self, layerList, useInheritance, onlyWithElements):
        filterList = []
        if onlyWithElements:
            lyrsWithElements = self.abstractDb.getLayersWithElements(layerList)
        else:
            lyrsWithElements = layerList
        if useInheritance:
            finalList = self.abstractDb.getLayersFilterByInheritance(lyrsWithElements)
        else:
            finalList = layerList
        return finalList

    def load(self, layerList, useQml = False, uniqueLoad = False, useInheritance = False, stylePath = None, onlyWithElements = False):
        '''
        1. Get loaded layers
        2. Filter layers;
        3. Load domains;
        4. Get Aux Dicts;
        5. Build Groups;
        6. Load Layers;
        '''
        #1. Get Loaded Layers
        loadedLayers = self.iface.legendInterface().layers()
        loadedGroups = self.iface.legendInterface().groups()
        #4. Filter Layers:
        filteredLayerList = self.filterLayerList(layerList, useInheritance, onlyWithElements)
        #2. Load Domains
        dbGroup = self.getDatabaseGroup(loadedGroups)
        domainGroup = self.createGroup(loadedGroups, self.tr("Domains"), dbGroup)
        domLayerDict = self.loadDomains(filteredLayerList, loadedLayers,domainGroup)
        #3. Get Aux dicts
        geomDict = self.abstractDb.getGeomDict()
        domainDict = self.abstractDb.getDomainDict()
        constraintDict = self.abstractDb.getCheckConstraintDict()
        multiColumnsDict = self.abstractDb.getMultiColumnsDict()
        lyrDict = self.getLyrDict(filteredLayerList)
        #4. Build Groups
        groupDict = self.prepareGroups(loadedGroups, dbGroup, filteredLayerList)
        #5. load layers
        for prim in lyrDict.keys():
            for cat in lyrDict[prim].keys():
                self.loadLayer(lyrDict[prim][cat],groupDict[prim][cat], useInheritance, useQml,uniqueLoad,stylePath,geomDict,domainDict,constraintDict,multiColumnsDict)
# 
#             vlayer.loadNamedStyle(vlayerQml, False)
#             attrList = vlayer.pendingFields()
#             for field in attrList:
#                 i = vlayer.fieldNameIndex(field.name())
#                 if vlayer.editorWidgetV2(i) == 'ValueRelation':
#                     groupList = iface.legendInterface().groups()
#                     groupRelationshipList = iface.legendInterface().groupLayerRelationship()
#                     if database not in groupList:
#                         idx = iface.legendInterface().addGroup(database, True,-1)
#                         domainIdGroup = iface.legendInterface().addGroup(self.tr("Dominios"), True, idx)
#                     else:
#                         idx = groupList.index(database)
#                         if "Dominios" not in groupList[idx::]:
#                             domainIdGroup = iface.legendInterface().addGroup(self.tr("Dominios"), True, idx)
#                         else:
#                             domainIdGroup = groupList[idx::].index("Dominios")
#     
#                     valueRelationDict = vlayer.editorWidgetV2Config(i)
#                     domainTableName = valueRelationDict['Layer']
#                     loadedLayers = iface.legendInterface().layers()
#                     domainLoaded = False
#                     for ll in loadedLayers:
#                         if ll.name() == domainTableName:
#                             candidateUri = QgsDataSourceURI(ll.dataProvider().dataSourceUri())
#                             if host == candidateUri.host() and database == candidateUri.database() and port == int(candidateUri.port()):
#                                 domainLoaded = True
#                                 domLayer = ll
#                     if not domainLoaded:
#                         uri = "dbname='%s' host=%s port=%s user='%s' password='%s' key=code table=\"dominios\".\"%s\" sql=" % (database, host, port, user, password, domainTableName)
#                         #TODO Load domain layer into a group
#                         domLayer = iface.addVectorLayer(uri, domainTableName, self.provider)
#                         iface.legendInterface().moveLayer(domLayer, domainIdGroup)
#                     valueRelationDict['Layer'] = domLayer.id()
#                     vlayer.setEditorWidgetV2Config(i,valueRelationDict)
#             self.qmlLoaded.emit()
#         
#         if stylePath:
#             fullPath = self.getStyle(stylePath, self.qmlName)
#             if fullPath:
#                 vlayer.applyNamedStyle(fullPath)
# 
#         iface.legendInterface().moveLayer(vlayer, idSubgrupo)
#             
#         if not vlayer.isValid():
#             QgsMessageLog.logMessage(vlayer.error().summary(), "DSG Tools Plugin", QgsMessageLog.CRITICAL)
#         vlayer = self.createMeasureColumn(vlayer)
#         return vlayer

    def loadLayer(self, lyrName, idSubgrupo, useInheritance, useQml, uniqueLoad,stylePath,geomDict,domainDict,constraintDict,multiColumnsDict):
        if uniqueLoad:
            lyr = self.checkLoaded(lyrName)
            if lyr:
                return lyr
        if useInheritance:
            sql = ''
        schema = geomDict['tablePerspective'][lyrName]['schema']
        geomColumn = geomDict['tablePerspective'][lyrName]['geometryColumn']
        crs =  geomDict['tablePerspective'][lyrName]['srid']
        sql = self.abstractDb.gen.loadLayerFromDatabase(schema+'.'+lyrName)
        self.setDataSource(schema, layer, geomColumn, sql)

        vlayer = iface.addVectorLayer(self.uri.uri(), lyrName, self.provider)
        vlayer.setCrs(crs)
        if useQml:
            qmldir = ''
            try:
                qmldir = self.abstractDb.getQmlDir()
            except Exception as e:
                self.problemOccurred.emit(self.tr('A problem occurred! Check log for details.'))
                QgsMessageLog.logMessage(e.args[0], "DSG Tools Plugin", QgsMessageLog.CRITICAL)
                return None
            vlayerQml = os.path.join(qmldir, self.qmlName+'.qml')
        else:
            vlayer = self.setDomainsAndRestrictions(vlayer, domainDict, constraintDict, multiColumnsDict)
        if stylePath:
            fullPath = self.getStyle(stylePath, self.qmlName)
            if fullPath:
                vlayer.applyNamedStyle(fullPath)
        iface.legendInterface().moveLayer(vlayer, idSubgrupo)   
        if not vlayer.isValid():
            QgsMessageLog.logMessage(vlayer.error().summary(), "DSG Tools Plugin", QgsMessageLog.CRITICAL)
        vlayer = self.createMeasureColumn(vlayer)

    def getDomainsFromDb(self,layerList, loadedLayers):
        domainDict = self.abstractDb.getDomainDict()
        domainList = []
        keys = domainDict.keys()
        for lyr in layerList:
            if lyr in keys:
                for attr in domainDict[lyr]['columns']:
                    dom = domainDict[lyr]['columns']['references']
                    if dom not in domainList:
                        domainList.append(dom)
        return domainList

    def getDomainsToBeLoaded(self, layerList, loadedLayers):
        domains = self.getDomainsFromDb(layerList)
        loadedDomains = []
        for domain in domains:
            domLyr = self.checkLoaded(domain, loadedLayers)
            if domLyr:
                loadedDomains.append(domLyr.name())
        domainsToBeLoaded = []
        for domain in domains:
            if domain not in loadedDomains:
                domainsToBeLoaded.append(domain)
        return domainsToBeLoaded

    def loadDomains(self,layerList, loadedLayers, domainGroup):
        domainsToBeLoaded = self.getDomainsToBeLoaded(layerList, loadedLayers)
        domainsToBeLoaded.sort(reverse=True)
        domLayerDict = dict()
        for domainTableName in domainsToBeLoaded:
            uri = "dbname='%s' host=%s port=%s user='%s' password='%s' key=code table=\"dominios\".\"%s\" sql=" % (self.database, self.host, self.port, self.user, self.password, domainTableName)
            domLayer = iface.addVectorLayer(uri, domainTableName, self.provider)
            domLayerDict[domainTableName]=domLayer
            iface.legendInterface().moveLayer(domLayer, domainIdGroup)
        return domLayerDict

    def getStyleFromDb(self, edgvVersion, className):
        return self.abstractDb.getLyrStyle(edgvVersion,className)
    
    def isLoaded(self,lyr):
        return False

    def setDomainsAndRestrictions(self, lyr, domainDict,constraintDict,multiColumnsDict):
        lyrAttributes = lyr.pendingFields()
        for i in len(lyrAttributes):
            if lyrAttributes[i] == 'id' or 'id_' in lyrAttributes[i]:
                pass
            else:
                if lyr in domainDict.keys():
                    if lyrAttributes[i] in domainDict[lyr]['columns'].keys():
                        refTable = domainDict[lyr]['columns'][attr]['references']
                        refPk = domainDict[lyr]['columns'][attr]['refPk']
                        filter = 
        return lyr
