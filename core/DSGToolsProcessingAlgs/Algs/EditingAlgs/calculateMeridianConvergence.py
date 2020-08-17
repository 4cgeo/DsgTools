# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""
from builtins import str
import os
from qgis.PyQt import uic
from qgis.PyQt.QtWidgets import QDialog
from qgis.core import QgsCoordinateReferenceSystem, QgsDistanceArea, QgsCoordinateTransform, QgsPointXY, QgsProject

import math

from PyQt5.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsFeature,
                       QgsField,
                       QgsFields,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterField,
                       Qgis,
                       )
import processing
import requests
import json
from datetime import datetime
from qgis.utils import iface
from qgis.gui import QgsMessageBar
#from qgis.PyQt.QtWidgets import *
#from qgis.core import Qgis
#from qgis.PyQt.QtCore import QVariant

class ConvergenciaProcessingAlgorithm(QgsProcessingAlgorithm):
    """
    Este scripy calcula a CONVERGENCIA MERIDIANA.
    """
    INPUT = 'INPUT'
    #OUTPUT = 'OUTPUT'
    CONV_MER = 'convergencia_meridiana'
    #VAR_ANUAL = 'var_anual'

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ConvergenciaProcessingAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'convergencia_meridiana'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Convergencia Meridiana')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('scripts')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'scripts'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr("Calcula a Convergencia de Meridiano. Crie o campo 'conv_mer' antes de rodar a ferramenta. Produzido pela Seção de Desenvolvimento do 4º CGEO")

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        # We add the input vector features source. It can have any kind of
        # geometry.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Selecione a Camada'),
                [QgsProcessing.TypeVectorPolygon]
            )
        )
        
        self.addParameter(
            QgsProcessingParameterField(
                self.CONV_MER,
                'Selecione o campo para Convergência Meridiana:',
                '',
                self.INPUT))
               
        
    def calculateConvergence2(self, longitude, latitude, a, b):
        """Calculates the meridian convergence
        """
        centralMeridian = int(abs(longitude)/6)*6 + 3
        if longitude < 0:
            centralMeridian = centralMeridian*(-1)

        deltaLong =  centralMeridian - longitude 

        p = 0.0001*( deltaLong*3600 )

        xii = math.sin(math.radians(latitude))*math.pow(10, 4)

        e2 = math.sqrt(a*a - b*b)/b

        c5 = math.pow(math.sin(math.radians(1/3600)), 4)*math.sin(math.radians(latitude))*math.pow(math.cos(math.radians(latitude)), 4)*(2 - math.pow(math.tan(math.radians(latitude)), 2))*math.pow(10, 20)/15

        xiii = math.pow(math.sin(math.radians(1/3600)), 2)*math.sin(math.radians(latitude))*math.pow(math.cos(math.radians(latitude)), 2)*(1 + 3*e2*e2*math.pow(math.cos(math.radians(latitude)), 2) + 2*math.pow(e2, 4)*math.pow(math.cos(math.radians(latitude)), 4))*math.pow(10, 12)/3

        cSeconds = xii*p + xiii*math.pow(p, 3) + c5*math.pow(p, 5)

        c = cSeconds/3600

        return c

    def getSemiMajorAndSemiMinorAxis(self):
        """Obtains the semi major axis and semi minor axis from the used ellipsoid
        """
        currentLayer = iface.activeLayer()
        distanceArea = QgsDistanceArea()
        distanceArea.setEllipsoid(currentLayer.crs().ellipsoidAcronym())
        a = distanceArea.ellipsoidSemiMajor()
        b = distanceArea.ellipsoidSemiMinor()
        
        return (a,b)
 
    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        source = self.parameterAsSource(
            parameters,
            self.INPUT,
            context
        )
        
        conv_mer = self.parameterAsString(
            parameters,
            self.CONV_MER,
            context)
       
        
        # Ler todas as feições selecionadas
        layer = iface.activeLayer()
        selection = layer.selectedFeatures()
        for f in selection:
            geom = f.geometry() #le a geometria da feição
            
            
            #MI = f.attribute('mi_1') #le o atributo da feição
            #print("MI: ", MI)
            #print('Area: ', geom.area())
            #print('Perimeter: ', geom.length())
            coord = str(geom.centroid()) #pega a coordenada do centro e converte em string
            t = coord.split() #Quebra a string nos espaços 
            
            long_1 =t[2] #recupera a longitude na lista 
            long2 = float(long_1[1:]) #retira o ( da string e converte para float
            #verifica se coordenada long esta em W ou E
            if long2 <= 0:
                lon1Hemisphere = "W"
                print(lon1Hemisphere)
                long3 = long_1[2:]
            else:
                lon1Hemisphere = "E"
                print(lon1Hemisphere)
                long3 = long_1[1:]

            lat_1 =t[3] #recupera a latitude na lista 
            lat2 = float(lat_1[:-2]) #retira o ) da string
            #verifica se coordenada lat esta em N ou S
            if lat2 <= 0:
                lat1Hemisphere = "S"
                print(lat1Hemisphere)
                lat3 = lat_1[1:-3]
            else:
                lat1Hemisphere = "N"
                print(lat1Hemisphere)
                lat3 = lat_1[:-3]
            
            (a,b) = self.getSemiMajorAndSemiMinorAxis()

            texto_conv =  str(self.calculateConvergence2(float(long3), float(lat3), a, b))
            
            feedback.pushInfo('Escrevendo os atributos!')
            conv_mer = f.fields().indexFromName("conv_mer")
            layer.startEditing()
            layer.changeAttributeValue(f.id(), conv_mer, texto_conv)
            layer.commitChanges()
          
    iface.messageBar().pushMessage("Calculando a Convergência Meridiana", level=Qgis.Success, duration=20)
        