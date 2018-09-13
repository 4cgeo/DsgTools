# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DsgTools
                                 A QGIS plugin
 Brazilian Army Cartographic Production Tools
                             -------------------
        begin                : 2018-09-13
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

from qgis.PyQt.QtWidgets import QWidget

from DsgTools.gui.CustomWidgets.DatabaseConversionWidgets.SupportedDrivers.postgisWidget import PostgisWidget
from DsgTools.gui.CustomWidgets.DatabaseConversionWidgets.SupportedDrivers.spatialiteWidget import SpatialiteWidget
from DsgTools.gui.CustomWidgets.DatabaseConversionWidgets.SupportedDrivers.shapefileWidget import ShapefileWidget
from DsgTools.gui.CustomWidgets.DatabaseConversionWidgets.SupportedDrivers.geopackageWidget import GeopackageWidget
from DsgTools.core.dsgEnums import DsgEnums

class DatasourceSelectionWidgetFactory(QWidget):
    """
    Class designed to prepare each selection widget to be added to a widget container.
    """

    def getSelectionWidget(source, parent=None):
        """
        Gets selection widget to be returned to user as selectionWidget attribute.
        :param parent: (QWidget) widget parent to newly instantiated selection widget.
        :param source: (DsgEnum.int) driver enum to have its widget produced.
        :return: (QWidget) selection widget for selected driver.
        """
        sourceDict = {
            DsgEnums.NoDriver : lambda : None,
            DsgEnums.PostGIS : lambda : PostgisWidget(parent=parent),
            DsgEnums.NewPostGIS : lambda : PostgisWidget(parent=parent),
            DsgEnums.SpatiaLite : lambda : SpatialiteWidget(parent=parent),
            DsgEnums.NewSpatiaLite : lambda : SpatialiteWidget(parent=parent),
            DsgEnums.Shapefile : lambda : ShapefileWidget(parent=parent),
            DsgEnums.NewShapefile : lambda : ShapefileWidget(parent=parent),
            DsgEnums.Geopackage : lambda : GeopackageWidget(parent=parent),
            DsgEnums.NewGeopackage : lambda : GeopackageWidget(parent=parent)
        }
        return sourceDict[source]()
