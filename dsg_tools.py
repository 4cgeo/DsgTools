# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DsgTools
                                 A QGIS plugin
 Brazilian Army Cartographic Production Tools
                              -------------------
        begin                : 2014-10-09
        git sha              : $Format:%H$
        copyright            : (C) 2014 by Luiz Andrade - Cartographic Engineer @ Brazilian Army
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
from PyQt4.QtCore import *
from PyQt4.QtGui import *
# Initialize Qt resources from file resources.py
import resources_rc
# Import the code for the dialog
from dsg_tools_dialog import DsgToolsDialog
import os.path


class DsgTools:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'DsgTools_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = DsgToolsDialog()

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&DSG Tools')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'DsgTools')
        self.toolbar.setObjectName(u'DsgTools')
        
        self.dsgTools = None
        self.menuBar = self.iface.mainWindow().menuBar()

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('DsgTools', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the InaSAFE toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action
    
    def addMenu(self, parent, name, title, icon_path):
        child = QMenu(parent)
        child.setObjectName(name)
        child.setTitle(self.tr(title))
        child.setIcon(QIcon(icon_path))
        parent.addMenu(child)
        return child
        
    def createToolButton(self, parent, text):
        button = QToolButton(parent)
        button.setObjectName(text)
        button.setToolButtonStyle(Qt.ToolButtonIconOnly)
        button.setPopupMode(QToolButton.MenuButtonPopup)
        parent.addWidget(button)
        return button

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""
        
        self.dsgTools = QMenu(self.iface.mainWindow())
        self.dsgTools.setObjectName(u'DsgTools')
        self.dsgTools.setTitle(self.tr(u'DSG Tools'))

        self.menuBar.insertMenu(self.iface.firstRightStandardMenu().menuAction(), self.dsgTools)

        #Sub menus
        database = self.addMenu(self.dsgTools, u'database', u'Database Tools',':/plugins/DsgTools/database.png')
        layers = self.addMenu(self.dsgTools, u'layers', u'Layer Tools',':/plugins/DsgTools/layers.png')
        complex = self.addMenu(self.dsgTools, u'complex', u'Complex Tools',':/plugins/DsgTools/complex.png')

        #QToolButtons
        self.databaseButton = self.createToolButton(self.toolbar, u'DatabaseTools')
        self.layerButton = self.createToolButton(self.toolbar, u'LayerTools')
        self.complexButton = self.createToolButton(self.toolbar, u'ComplexTools')

        icon_path = ':/plugins/DsgTools/spatialite.png'
        action = self.add_action(
            icon_path,
            text=self.tr(u'Create Spatialite'),
            callback=self.createSpatialiteDatabase,
            parent=database,
            add_to_menu=False,
            add_to_toolbar=False)
        database.addAction(action)
        self.databaseButton.addAction(action)
        self.databaseButton.setDefaultAction(action)
              
        icon_path = ':/plugins/DsgTools/postgis.png'
        action = self.add_action(
            icon_path,
            text=self.tr(u'Create PostGIS'),
            callback=self.createPostGISDatabase,
            parent=database,
            add_to_menu=False,
            add_to_toolbar=False)
        database.addAction(action)
        self.databaseButton.addAction(action)

        icon_path = ':/plugins/DsgTools/category.png'
        action = self.add_action(
            icon_path,
            text=self.tr(u'Load by Category'),
            callback=self.loadByCategory,
            parent=database,
            add_to_menu=False,
            add_to_toolbar=False)
        layers.addAction(action)
        self.layerButton.addAction(action)
        self.layerButton.setDefaultAction(action)
 
        icon_path = ':/plugins/DsgTools/class.png'
        action = self.add_action(
            icon_path,
            text=self.tr(u'Load by Class'),
            callback=self.loadByClass,
            parent=database,
            add_to_menu=False,
            add_to_toolbar=False)
        layers.addAction(action)
        self.layerButton.addAction(action)
 
        icon_path = ':/plugins/DsgTools/frame.png'
        action = self.add_action(
            icon_path,
            text=self.tr(u'Create Frame'),
            callback=self.createFrame,
            parent=database,
            add_to_menu=False,
            add_to_toolbar=False)
        layers.addAction(action)
        self.layerButton.addAction(action)

        icon_path = ':/plugins/DsgTools/complex.png'
        action = self.add_action(
            icon_path,
            text=self.tr(u'Create Complex'),
            callback=self.createComplex,
            parent=database,
            add_to_menu=False,
            add_to_toolbar=False)
        complex.addAction(action)
        self.complexButton.addAction(action)
        self.complexButton.setDefaultAction(action)

        icon_path = ':/plugins/DsgTools/inspect.png'
        action = self.add_action(
            icon_path,
            text=self.tr(u'Inspect Complex'),
            callback=self.inspectComplex,
            parent=database,
            add_to_menu=False,
            add_to_toolbar=False)
        complex.addAction(action)
        self.complexButton.addAction(action)
        
    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&DSG Tools'),
                action)
            self.iface.removeToolBarIcon(action)
            
        if self.dsgTools is not None:
            self.menuBar.removeAction(self.dsgTools.menuAction())

    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass
        
    def createSpatialiteDatabase(self):
        self.databaseButton.setDefaultAction(self.toolbar.sender())
        pass

    def createPostGISDatabase(self):
        self.databaseButton.setDefaultAction(self.toolbar.sender())
        pass

    def loadByCategory(self):
        self.layerButton.setDefaultAction(self.toolbar.sender())
        pass

    def loadByClass(self):
        self.layerButton.setDefaultAction(self.toolbar.sender())
        pass

    def createFrame(self):
        self.layerButton.setDefaultAction(self.toolbar.sender())
        pass

    def createComplex(self):
        self.complexButton.setDefaultAction(self.toolbar.sender())
        pass

    def inspectComplex(self):
        self.complexButton.setDefaultAction(self.toolbar.sender())
        pass

