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

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProject,
                       QgsVectorLayer,
                       QgsProcessingContext,
                       QgsProcessingException,
                       QgsProcessingParameterType,
                       QgsProcessingMultiStepFeedback,
                       QgsProcessingParameterDefinition,
                       QgsProcessingParameterFeatureSink)

from DsgTools.core.DSGToolsProcessingAlgs.algRunner import AlgRunner
from DsgTools.core.GeometricTools.spatialRelationsHandler import SpatialRelationsHandler
from DsgTools.core.DSGToolsProcessingAlgs.Algs.ValidationAlgs.validationAlgorithm import ValidationAlgorithm

class EnforceSpatialRulesAlgorithm(ValidationAlgorithm):
    RULES_SET = "RULES_SET"
    POINT_FLAGS = "POINT_FLAGS"
    LINE_FLAGS = "LINE_FLAGS"
    POLYGON_FLAGS = "POLYGON_FLAGS"
    # a map to canvas layers, reset every time alg is asked to be run: make
    # reusage of layers available for each cycle, reducing re-reading time 
    __layers = dict()
    # the order here matters and must be denial right after and the same as in
    # the wrapper widget
    __predicates = [
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "equals"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "is not equals"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "disjoint"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "is not disjoint"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "intersects"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "does not intersect"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "touches"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "does not touch"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "crosses"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "does not cross"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "within"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "is not within"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "overlaps"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "does not overlap"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "contains"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "does not contain"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "covers"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "does not cover"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "covered by"),
            QCoreApplication.translate("EnforceSpatialRulesAlgorithm", "is not covered by")
        ]
    EQUALS, NOTEQUALS, DISJOINT, NOTDISJOINT, INTERSECTS, NOTINTERSECTS, \
        TOUCHES, NOTTOUCHES, CROSSES, NOTCROSSES, WITHIN, NOTWITHIN, OVERLAPS, \
        NOTOVERLAPS, CONTAINS, NOTCONTAINS, COVERS, NOTCOVERS, COVEREDBY, \
        NOTCOVEREDBY = range(len(__predicates))


    def initAlgorithm(self, config):
        """
        Parameter setting.
        """
        spatialRulesSetter = ParameterSpatialRulesSet(
            self.RULES_SET,
            description=self.tr('Spatial Rules Set')
        )
        spatialRulesSetter.setMetadata({
            'widget_wrapper' : 'DsgTools.gui.ProcessingUI.enforceSpatialRuleWrapper.EnforceSpatialRuleWrapper'
        })
        self.addParameter(spatialRulesSetter)

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.POINT_FLAGS,
                self.tr('{0} flags (points)').format(self.displayName())
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.LINE_FLAGS,
                self.tr('{0} flags (lines)').format(self.displayName())
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.POLYGON_FLAGS,
                self.tr('{0} flags (polygons)').format(self.displayName())
            )
        )

    def parameterAsSpatialRulesSet(self, parameters, name, context):
        return parameters[name]

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'enforcespatialrules'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Enforce Spatial Rules')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('Quality Assurance Tools (Identification Processes)')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'DSGTools: Quality Assurance Tools (Identification Processes)'

    def tr(self, string):
        return QCoreApplication.translate('EnforceSpatialRulesAlgorithm', string)

    def createInstance(self):
        return EnforceSpatialRulesAlgorithm()

    def setupLayer(self, layername, expression, context=None, feedback=None):
        """
        Reads layer from canvas and applies the filtering expression.
        :param layername: (str) layer to be setup.
        :param expression: (str) filtering expression to be applied to the
                           target layer.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) vector layer ready to be set up.
        """
        if layername not in self.__layers:
            # by default, it is assumed layer names are unique on canvas
            vl = QgsProject.instance().mapLayersByName(layername)
            if not vl:
                raise QgsProcessingException(
                    self.tr("Layer {l} was not found!").format(l=layername)
                )
            # layer caching happens here
            self.__layers[layername] = vl[0]
        vl = self.__layers[layername]
        if expression:
            vl = AlgRunner().runFilterExpression(
                vl,
                expression,
                context or QgsProcessingContext(),
                outputLyr='memory:',
                feedback=feedback
            )
        return vl

    def disjoint(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A is disjoint to layer B.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A are not disjoint to features
                 of layer B.
        """
        # disjoint is kept to keep parallel with other methods. it only counts
        # completely dissociated features (obviously)
        return QgsVectorLayer()

    def contains(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A contains layer B, respecting its cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A not contanining features of
                 layer B.
        """
        return QgsVectorLayer()

    def notContains(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A does not contain layer B, respecting its cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A containing features of layer B.
        """
        return QgsVectorLayer()

    def intersects(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A intersects layer B, respecting its cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A not intersecting features of
                 layer B.
        """
        return QgsVectorLayer()

    def notIntersects(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A does not intersect layer B, respecting its
        cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A intersecting features of
                 layer B.
        """
        return QgsVectorLayer()

    def isContained(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A is contained by layer B, respecting its cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A not contained by features of
                 layer B.
        """
        return QgsVectorLayer()

    def isNotContained(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A is not contained by layer B, respecting its
        cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A being contained by layer B.
        """
        return QgsVectorLayer()

    def touches(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A touches layer B, respecting its cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A not touching layer B.
        """
        return QgsVectorLayer()

    def notTouches(self, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Checks if layer A does not touch layer B, respecting its cardinality.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences features of layer A touching layer B.
        """
        return QgsVectorLayer()

    def verifyTopologicalRelation(self, predicate, layerA, layerB, cardinality, context=None, feedback=None):
        """
        Verifies a given topological relation between two layers, returning a
        vector layer pointing where these predicates have happened.
        :param predicate: (int) topological relation to be tested.
        :param layerA: (QgsVectorLayer) reference layer to be tested.
        :param layerB: (QgsVectorLayer) the other layer to be tested.
        :param cardinality: (str) a string representing the upper and lower
                            bounds of occurrences an event can happen between
                            the target layers.
        :param context: (QgsProcessingContext) environment context in which 
                        layer is retrieved and setup.
        :param feedback: (QgsProcessingFeedback) QGIS progress tracking
                         component.
        :return: (QgsVectorLayer) layer containing features representing the
                 occurrences of the given test (their flags).
        """
        handler = SpatialRelationsHandler()        
        if predicate not in range(len(self.availablePredicates())):
            raise QgsProcessingException(
                self.tr("Unsupported relation {n}").format(n=predicate)
            )
        return {
            self.EQUALS : handler.isEqual,
            self.NOTEQUALS : handler.disjoint,
            self.DISJOINT : handler.disjoint,
            self.NOTDISJOINT : handler.disjoint,
            self.INTERSECTS : handler.disjoint,
            self.NOTINTERSECTS : handler.disjoint,
            self.TOUCHES : handler.disjoint,
            self.NOTTOUCHES : handler.disjoint,
            self.CROSSES : handler.disjoint,
            self.NOTCROSSES : handler.disjoint,
            self.WITHIN : handler.disjoint,
            self.NOTWITHIN : handler.disjoint,
            self.OVERLAPS : handler.contains,
            self.NOTOVERLAPS : handler.notContains,
            self.CONTAINS : handler.intersects,
            self.NOTCONTAINS : handler.notIntersects,
            self.NOTCOVERS : handler.isContained,
            self.ISNOTCONTAINED : handler.isNotContained,
            self.COVEREDBY : handler.touches,
            self.NOTCOVEREDBY : handler.notTouches
        }[predicate](layerA, layerB, cardinality, context, feedback)

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        rules = self.parameterAsSpatialRulesSet(
            parameters, self.RULES_SET, context
        )
        multiStepFeedback = QgsProcessingMultiStepFeedback(len(rules), feedback)
        multiStepFeedback.setCurrentStep(0)
        for current, rule in enumerate(rules):
            if feedback.isCanceled() or multiStepFeedback.isCanceled():
                return {}
            name = rule["name"]
            multiStepFeedbackInner = QgsProcessingMultiStepFeedback(
                2, multiStepFeedback
            )
            out = self.verifyTopologicalRelation(
                rule["predicate"],
                self.setupLayer(
                    rule["layer_a"], rule["filter_a"], context, None
                ),
                self.setupLayer(
                    rule["layer_b"], rule["filter_b"], context, None
                ),
                rule["cardinality"],
                context,
                multiStepFeedbackInner
            )
            multiStepFeedbackInner.setCurrentStep(1)
            if out.featureCount() > 0:
                # raise flags
                pass
            multiStepFeedbackInner.setCurrentStep(2)
            multiStepFeedback.setCurrentStep(current + 1)
        # dont forget to clear cached layers
        self.__layers = dict()
        return {}

class ParameterSpatialRulesSetType(QgsProcessingParameterType):

    def __init__(self):
        super().__init__()

    def create(self, name):
        return ParameterSpatialRulesSet(name)

    def metadata(self):
        return {'widget_wrapper': 'DsgTools.gui.ProcessingUI.enforceSpatialRuleWrapper.EnforceSpatialRuleWrapper'}

    def name(self):
        return QCoreApplication.translate('Processing', 'Spatial Rules Set')

    def id(self):
        return 'spatial_rules_set_type'

    def description(self):
        return QCoreApplication.translate('Processing', 'Set of spatial rules. Used on Spatial Rules Checker.')

class ParameterSpatialRulesSet(QgsProcessingParameterDefinition):

    def __init__(self, name, description=''):
        super().__init__(name, description)

    def clone(self):
        copy = ParameterSpatialRulesSet(self.name(), self.description())
        return copy

    def type(self):
        return self.typeName()

    @staticmethod
    def typeName():
        return 'spatial_rules_set'

    def checkValueIsAcceptable(self, value, context=None):
        return True

    def valueAsPythonString(self, value, context):
        return str(value)

    def asScriptCode(self):
        raise NotImplementedError()

    @classmethod
    def fromScriptCode(cls, name, description, isOptional, definition):
        raise NotImplementedError()
