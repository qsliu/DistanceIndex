# -*- coding: utf-8 -*-
"""
Program Author: Qingsong Liu
Email: qliu20@kent.edu

This program was implemented based on Clark, Philip J, and Francis C Evans. 1954.

Reference:
Clark, Philip J, and Francis C Evans. 1954.
“Distance to Nearest Neighbor as a Measure of Spatial Relationships in Populations.”
Ecology 35 (4): 445–53. https://doi.org/10.2307/1931034.
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtCore import QDate,QDateTime,QTime
from qgis.core import (QgsProcessing,
                       QgsSpatialIndex,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFileDestination,
                       QgsProcessingParameterField,
                       QgsProcessingOutputNumber,
                       QgsFeatureRequest,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterCrs,
                       QgsProcessingParameterDistance,
                       QgsCoordinateReferenceSystem,
                       QgsProject,
                       QgsCoordinateTransform,
                       QgsDistanceArea)
from qgis import processing
from rtree import index
import math,tempfile,os
import codecs,re
class NNIAlgorithm(QgsProcessingAlgorithm):

    INPUT_LAYER = 'INPUT_LAYER'
    OUTPUT_HTML_FILE = 'OUTPUT_HTML_FILE'
    DISTANCE_UNIT = "DISTANCE_UNIT"
    NUM_TIME_UNITS = "NUM_TIME_UNITS"
    NNI_OUT = "NNI_RES"
    CRS = "CRS"
    BASE_LAYER = "BASE_LAYER"


    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return NNIAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Nearest Neighborhood Index'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Nearest Neighborhood Index')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('Distance Indicators')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'DistanceIndicators'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr("Nearest Neighborhood Index")

    def initAlgorithm(self, config=None):

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_LAYER,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorPoint]
            )
        )

        self.addParameter(QgsProcessingParameterDistance(self.DISTANCE_UNIT,self.tr("Distance defined as 1 unit in the calcualtion"), 1,
                                                         self.CRS, False, 0))
        # "ProjectCrs" make sure
        self.addParameter(QgsProcessingParameterCrs(self.CRS, self.tr('Projection used in the calculation'), 'ProjectCrs'))

        self.addParameter(QgsProcessingParameterFeatureSource(self.BASE_LAYER, self.tr('Base layer'),[QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterFileDestination(self.OUTPUT_HTML_FILE, self.tr('Statistics'),
                                                                self.tr('HTML files (*.html)'), None, True))

        self.addOutput(QgsProcessingOutputNumber(self.NNI_OUT, self.tr('NNI')))

    def distance(self,x1,y1,z1, x2,y2,z2):
        return math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT_LAYER, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_LAYER))

        # if you change the unit type, the program will automatically change for you
        numDistancePerUnit = self.parameterAsDouble(parameters, self.DISTANCE_UNIT, context)

        target_crs = self.parameterAsCrs(parameters, self.CRS, context)
        baselayer = self.parameterAsSource(parameters, self.BASE_LAYER, context)
        # Minimum Bounding Rectangular
        output_file = self.parameterAsFileOutput(parameters, self.OUTPUT_HTML_FILE, context)

        # Send some information to the user
        feedback.pushInfo('Start reading layer information...')

        input_crs = source.sourceCrs()
        if input_crs != target_crs:
            transtotarget = QgsCoordinateTransform(input_crs,target_crs,QgsProject.instance())

        index = QgsSpatialIndex(source.getFeatures(QgsFeatureRequest().setSubsetOfAttributes([]).
                                                   setDestinationCrs(target_crs, context.transformContext())))
        distance = QgsDistanceArea()
        distance.setSourceCrs(target_crs, context.transformContext())
        distance.setEllipsoid(context.ellipsoid())

        # calculate the obseved NN
        request = QgsFeatureRequest().setSubsetOfAttributes([])
        features = source.getFeatures(request)
        count = source.featureCount()
        r_obs = 0.0
        total = 100.0 / count if count else 0
        for current, f in enumerate(features):
            if feedback.isCanceled():
                break
            src = f.geometry().boundingBox().center()

            neighbors = index.nearestNeighbor(src, 2)
            neighbors.remove(f.id())
            ft = next(source.getFeatures(
                                QgsFeatureRequest().setFilterFid(neighbors[0]).
                                             setSubsetOfAttributes([], source.fields()).
                                             setDestinationCrs(target_crs, context.transformContext())
                                             )
                      )
            closest = ft.geometry().boundingBox().center()
            Dist = src.distance( closest ) / numDistancePerUnit
            r_obs += Dist
        r_obs = r_obs/count

        # volume calculation
        request = QgsFeatureRequest().setSubsetOfAttributes([])
        features = baselayer.getFeatures(request)

        baselayer_crs = baselayer.sourceCrs()
        if baselayer_crs != target_crs:
            transtotarget = QgsCoordinateTransform(baselayer_crs,target_crs,QgsProject.instance())
        da = QgsDistanceArea()
        da.setSourceCrs(target_crs, context.transformContext())
        sum_area = 0.0
        for current, f in enumerate(features):
            g = f.geometry()
            if baselayer_crs != target_crs:
                g.transform(transtotarget)
            g_area = da.measureArea(g)
            sum_area += g_area / numDistancePerUnit / numDistancePerUnit
        # calculate theory NN
        # rho = count/sum_area
        # r_e = 0.5/rho**(1/2)
        r_e = 0.5 * sum_area ** (1 / 2)/ count ** (1 / 2)
        # calculate the final nni
        res = r_obs/r_e
        sd_r_e = 0.26136 * sum_area ** (1 / 2)/(count * count)**(1/2)
        z_score = (r_obs - r_e)/sd_r_e

        data = []
        data.append(self.tr('Nearest Neighborhood Index (NNI):  {}').format(res))
        data.append(self.tr('r_obs:  {}').format(r_obs))
        data.append(self.tr('r_e:  {}').format(r_e))
        data.append(self.tr('Standard deviation of r_e:  {}').format(sd_r_e))
        data.append(self.tr('Nearest Neighborhood Index Z-score:  {}').format(z_score))

        results = {}
        results[self.NNI_OUT] = res

        if output_file:
            self.createHTML(output_file, data)
            results[self.OUTPUT_HTML_FILE] = output_file

        return results

    def createHTML(self, outputFile, algData):
        with codecs.open(outputFile, 'w', encoding='utf-8') as f:
            f.write('<html><head>\n')
            f.write('<meta http-equiv="Content-Type" content="text/html; \
                    charset=utf-8" /></head><body>\n')
            for s in algData:
                f.write('<p>' + str(s) + '</p>\n')
            f.write('</body></html>\n')
