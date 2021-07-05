# -*- coding: utf-8 -*-
"""
Program Author: Qingsong Liu
Email: qliu20@kent.edu

This program was implemented based on Jay Lee's article.

Reference:
2021 Lee, J., S. Li, S. Wang, J. Wang, and J. Li. Spatio-Temporal Nearest neighbor Index for Measuring Space-Time Clustering among Geographic Events. Papers in Applied Geography. https://doi.org/10.1080/23754931.2020.1810112.
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtCore import QDate,QDateTime,QTime
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
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
class STNNIAlgorithm(QgsProcessingAlgorithm):

    INPUT_LAYER = 'INPUT_LAYER'
    DATE_FIELD_NAME = 'DATE_FIELD_NAME'
    OUTPUT_HTML_FILE = 'OUTPUT_HTML_FILE'
    DISTANCE_UNIT = "DISTANCE_UNIT"
    NUM_TIME_UNITS = "NUM_TIME_UNITS"
    TIME_UNITS = "TIME_UNITS"
    STNNI_OUT = "STNNI_RES"
    CRS = "CRS"
    BASE_LAYER = "BASE_LAYER"


    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return STNNIAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Space-time Nearest Neighborhood Index'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Space-time Nearest Neighborhood Index')

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
        return self.tr("Space-time Nearest Neighborhood Index")

    def initAlgorithm(self, config=None):
        self.time_TIME_UNITS = [self.tr("Year"),
                           self.tr("Month"),
                           self.tr("Day"),
                           self.tr("Hour"),
                           self.tr("Minute"),
                           self.tr("Second")]

        # We add the input vector features source. It can have any kind of
        # geometry.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_LAYER,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorPoint]
            )
        )


        self.addParameter(QgsProcessingParameterField(self.DATE_FIELD_NAME,
                                                      self.tr('Datetime info field (must in yyyy-mm-dd hh:mm:ss format)'),
                                                      None, self.INPUT_LAYER, QgsProcessingParameterField.String))

        self.addParameter(QgsProcessingParameterDistance(self.DISTANCE_UNIT,self.tr("Distance defined as 1 unit in the calcualtion"), 1,
                                                         self.CRS, False, 0))
        # "ProjectCrs" make sure
        self.addParameter(QgsProcessingParameterCrs(self.CRS, self.tr('Projection used in the calculation'), 'ProjectCrs'))

        self.addParameter(QgsProcessingParameterNumber(self.NUM_TIME_UNITS,
                                                       self.tr('Number of time unit'), type=QgsProcessingParameterNumber.Double, minValue=0, defaultValue=1))

        self.addParameter(QgsProcessingParameterEnum(self.TIME_UNITS,
                                                     self.tr('Calculation Time Unit'), options=self.time_TIME_UNITS, defaultValue=2))

        self.addParameter(QgsProcessingParameterFeatureSource(self.BASE_LAYER, self.tr('Base layer'),[QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterFileDestination(self.OUTPUT_HTML_FILE, self.tr('Statistics'),
                                                                self.tr('HTML files (*.html)'), None, True))

        self.addOutput(QgsProcessingOutputNumber(self.STNNI_OUT, self.tr('STNNI')))
        self.checkedtime = False
        self.timepattern = None

    def string2datetime(self,value):
        #only check the first record
        if not self.checkedtime:
            pattern1 = re.compile("\s*\d{4}-\d{2}-\d{2}\s*\d{2}:\d{2}:\d{2}\s*")
            pattern2 = re.compile("\s*\d{4}-\d{2}-\d{2}\s*")
            pattern3 = re.compile("\s*\d{2}:\d{2}:\d{2}\s*")
            if pattern1.match(value):
                self.timepattern = "datetime"
            elif pattern2.match(value):
                self.timepattern = "date"
            elif pattern3.match(value):
                self.timepattern = "time"
            else:
                self.timepattern = None
            self.checkedtime = True
        if self.timepattern == "datetime":
            return QDateTime.fromString(value,"yyyy-MM-dd hh:mm:ss")
        elif self.timepattern == "date":
            return QDateTime(QDate.fromString(value,"yyyy-MM-dd"))
        elif self.timepattern == "time":
            now = QDateTime.currentDateTime()
            return QDateTime(now.date(), QTime.fromString(value, "hh:mm:ss"))
        else:
            return None

    def convertTime(self, cur, basetime, unit, numperunit):
        res = None
        if unit == "Year":
            diff = cur.date().year() - basetime.date().year()
            res = diff / numperunit
        elif unit == "Month":
            diff = (cur.date().year() - basetime.date().year()-1)*12
            diff += cur.date().month()
            diff += 12 - basetime.date().month()
            res = diff / numperunit
        elif unit == "Day":
            diff = basetime.date().daysTo(cur.date())
            res = diff / numperunit
        elif unit == "Hour":
            diff = (cur.toSecsSinceEpoch() - basetime.toSecsSinceEpoch())/60.0/60.0
            res = diff / numperunit
        elif unit == "Minute":
            diff = (cur.toSecsSinceEpoch() - basetime.toSecsSinceEpoch())/60.0
            res = diff / numperunit
        elif unit == "Second":
            diff = cur.toSecsSinceEpoch() - basetime.toSecsSinceEpoch()
            res = diff / numperunit
        else:
            res = None
        return res

    def distance(self,x1,y1,z1, x2,y2,z2):
        return math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT_LAYER, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_LAYER))

        TimeUnitType = self.parameterAsEnum(parameters, self.TIME_UNITS, context)
        numTimePerUnit = self.parameterAsDouble(parameters, self.NUM_TIME_UNITS, context)

        # if you change the unit type, the program will automatically change for you
        numDistancePerUnit = self.parameterAsDouble(parameters, self.DISTANCE_UNIT, context)

        field_name = self.parameterAsString(parameters, self.DATE_FIELD_NAME, context)
        fieldidx = source.fields().lookupField(field_name)

        target_crs = self.parameterAsCrs(parameters, self.CRS, context)
        baselayer = self.parameterAsSource(parameters, self.BASE_LAYER, context)
        # Minimum Bounding Rectangular
        output_file = self.parameterAsFileOutput(parameters, self.OUTPUT_HTML_FILE, context)

        request = QgsFeatureRequest().setSubsetOfAttributes([field_name],source.fields())
        features = source.getFeatures(request)
        count = source.featureCount()

        # Send some information to the user
        feedback.pushInfo('Start reading layer information...')

        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / count if count else 0
        xlist = []
        ylist = []
        zlist = []
        zmin = float('inf')
        zminidx = None
        fids = []
        fids2index = {}

        input_crs = source.sourceCrs()
        if input_crs != target_crs:
            transtotarget = QgsCoordinateTransform(input_crs,target_crs,QgsProject.instance())

        for current, feature in enumerate(features):
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break
            fid = feature.id()
            fids.append(fid)
            fids2index[fid] = current

            cur_geo = feature.geometry()
            if input_crs != target_crs:
                cur_geo.transform(transtotarget)
            xlist.append(cur_geo.asPoint().x())
            ylist.append(cur_geo.asPoint().y())

            zattr = feature.attributes()[fieldidx]

            c_zattr = self.string2datetime(zattr)
            if c_zattr.date().year()==0:
                raise QgsProcessingException(
                    self.tr('Please check time field {} with value of {}').format(field_name,zattr))
            if c_zattr==None:
                raise QgsProcessingException(
                    self.tr('Please check time field {} with value of {}').format(field_name,zattr))
            zlist.append(c_zattr)
            zSeconds = c_zattr.toSecsSinceEpoch()
            if zmin > zSeconds:
                zmin = zSeconds
                zminidx = current
            # Update the progress bar
            feedback.setProgress(int(current * total))

        feedback.pushInfo('Convert space and time coordinate to defined units...')
        basedatetime = zlist[zminidx]
        xlist = [ele / numDistancePerUnit for ele in xlist]
        ylist = [ele / numDistancePerUnit for ele in ylist]
        zlist = [self.convertTime(ele,basedatetime,self.time_TIME_UNITS[TimeUnitType],numTimePerUnit) for ele in zlist]

        feedback.pushInfo('Construct spatial index...')
        # construct spatial index
        p = index.Property()
        p.dat_extension = 'data'
        p.idx_extension = 'index'
        fd = tempfile.gettempdir()
        p.dimension = 3
        dx3d = index.Index(os.path.join(fd,'3d_index'), properties=p)
        for current, f in enumerate(fids):
            dx3d.insert(f, (xlist[current],ylist[current],zlist[current],xlist[current],ylist[current],zlist[current]))

        feedback.pushInfo('calculate the obseved STNNI...')
        # calculate the obseved NN
        r_obs = 0.0
        for current,f in enumerate(fids):
            nn = list(dx3d.nearest((xlist[current],ylist[current],zlist[current],xlist[current],ylist[current],zlist[current]), 2))
            nn.remove(current)

            r_obs += self.distance(xlist[current],ylist[current],zlist[current],
                                   xlist[fids2index[nn[0]]],ylist[fids2index[nn[0]]],zlist[fids2index[nn[0]]])
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
        z_max = max(zlist)
        volume = sum_area * z_max
        if volume == 0: volume = 1
        # rho = count/volume
        # r_e = 0.55396/rho**(1/3)
        r_e = 0.55396 * volume** (1 / 3) / count ** (1 / 3)
        # calculate the final space-time nni
        res = r_obs/r_e
        sd_r_e = 0.04054**(1/2)*volume** (1 / 3)/count**(1/3)
        z_score = (r_obs - r_e)/sd_r_e

        feedback.pushInfo('Printing output STNNI...\n\n')
        data = []
        data.append(self.tr('Time information field: {}').format(field_name))
        data.append(self.tr('Space-time Nearest Neighborhood Index:  {}').format(res))
        data.append(self.tr('r_obs:  {}').format(r_obs))
        data.append(self.tr('r_e:  {}').format(r_e))
        data.append(self.tr('Standard deviation of r_e:  {}').format(sd_r_e))
        data.append(self.tr('Space-time Nearest Neighborhood Index Z-score:  {}').format(z_score))

        results = {}
        results[self.STNNI_OUT] = res

        if output_file:
            self.createHTML(output_file, data)
            results[self.OUTPUT_HTML_FILE] = output_file
        del dx3d
        if os.path.exists(os.path.join(fd,'3d_index.index')):
            os.remove(os.path.join(fd,'3d_index.index'))
        if os.path.exists(os.path.join(fd, '3d_index.data')):
            os.remove(os.path.join(fd, '3d_index.data'))
        return results

    def createHTML(self, outputFile, algData):
        with codecs.open(outputFile, 'w', encoding='utf-8') as f:
            f.write('<html><head>\n')
            f.write('<meta http-equiv="Content-Type" content="text/html; \
                    charset=utf-8" /></head><body>\n')
            for s in algData:
                f.write('<p>' + str(s) + '</p>\n')
            f.write('</body></html>\n')
