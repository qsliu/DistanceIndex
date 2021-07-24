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
                       QgsUnitTypes,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFileDestination,
                       QgsProcessingParameterField,
                       QgsProcessingOutputNumber,
                       QgsFeatureRequest,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterString,
                       QgsProcessingParameterBoolean,
                       QgsProject,
                       QgsCoordinateTransform,
                       QgsDistanceArea)
from rtree import index
import tempfile,os
import codecs,re
class STNNIAlgorithm(QgsProcessingAlgorithm):

    INPUT_LAYER = 'INPUT_LAYER'
    DATE_FIELD_NAME = 'DATE_FIELD_NAME'
    DATE_FIELD_FORMAT = "DATE_FIELD_FORMAT"
    OUTPUT_HTML_FILE = 'OUTPUT_HTML_FILE'
    DISTANCE_UNIT = "DISTANCE_UNIT"
    NUM_TIME_UNITS = "NUM_TIME_UNITS"
    TIME_UNITS = "TIME_UNITS"
    STNNI_OUT = "STNNI_RES"
    CRS = "CRS"
    BASE_LAYER = "BASE_LAYER"
    BOOL_USE_PRO_CRS = "BOOL_USE_PRO_CRS"


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

        self.addParameter(QgsProcessingParameterString(self.DATE_FIELD_FORMAT,
                                                       self.tr('Datetime format: e.g., yyyy-mm-dd hh:mm:ss'),optional=True))

        self.addParameter(QgsProcessingParameterNumber(self.DISTANCE_UNIT,
                                                       self.tr('Define distance unit in meters (How many meters per unit?)'),
                                                       type=QgsProcessingParameterNumber.Double,
                                                       minValue=1.0,
                                                       defaultValue=1.0))

        self.addParameter(QgsProcessingParameterNumber(self.NUM_TIME_UNITS,
                                                       self.tr('Number of time unit'), type=QgsProcessingParameterNumber.Double, minValue=0, defaultValue=1))

        self.addParameter(QgsProcessingParameterEnum(self.TIME_UNITS,
                                                     self.tr('Calculation Time Unit'), options=self.time_TIME_UNITS, defaultValue=2))

        self.addParameter(QgsProcessingParameterFeatureSource(self.BASE_LAYER, self.tr('Base layer'),[QgsProcessing.TypeVectorPolygon]))

        self.addParameter(QgsProcessingParameterBoolean(self.BOOL_USE_PRO_CRS,
                                                        self.tr('Use on-the-fly map Projection in calculation'),
                                                        defaultValue=False))
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

    def processAlgorithm(self, parameters, context, feedback):

        source = self.parameterAsSource(parameters, self.INPUT_LAYER, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_LAYER))

        TimeUnitType = self.parameterAsEnum(parameters, self.TIME_UNITS, context)
        numTimePerUnit = self.parameterAsDouble(parameters, self.NUM_TIME_UNITS, context)

        # if you change the unit type, the program will automatically change for you
        numDistancePerUnit = self.parameterAsDouble(parameters, self.DISTANCE_UNIT, context)

        field_name = self.parameterAsString(parameters, self.DATE_FIELD_NAME, context)
        self.field_format = self.parameterAsString(parameters,self.DATE_FIELD_FORMAT,context)

        fieldidx = source.fields().lookupField(field_name)

        baselayer = self.parameterAsSource(parameters, self.BASE_LAYER, context)
        baselayer_crs = baselayer.sourceCrs()

        bool_prj_crs = self.parameterAsBoolean(parameters,self.BOOL_USE_PRO_CRS,context)

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
        geoList = {}
        xlist = []
        ylist = []
        zlist = []
        zmin = float('inf')
        zminidx = None
        fids = []
        fids2index = {}

        input_crs = source.sourceCrs()
        prj_crs = QgsProject.instance().crs()

        if bool_prj_crs:
            if prj_crs.isGeographic():
                raise QgsProcessingException('Only projection maps are supported now, '
                                         + 'please set on-the-fly map Projection for the project!')
        else:
            # if not use the on-the-fly projection (project projection), and the input point layer have different CRS compare to base layer, then there is an error.
            if input_crs != baselayer_crs:
                raise QgsProcessingException(
                self.tr('Point layer and polygon layer coordinate systems do not coincide \n and the STNNI cannot be calculated.\n'
                        +'You can solve it by check Use on-the-fly map Projection in calculation!'))
            else:
                if input_crs.isGeographic():
                    raise QgsProcessingException('Only projection maps are supported now, not geographic coordinate systems.'
                                             + "you could set on-the-fly map Projection for the project or set projection for each layer!")
        transtotarget = None
        if bool_prj_crs ==True and input_crs != prj_crs:
            transtotarget = QgsCoordinateTransform(input_crs,prj_crs,QgsProject.instance())

        for current, feature in enumerate(features):
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break
            fid = feature.id()
            fids.append(fid)
            fids2index[fid] = current

            cur_geo = feature.geometry()
            if transtotarget!=None:
                cur_geo.transform(transtotarget)
            geoList[fid] = cur_geo
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
        zlist = [self.convertTime(ele,basedatetime,self.time_TIME_UNITS[TimeUnitType],numTimePerUnit) for ele in zlist]

        feedback.pushInfo('Construct spatial index...')
        # construct spatial index
        p = index.Property()
        p.dat_extension = 'data'
        p.idx_extension = 'index'
        fd = tempfile.gettempdir()
        p.dimension = 3
        dx3d = index.Index(os.path.join(fd,'3d_index'), properties=p,overwrite=True)
        for current, f in enumerate(fids):
            dx3d.insert(f, (xlist[fids2index[f]],ylist[fids2index[f]],zlist[fids2index[f]],
                            xlist[fids2index[f]],ylist[fids2index[f]],zlist[fids2index[f]]))

        feedback.pushInfo('calculate the obseved STNNI...')
        # calculate the obseved NN
        r_obs = 0.0
        distance = QgsDistanceArea()
        cur_prj = input_crs
        if bool_prj_crs and input_crs!=prj_crs:
            cur_prj = prj_crs

        distance.setSourceCrs(cur_prj, context.transformContext())
        distance.setEllipsoid(cur_prj.ellipsoidAcronym ())

        o_units = distance.lengthUnits()
        o2meter_factor = QgsUnitTypes.fromUnitToUnitFactor(o_units, QgsUnitTypes.DistanceMeters)

        for current,f in enumerate(fids):
            nn = list(dx3d.nearest((xlist[fids2index[f]],ylist[fids2index[f]],zlist[fids2index[f]],
                                    xlist[fids2index[f]],ylist[fids2index[f]],zlist[fids2index[f]]), 2))
            nn.remove(f)

            if cur_prj.isGeographic():
                xydis = distance.measureLine(geoList[f].asPoint(), geoList[nn[0]].asPoint())*o2meter_factor/numDistancePerUnit
            else:
                xydis = geoList[f].asPoint().distance(geoList[nn[0]].asPoint()) * o2meter_factor / numDistancePerUnit

            r_obs += (xydis**2 + (zlist[current]-zlist[fids2index[nn[0]]])**2)**(1/2)

        r_obs = r_obs/count

        # volume calculation
        request = QgsFeatureRequest().setSubsetOfAttributes([])
        features = baselayer.getFeatures(request)


        if baselayer_crs != prj_crs:
            trans2 = QgsCoordinateTransform(baselayer_crs,prj_crs,QgsProject.instance())
        da = QgsDistanceArea()
        if bool_prj_crs:
            da.setSourceCrs(prj_crs, context.transformContext())
            da.setEllipsoid(prj_crs.ellipsoidAcronym())
        else:
            da.setSourceCrs(baselayer_crs, context.transformContext())
            da.setEllipsoid(baselayer_crs.ellipsoidAcronym())

        sum_area = 0.0
        da.areaUnits()
        o_units = da.areaUnits()
        o2squaremeter_factor = QgsUnitTypes.fromUnitToUnitFactor(o_units, QgsUnitTypes.AreaSquareMeters)
        for current, f in enumerate(features):
            g = f.geometry()
            if bool_prj_crs==True and baselayer_crs != prj_crs:
                g.transform(trans2)
            g_area = da.measureArea(g) * o2squaremeter_factor
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
        data.append(self.tr('Python command \n{}').format(self.asPythonCommand(parameters,context)))
        data.append(self.tr('Point layer: {}').format(source))
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
