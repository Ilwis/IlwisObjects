import unittest
import ilwis
import config
import numpy as np
import math
from os.path import abspath


def testExceptionCondition6(p, func, parm1, parm2, parm3, parm4, parm5, parm6, message):
      try:
           rc3 = func(parm1, parm2, parm3, parm4, parm5, parm6)
           p.isTrue(False, message)
      except ilwis.IlwisException as ex:
           p.isTrue(True, message)

def testExceptionCondition5(p, func, parm1, parm2, parm3, parm4, parm5, message):
      try:
           rc3 = func(parm1, parm2, parm3, parm4, parm5)
           p.isTrue(False, message)
      except ilwis.IlwisException as ex:
           p.isTrue(True, message)

def testExceptionCondition4(p, func, parm1, parm2, parm3, parm4, message):
      try:
           rc3 = func(parm1, parm2, parm3, parm4)
           p.isTrue(False, message)
      except ilwis.IlwisException as ex:
           p.isTrue(True, message)

def testExceptionCondition3(p, func, parm1, parm2, parm3, message):
      try:
           rc3 = func(parm1, parm2, parm3)
           p.isTrue(False, message)
      except ilwis.IlwisException as ex:
           p.isTrue(True, message)

def testExceptionCondition2(p, func, parm1, parm2, message):
      try:
           rc3 = func(parm1, parm2)
           p.isTrue(False, message)
      except ilwis.IlwisException as ex:
           p.isTrue(True, message)

def testExceptionCondition1(p, func, parm1, message):
      try:
           rc3 = func(parm1)
           p.isTrue(False, message)
      except ilwis.IlwisException as ex:
           p.isTrue(True, message)           

class BaseTest(unittest.TestCase):
    def prepare(self, testdir):
        try:
            ilwis.disconnectIssueLogger()
            workingdir = abspath('.') + '/' + testdir
            ilwis.Engine.setWorkingCatalog(workingdir)
            ilwis.connectIssueLogger()
        except ilwis.IlwisException:
            ilwis.connectIssueLogger()
            self.skipTest("could not set working directory!")

    def isEqual(self, str1, str2, msg):
        config.testCount += 1 
        result = 'SUCCESS'
        if (str1 != str2):
           config.fails.append(str(config.testCount) + ":" + self.decoration + " : "+ msg)
           result = 'FAIL'

        print(f'{config.testCount:5} {msg:65}  {result}')

    def isAlmostEqualNum(self, num1, num2, delta, msg) :
        config.testCount += 1 
        result = 'SUCCESS'
        if (abs(num1 - num2) > delta):
           config.fails.append(str(config.testCount) + ":" + self.decoration + " : "+ msg)
           result = 'FAIL'

        print(f'{config.testCount:5} {msg:65}  {result}')

    def isAlmostEqualEnvelope(self, env1, env2, delta, msg):
        minx1 = env1.minCorner().x
        miny1 = env1.minCorner().y
        maxx1 = env1.maxCorner().x
        maxy1 = env1.maxCorner().y 
        minx2 = env2.minCorner().x
        miny2 = env2.minCorner().y
        maxx2 = env2.maxCorner().x
        maxy2 = env2.maxCorner().y 

        config.testCount += 1 
        result = 'SUCCESS'
        if (abs(minx1 - minx2) > delta or abs(maxx1 - maxx2) > delta or abs(miny1 - miny2) > delta or abs(maxy1 - maxy2) > delta):
           config.fails.append(str(config.testCount) + ":" + self.decoration + " : "+ msg)
           result = 'FAIL'
        
        print(f'{config.testCount:5} {msg:65}  {result}')



    def isTrue(self, b, msg):
        config.testCount += 1 
        result = 'FAIL'
        if (b):
           result = 'SUCCESS'
        else:
            config.fails.append(str(config.testCount) + ":" + self.decoration + ":"+ msg)         
        print(f'{config.testCount:5} {msg:65}  {result}')

    def isFalse(self, b, msg):
        config.testCount += 1 
        result = 'SUCCESS'
        if (b):
            config.fails.append(str(config.testCount) + ":" +self.decoration + ":"+ msg)
            result = 'FAIL'
        print(f'{config.testCount:5} {msg:65}  {result}') 
                     

    def decorateFunction(self, mod, fn) :
        self.decoration = mod + " ==> " + fn 
        print("\n" + self.decoration + "\n")

    def createEmptySmallNumericRaster(self):
        grf = ilwis.GeoReference("epsg:4326", ilwis.Envelope("0 25 30 60") , ilwis.Size(15,12))
        dfNum = ilwis.DataDefinition(ilwis.NumericDomain("code=value"), ilwis.NumericRange(0.0, 1000000.0, 1.0))
        rc = ilwis.RasterCoverage()
        rc.setGeoReference(grf)
        rc.setDataDef(dfNum)

        return rc

    def createThematicDomain(self):
        tr = ilwis.ThematicRange()
        tr.add("grass", "1", "mostely green")
        tr.add("stone", "2", "greyish")
        tr.add("houses", "3", "mixed colors")
        tr.add("water", "4", "blueish")

        td = ilwis.ItemDomain(tr) 

        return td       

    def createEmptySmallThematicRaster(self):
        
        td = self.createThematicDomain()
      
     
        self.dfThematic = ilwis.DataDefinition(td)
        grf = ilwis.GeoReference("epsg:4326", ilwis.Envelope("0 25 30 60") , ilwis.Size(15,12))
        tbl = self.createKeyedTestTable()
        rc = ilwis.RasterCoverage()
        rc.setDataDef(td)
        rc.setGeoReference(grf)
        rc.setAttributes(tbl, "items")
   
        return rc        

    def createSmallThematicRaster1Layer(self):
        rc = self.createEmptySmallThematicRaster()
        rc.setSize(ilwis.Size(15,12))
        baseSize = 15 * 12
        array1 = np.empty(baseSize, dtype = int)

        for i in range(len(array1)):
            array1[i] = i % 3

        rc.array2raster(array1)
     
        return rc    

    def createSmallNumericRaster1Layer(self,offset=0):
        rc = self.createEmptySmallNumericRaster()
        rc.setSize(ilwis.Size(15,12))
        baseSize = 15 * 12
        array1 = np.empty(baseSize, dtype = np.int64)

        for i in range(len(array1)):
            array1[i] = i * 10 + offset

        array1[9] = ilwis.Const.iUNDEF # pos 8,0 is undefined
        rc.array2raster(array1)
     
        return rc

    def createSmallNumericRaster3Layers(self, alternate=0):
        rc = self.createEmptySmallNumericRaster()
        rc.setSize(ilwis.Size(15,12,3))
        baseSize = 15 * 12

        array1 = np.empty(baseSize, dtype = np.int64)
        array2 = np.empty(baseSize, dtype = np.int64)
        array3 = np.empty(baseSize, dtype = np.int64)

        

        if ( alternate == 0):
            for i in range(len(array1)):
                array1[i] = i * 10
                array2[i] = (i + baseSize) * 10
                array3[i] = (i + 2*baseSize) * 10 
        if ( alternate == 1):                
            for i in range(len(array1)):
                array1[i] = i * 10 + 10 * math.sin(math.radians(i*10))
                array2[i] = (i + baseSize + 2 * math.sin(math.radians(i*10))) * 10
                array3[i] = (i + 2*baseSize + 3 * + 2 * math.sin(math.radians(i*10))) * 10
            #print( array1)
        array1[5 * 6] = ilwis.Const.iUNDEF # pos 0,2 is undefined            
        rc.array2raster(array1, 0)            
        rc.array2raster(array2, 1)
        rc.array2raster(array3, 2)                                

        return rc  

    def createEmptyTestTable(self):
        tbl = ilwis.Table()
        tbl.addColumn("ints", "integer")
        tbl.addColumn("floats", "value")
        cdef = ilwis.ColumnDefinition("items", self.createThematicDomain(), 2)
        tbl.addColumn(cdef)
        tbl.addColumn("strings1","text")
        tbl.addColumn("strings2","text")
        tbl.addColumn("strings3","text")
        tbl.addColumn("intt", "integer")

        return tbl

    def createTestTable(self):
        tbl = self.createEmptyTestTable()

        tbl.setCell("ints", 0 , 22)
        tbl.setCell("floats", 0, 34.987)
        tbl.setCell("items", 0, "stone")
        tbl.setCell("strings1", 0, "aap")
        tbl.setCell("strings2", 0, "100")
        tbl.setCell("strings3", 0, "houses")

        tbl.setCell("ints", 1 , 72)
        tbl.setCell("floats", 1, 114.6)
        tbl.setCell("items", 1, "water")
        tbl.setCell("strings1", 1, "noot")
        tbl.setCell("strings2", 1, "200")
        tbl.setCell("strings3", 1, "water")

        tbl.setCell("ints", 2 , 190)
        tbl.setCell("floats", 2, 13.6)
        tbl.setCell("items", 2, "houses")
        tbl.setCell("strings1", 2, "mies")
        tbl.setCell("strings2", 2, "300")
        tbl.setCell("strings3", 2, "houses")

        tbl.setCell("ints", 3 , 77)
        tbl.setCell("floats", 3, 10.12)
        tbl.setCell("items", 3, "houses") 
        tbl.setCell("strings1", 3, "wim")  
        tbl.setCell("strings2", 3, "400") 
        tbl.setCell("strings3", 3, "stone")

        tbl.setCell("ints", 4 , 309)
        tbl.setCell("floats", 4, 40.12)
        tbl.setCell("items", 4, "houses")  
        tbl.setCell("strings1", 4, "zus")  
        tbl.setCell("strings2", 4, "500") 
        tbl.setCell("strings3", 4, "notvalid") 

        tbl.setCell("ints", 5 , 309)
        tbl.setCell("floats", 5, 477.23)
        tbl.setCell("items", 5, "water")                
        tbl.setCell("strings1", 5, "jet")  
        tbl.setCell("strings2", 5, "600")  
        tbl.setCell("strings3", 5, "grass") 

        return tbl

    def createKeyedTestTable(self):
        tbl = self.createEmptyTestTable()

        #create attribute table; items column must be the key so no duplicates allowed
        tbl.setCell("ints", 0 , 22)
        tbl.setCell("floats", 0, 34.987)
        tbl.setCell("items", 0, "stone")
        tbl.setCell("strings1", 0, "aap")
        tbl.setCell("strings2", 0, "100")
        tbl.setCell("strings3", 0, "houses")

        tbl.setCell("ints", 1 , 72)
        tbl.setCell("floats", 1, 114.6)
        tbl.setCell("items", 1, "water")
        tbl.setCell("strings1", 1, "noot")
        tbl.setCell("strings2", 1, "200")
        tbl.setCell("strings3", 1, "water")

        tbl.setCell("ints", 2 , 190)
        tbl.setCell("floats", 2, 13.6)
        tbl.setCell("items", 2, "houses")
        tbl.setCell("strings1", 2, "mies")
        tbl.setCell("strings2", 2, "300")
        tbl.setCell("strings3", 2, "houses")        

        tbl.setCell("ints", 3 , 190)
        tbl.setCell("floats", 3, 13.6)
        tbl.setCell("items", 3, "grass")
        tbl.setCell("strings1", 3, "mies")
        tbl.setCell("strings2", 3, "300")
        tbl.setCell("strings3", 3, "houses")

        return tbl   

    def createFeatureCoverage(self):
        fcNew = ilwis.FeatureCoverage()
        csy = ilwis.CoordinateSystem("code=epsg:4326") # create coordinate system
        fcNew.setCoordinateSystem(csy)
        fcNew.setEnvelope(ilwis.Envelope('10 30 40 70')) 

        fcNew.addAttribute("ints", "integer")
        fcNew.addAttribute("floats", "value")
        cdef = ilwis.ColumnDefinition("items", self.createThematicDomain(), 2)
        fcNew.addAttribute(cdef)
        fcNew.addAttribute("strings1","text")
        fcNew.addAttribute("strings2","text")

        feature = fcNew.newFeature('Polygon((35.9 36.5, 35.9 38.6, 40.5 38.6, 40.5 36.5, 35.9 36.5))')
        feature.setAttribute('items', 'grass')
        feature.setAttribute('ints', 120)
        feature.setAttribute('floats', 23.89)
        feature.setAttribute('strings1', 'Aap')
        feature.setAttribute('strings2', '300')

        feature = fcNew.newFeature('Polygon((34.7 36.0, 34.7 38.0, 40.0 38.0, 40.0 36.0, 34.7 36.0))')
        feature.setAttribute('items', 'houses')
        feature.setAttribute('ints', 20)
        feature.setAttribute('floats', 1020.67)
        feature.setAttribute('strings1', 'Noot')
        feature.setAttribute('strings2', '150')

        feature = fcNew.newFeature('Polygon((24.7 56.0, 27.7 38.0, 31.0 38.0, 27.0 46.0, 24.7 56.0))')
        feature.setAttribute('items', 'water')
        feature.setAttribute('ints', 1120)
        feature.setAttribute('floats', 0.58)
        feature.setAttribute('strings1', '?')
        feature.setAttribute('strings2', '450')

        feature = fcNew.newFeature('Point(25 60)')
        feature.setAttribute('items', 'grass')
        feature.setAttribute('ints', 7120)
        feature.setAttribute('floats', 2.58)
        feature.setAttribute('strings1', 'no idead')
        feature.setAttribute('strings2', '1020')

        return fcNew

    def arrayValues(self, rc): # this is purely voor viewing the content of a test small map
        print('\n')
        xs = rc.size().xsize
        ys = rc.size().ysize
        arr = np.fromiter(ilwis.PixelIterator(rc), dtype=np.float64)
        for y in range(ys):
            print('\n')
            for x in range(xs):
                print(arr[y * xs + x], end=' ')
        print('\n')           


    


           


