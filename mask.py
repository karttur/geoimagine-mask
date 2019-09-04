'''
Created on 9 Nov 2018

@author: thomasgumbricht
'''


import numpy as np
from geoimagine.ktnumba import SingleMask, AddToMask, SetMask
#import garbage collect
import gc

def MultiBandMasking(srcDatumLayerD, dstDatumLayerD):
    maskNull = 1
    dstComp = list(dstDatumLayerD.keys())[0]
    MASK = np.zeros_like( dstDatumLayerD[dstComp].layer.NPBAND ) 
    for srcComp in srcDatumLayerD:
        '''
        print ('srcComp',srcComp)
        print (srcDatumLayerD[srcComp].layer.NPBAND)
        '''
        srcNull = srcDatumLayerD[srcComp].layer.cellnull

        MASK = AddToMask(MASK, srcDatumLayerD[srcComp].layer.NPBAND, maskNull, srcNull)
    for dstComp in dstDatumLayerD:
        outNull = dstDatumLayerD[dstComp].comp.cellnull
        dstDatumLayerD[dstComp].layer.NPBAND = SetMask(MASK, dstDatumLayerD[dstComp].layer.NPBAND, maskNull, outNull)
        '''
        dstLayerD[dstComp].layer.NPBAND[maskA == 1] = outNull   
        maskA[(srcLayerD[srcKey].layer.NPBAND == inNull) | (maskA == 1)] = 1
        for dstComp in dstLayerD:
            outNull = dstLayerD[dstComp].comp.cellnull
            dstLayerD[dstComp].layer.NPBAND[maskA == 1] = outNull
        '''
        
def SingleBandMasking(srcLayer, dstLayer):
    maskNull = 1
    srcNull = srcLayer.layer.cellnull
    dstNull = dstLayer.comp.cellnull
    MASK = np.zeros_like( srcLayer.layer.NPBAND ) 
    MASK = SingleMask(MASK, srcLayer.layer.NPBAND, maskNull, srcNull)
    MASK[(srcLayer.layer.NPBAND == srcNull)] = maskNull
    dstLayer.layer.NPBAND = SetMask(MASK, dstLayer.layer.NPBAND, maskNull, dstNull)
    
class ProcessMasking:
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        if self.process.proc.processid == 'createstaticmaskancillary':
            self._CreateStaticMask()
        elif self.process.proc.processid == 'applystaticmaskancillary':
            self._ApplyStaticMask()
        
    

    def _CreateStaticMask(self):
 
        for locus in self.process.srcLayerD:       
            self._GetDatumComp(locus)
            firstlayer = True
            if self.process.dstLayerD[locus][self.dstDatum][self.dstcomp]._Exists() and not self.process.overwrite:
                self.session._InsertLayer(self.process.dstLayerD[locus][self.dstDatum][self.dstcomp], self.process.overwrite, self.process.delete)
                print ('mask ready',self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].FPN)
                continue 
            
            dstNull = self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].comp.cellnull
            print ('dstNull', dstNull)
            if dstNull == 1:
                exit('Mask values can not be = 1')
            if self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].comp.celltype.lower() not in ['byte','uint8']:
                print (self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].comp.celltype)
                exit('Mask celltype must be byte/uint8')   
            for srcDatum in self.srcDatumL:
                if not (self.process.srcLayerD[locus][srcDatum][self.srccomp]):
                    print ('Src composition missing',srcDatum)
                    continue
                self.process.srcLayerD[locus][srcDatum][self.srccomp].ReadRasterLayer() 
                if firstlayer:
                    firstDatum = srcDatum
                    srcNull = self.process.srcLayerD[locus][srcDatum][self.srccomp].layer.cellnull
                    MASK = np.ones_like( self.process.srcLayerD[locus][srcDatum][self.srccomp].layer.NPBAND, dtype=np.uint8 )
                    MASK[(self.process.srcLayerD[locus][srcDatum][self.srccomp].layer.NPBAND == srcNull)] = dstNull
                    firstlayer = False
                else:
                    srcNull = self.process.srcLayerD[locus][srcDatum][self.srccomp].layer.cellnull
                    print ('srcNull',srcNull)
                    if self.process.params.maskanynull:
                        pass
                        #MASK[(self.process.srcLayerD[locus][srcDatum][self.srccomp].layer.NPBAND == srcNull)] = dstNull
                    else:
                        pass
                        #MASK[(self.process.srcLayerD[locus][srcDatum][self.srccomp].layer.NPBAND != srcNull)] = 1
                    #Close the layer
                    self.process.srcLayerD[locus][srcDatum][self.srccomp] = None
                    #garbage collect to free memory
                    gc.collect()
              
            self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].layer = lambda:None
            self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].layer.NPBAND = MASK
            #copy the geoformat from the src layer
            self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][firstDatum][self.srccomp].layer)
            #write the results
            self.process.dstLayerD[locus][self.dstDatum][self.dstcomp].CreateDSWriteRasterArray()
            #Register the layer
            self.session._InsertLayer(self.process.dstLayerD[locus][self.dstDatum][self.dstcomp], self.process.overwrite, self.process.delete)

    def _ApplyStaticMask(self):
 
        for locus in self.process.srcLayerD:
            self._GetStaticMask(locus)  
            for datum in self.process.srcLayerD[locus]:
                if self.process.dstLayerD[locus][datum][self.srccomp]._Exists() and not self.process.overwrite:
                    self.session._InsertLayer(self.process.dstLayerD[locus][datum][self.srccomp], self.process.overwrite, self.process.delete)
                    print ('layer already masked',self.process.dstLayerD[locus][datum][self.srccomp].FPN)
                    continue
                
                self.process.srcLayerD[locus][datum][self.srccomp].ReadRasterLayer() 
                srcNull = self.process.srcLayerD[locus][datum][self.srccomp].layer.cellnull
                BAND = self.process.srcLayerD[locus][datum][self.srccomp].layer.NPBAND
                
                BAND[self.MASK == self.maskNull] = srcNull
                self.process.dstLayerD[locus][datum][self.srccomp].layer = lambda:None
                self.process.dstLayerD[locus][datum][self.srccomp].layer.NPBAND = BAND
                #copy the geoformat from the src layer
                self.process.dstLayerD[locus][datum][self.srccomp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][self.srccomp].layer)
                #write the results
                self.process.dstLayerD[locus][datum][self.srccomp].CreateDSWriteRasterArray()
                #Register the layer
                self.session._InsertLayer(self.process.dstLayerD[locus][datum][self.srccomp], self.process.overwrite, self.process.delete)

                self.process.srcLayerD[locus][datum][self.srccomp] = None
                self.process.dstLayerD[locus][datum][self.srccomp] = None
                #garbage collect to free memory
                gc.collect()
            self._CloseMask()
        
    def _GetDatumComp(self,locus):
        self.srcDatumL = []
        dstDatumL = []
        for srcdatum in self.process.srcLayerD[locus]:
            self.srcDatumL.append(srcdatum)
            for srccomp in self.process.srcLayerD[locus][srcdatum]:
                self.srccomp = srccomp
        for dstdatum in self.process.dstLayerD[locus]:
            dstDatumL.append(dstdatum)
            for dstcomp in self.process.dstLayerD[locus][dstdatum]:
                self.dstcomp = dstcomp
        if len(dstDatumL) != 1:
            exitstr = 'ERORR in ProcessMask (mask.mask'
            exit(exitstr)
        else:
            self.dstDatum = dstDatumL[0]
            
    def _GetStaticMask(self,locus):
        datum = list(self.process.srcLayerD[locus].keys())[0]
        for srccomp in self.process.srcLayerD[locus][datum]:
            if self.process.srcLayerD[locus][datum][srccomp].comp.id == 'mask':
                self.mask = srccomp
            elif self.process.srcLayerD[locus][datum][srccomp].comp.id == 'layer':
                self.srccomp = srccomp
            else:
                exitstr = 'To apply a static mask the ids from the comps must be mask and layer'
                exit(exitstr)
        #Open the static mask to apply on all the layers in this locus
        self._OpenMask(locus,datum,self.mask)
                
    def _OpenMask(self,locus,datum,comp):
        self.process.srcLayerD[locus][datum][comp].ReadRasterLayer()
        self.MASK = self.process.srcLayerD[locus][datum][comp].layer.NPBAND
        self.maskNull = self.process.srcLayerD[locus][datum][comp].layer.cellnull
        
    def _CloseMask(self):
        self.MASK = None

        
                


                

        
        
