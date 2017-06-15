from rootpy.io import root_open
from rootpy.plotting import Canvas
from rootpy.interactive import wait
import numpy as np
import pprint
import ROOT

f = root_open("test_triggergeom_921_D17.root")
#f = root_open("test_triggergeom.root")                                                                                            
cells_tree = f.Get("hgcaltriggergeomtester/TreeCells")
TC_tree = f.Get("hgcaltriggergeomtester/TreeTriggerCells")

#################################                                                                                                   
## added for testing purpose (ONLY ONE WAFER)##                                                                                      
t = ROOT.TBrowser()
f_ = ROOT.TFile("test_triggergeom_921_D17.root")
cells_tree = f_.Get("hgcaltriggergeomtester/TreeCells")
'''cells_tree.Draw("y:x","layer==28 && subdet==3 && zside==1")                                                                      
cells_tree.SetMarkerColor(2)                                                                                                         
wafID = 121                                                                                                                          
parameter_line = "layer==28 && subdet==3 && zside==1 && wafer==" + str(wafID)                                                       
###print parameter_line                                                                                                              
cells_tree.Draw("y:x",parameter_line,"same")                                                                                        
wait(True)'''
############################################                                                                                         
#################TEST WHOLE MAP######################                                                                                
def GetEdgeWaferNeighbors(wafer_file_lines):
    edge_wafers_dict = {}
    for l in wafer_file_lines:
        line = l.split()
        if line[0] == '3' and '-1' in line[2:]:
            print line
            lista_bez_rupa = filter(lambda elem: elem != '-1', line[2:])
            edge_wafers_dict[line[1]] = lista_bez_rupa
    return edge_wafers_dict

## test written wafer neighbors in the file ##                                                                                      
wafer_file = open('wafer_map.txt','r')
wafer_file_lines = wafer_file.readlines()
edge_wafers_dict = GetEdgeWaferNeighbors(wafer_file_lines)
pp = pprint.PrettyPrinter(indent=4)
print "EDGE WAFERS"
pp.pprint(edge_wafers_dict)
wafer_file.close()

## visualize wafer neighbors ##                                                                                                      
cells_tree.Draw("y:x","layer==28 && subdet==3 && zside==1")
for wafID, waf_neighbors in edge_wafers_dict.iteritems():
    cells_tree.SetMarkerColor(2)
    #if int(wafID) < 10:
    #if wafID in ['5', '9', '147', '145', '79', '78', '262', '474', '475', '359', '358']:
    #if wafID in ['666', '665', '664', '663', '662', '661', '660', '499', '288', '4', '8', '548', '409', '408', '570', '263','121']:
    #if wafID in ['123', '59', '591', '590', '195', '312', '196', '628', '525', '3', '7', '523', '422', '449', '382', '383']:
    if wafID in ['100', '248', '171', '170', '289', '313', '2', '6', '500', '547', '569', '99', '222', '221', '659', '394']:
    #if len(waf_neighbors) == 4:
        parameter_line = "layer==28 && subdet==3 && zside==1 && wafer==" + str(wafID)
        print parameter_line
        cells_tree.Draw("y:x",parameter_line,"same")
        for n in waf_neighbors:
            cells_tree.SetMarkerColor(4)
            parameter_line = "layer==28 && subdet==3 && zside==1 && wafer==" + str(n)
            #print parameter_line
            cells_tree.Draw("y:x",parameter_line,"same")
wait(True)
#########################################                                                
