from rootpy.io import root_open
from root_numpy import tree2array
import numpy as np
import math as m
import pprint
from shapely.geometry import Polygon
import waferGeometry

f = root_open("test_triggergeom.root")
cells_tree = f.Get("hgcaltriggergeomtester/TreeCells") 
TC_tree = f.Get("hgcaltriggergeomtester/TreeTriggerCells")

cells_wafer_info = tree2array(cells_tree, branches = ['id','wafertype','wafer','layer','subdet','zside','x','y','cell','waferrow','wafercolumn'])
cells_tc_info = tree2array(TC_tree, branches = ['triggercell', 'c_id', 'c_cell', 'wafer', 'layer', 'subdet', 'zside'])
f.close()

#########################################################
## FUNCTIONS ##
def ExtractMappingCoordinates(koordx, koordy, d, x0, y0):
    distx = koordx - x0
    disty = koordy - y0
    correction_x = disty * m.tan(m.radians(30))
    map_i = (distx - correction_x)/(2*d)
    map_j_nominator = (float)(2*((float)(disty))*m.sqrt(3))/3
    map_j = map_j_nominator /(2*d)
    map_coord_i_j = "(" + str(round(map_i)) + "," + str(round(map_j)) + ")"
    map_coord_i_j_new = map_coord_i_j.replace('-0.', '0.')
    return map_coord_i_j_new

def ReconstructCellsVerticesList(wafertype, cell, koordx, koordy):
    wafer_geometry = waferGeometry.smallCellWafer if wafertype==1 else waferGeometry.largeCellWafer
    edge = wafer_geometry['cell_corner_size']
    xs = []
    ys = []
    diameter = edge/m.tan(m.radians(30))
    centerx = float(koordx)
    centery = float(koordy)
    # Shift center for half cells or corner cells
    if cell in wafer_geometry['half_cells_edge_left']:
        centerx -=  2.*m.sqrt(3)*edge/9.
    elif cell in wafer_geometry['half_cells_edge_topleft']:
        centerx -=  2.*m.sqrt(3)*edge/9.*m.cos(m.radians(60))
        centery +=  2.*m.sqrt(3)*edge/9.*m.sin(m.radians(60))
    elif cell in wafer_geometry['half_cells_edge_topright']:
        centerx +=  2.*m.sqrt(3)*edge/9.*m.cos(m.radians(60))
        centery +=  2.*m.sqrt(3)*edge/9.*m.sin(m.radians(60))
    elif cell in wafer_geometry['half_cells_edge_bottomright']:
        centerx +=  2.*m.sqrt(3)*edge/9.*m.cos(m.radians(60))
        centery -=  2.*m.sqrt(3)*edge/9.*m.sin(m.radians(60))
    elif cell in wafer_geometry['half_cells_edge_bottomleft']:
        centerx -=  2.*m.sqrt(3)*edge/9.*m.cos(m.radians(60))
        centery -=  2.*m.sqrt(3)*edge/9.*m.sin(m.radians(60))
    elif cell in wafer_geometry['half_cells_edge_right']:
        centerx +=  2.*m.sqrt(3)*edge/9.
    x = centerx - diameter/2.
    y = centery - edge/2.
    for angle in range(0, 360, 60):
        y += m.cos(m.radians(angle)) * edge
        x += m.sin(m.radians(angle)) * edge
        xs.append(x)
        ys.append(y)
    # Remove corners for half cells and corner cells
    if cell in wafer_geometry['half_cells_edge_left']:
        del xs[0]; del ys[0]
        del xs[-1]; del ys[-1]
    elif cell in wafer_geometry['half_cells_edge_topleft']:
        del xs[0]; del ys[0]
        del xs[0]; del ys[0]
    elif cell in wafer_geometry['half_cells_edge_topright']:
        del xs[1]; del ys[1]
        del xs[1]; del ys[1]
    elif cell in wafer_geometry['half_cells_edge_bottomright']:
        del xs[3]; del ys[3]
        del xs[3]; del ys[3]
    elif cell in wafer_geometry['half_cells_edge_bottomleft']:
        del xs[4]; del ys[4]
        del xs[4]; del ys[4]
    elif cell in wafer_geometry['half_cells_edge_right']:
        del xs[2]; del ys[2]
        del xs[2]; del ys[2]
    # Close the cell
    xs.append(xs[0])
    ys.append(ys[0])
    return xs,ys

def GetCellIDsFromWafer(cells_wafer_info, waferID, subdet):
    wafer_cellIDs = []
    for i in range(cells_wafer_info.shape[0]):
        if cells_wafer_info[i][2] == waferID:
            if cells_wafer_info[i][3] == 1:
                if cells_wafer_info[i][4] == subdet and cells_wafer_info[i][5] == -1:
                    elem = list([cells_wafer_info[i][0],(cells_wafer_info[i][6],cells_wafer_info[i][7]),cells_wafer_info[i][8]])
                    if elem not in wafer_cellIDs:
                        wafer_cellIDs.append(elem)
    print "min coordinates cell is " + str(min(wafer_cellIDs))
    return wafer_cellIDs

def GetMappingCellIDsWaferRing(wafer_ring_cellIDs, d, central_waferID):
    cellIDs_map = {}
    null_cell = 0
    temp = []
    for wafer, cells in wafer_ring_cellIDs.iteritems():
        if wafer == central_waferID:
            null_cell = min(cells)
            #print "referent cell is " + str(null_cell)
    for wafer, cells in wafer_ring_cellIDs.iteritems():
        temp = []
        for i in range(len(cells)):
            map_tuple = ExtractMappingCoordinates(cells[i][1][0], cells[i][1][1], d, null_cell[1][0], null_cell[1][1])
            temp.append(list([cells[i][0], map_tuple, cells[i][2]]))
        cellIDs_map[wafer] = temp
    return cellIDs_map

def ExtractNeighbourCellsRing(wafer_dict, neighbours, waferID, wafer_cellsIDs, wafer_cells_map, wafer_types_map):
    neighbour_map = {}
    temp = []
    vals_cellsxy = wafer_cellsIDs.values()
    vals = wafer_dict.values()
    mapa_cellsxy = [val for sublist in vals_cellsxy for val in sublist]
    mapa = [val for sublist in vals for val in sublist]
    for i in range (len(mapa)):
        if mapa[i] in wafer_dict[waferID]:
            trenutni_tuple_x = mapa[i][1][mapa[i][1].index("(") + 1:mapa[i][1].rindex(",")]
            trenutni_tuple_y = mapa[i][1][mapa[i][1].index(",") + 1:mapa[i][1].rindex(")")]
            temp = []
            for j in range (len(mapa)):
                if mapa[i] != mapa[j]:
                    if wafer_types_map[wafer_cells_map[mapa[i][0]]] == wafer_types_map[wafer_cells_map[mapa[j][0]]]:
                        for k in range(7):
                            searched_tuple = (neighbours[k][0], neighbours[k][1])
                            x = mapa[j][1][mapa[j][1].index("(") + 1:mapa[j][1].rindex(",")]
                            y = mapa[j][1][mapa[j][1].index(",") + 1:mapa[j][1].rindex(")")]
                            if float(x) == (float(trenutni_tuple_x) + searched_tuple[0]) and float(y) == (float(trenutni_tuple_y) + searched_tuple[1]):
                                temp.append(mapa[j])
                                neighbour_map[str(mapa[i])] = temp
                    else:
                        for k in range(len(mapa_cellsxy)):
                            if mapa_cellsxy[k][0] == mapa[i][0]:
                                cell_x = mapa_cellsxy[k][1][0]
                                cell_y = mapa_cellsxy[k][1][1]
                            elif mapa_cellsxy[k][0] == mapa[j][0]:
                                neighbor_x = mapa_cellsxy[k][1][0]
                                neighbor_y = mapa_cellsxy[k][1][1]
                        wafertype_central = wafer_types_map[wafer_cells_map[mapa[i][0]]]
                        x_central, y_central = ReconstructCellsVerticesList(wafertype_central, mapa[i][2], cell_x, cell_y)
                        cell1 = Polygon(zip(x_central, y_central))
                        wafertype_neighbor = wafer_types_map[wafer_cells_map[mapa[j][0]]]
                        x_neighbor, y_neighbor = ReconstructCellsVerticesList(wafertype_neighbor, mapa[j][2], neighbor_x, neighbor_y)
                        cell2 = Polygon(zip(x_neighbor, y_neighbor))
                        EPS = 0.05
                        if cell1.intersects(cell2) or cell1.buffer(EPS).intersects(cell2):
                            temp.append(mapa[j])
                            neighbour_map[str(mapa[i])] = temp
    return neighbour_map

def ExtractWaferInfoLayer1(cells_wafer_info, subdet):
    wafer_info = []
    wafers_with_types = {}
    first_layer = 1
    for i in range(cells_wafer_info.shape[0]):
        if cells_wafer_info[i][4] == subdet and cells_wafer_info[i][5] == -1:
            if cells_wafer_info[i][3] == first_layer:
                element = list([cells_wafer_info[i][2],(cells_wafer_info[i][9],cells_wafer_info[i][10])])
                if cells_wafer_info[i][1] == -1:
                    wafers_with_types[cells_wafer_info[i][2]] = 0
                else:
                    wafers_with_types[cells_wafer_info[i][2]] = cells_wafer_info[i][1]
                if element not in wafer_info:
                    wafer_info.append(element)
    return wafer_info, wafers_with_types

def ExtractWaferInfo(cells_wafer_info, subdet):
    wafer_info = []
    last_layer = 12
    if subdet == 3 : last_layer = 27
    for i in range(cells_wafer_info.shape[0]):
        if cells_wafer_info[i][4] == subdet and cells_wafer_info[i][5] == -1:
            if cells_wafer_info[i][3] == 1 or cells_wafer_info[i][3] == last_layer:
                element = list([cells_wafer_info[i][2],(cells_wafer_info[i][9],cells_wafer_info[i][10])])
                if element not in wafer_info:
                    wafer_info.append(element)
    return wafer_info

def ExtractCentralWafers(wafer_neighbors_map, wafer_types_map):
    central_wafers_with_config = []
    configuration_list = [] 
    for k in wafer_neighbors_map.keys():
        missing_neighbor = False
        config_string = str(wafer_types_map[k])
        for neighbor in wafer_neighbors_map[k]:
            if neighbor != -1:
                config_string += str(wafer_types_map[neighbor])
            else:
                missing_neighbor = True
                config_string += str(wafer_types_map[k])
        if config_string not in configuration_list and missing_neighbor == False:
            configuration_list.append(config_string)
            central_wafers_with_config.append((config_string, k)) 
    return central_wafers_with_config

def ExtractNeighbourWafers(wafer_info, neighbours):
    neighbours_map = {}
    neighbours_map_clean = {}
    temp = []
    temp2 = []
    for i in range(len(wafer_info)):
        waferID = wafer_info[i][0]
        wafer_x = wafer_info[i][1][0]
        wafer_y = wafer_info[i][1][1]
        temp = []
        temp2 = [0, 0, 0, 0, 0, 0]
        test_seq = [0, 0, 0, 0, 0, 0]
        for j in range(len(wafer_info)):
            if wafer_info[j] != waferID:
                for k in range(len(neighbours)):
                    neighbour_tuple = (neighbours[k][0], neighbours[k][1])
                    x = wafer_info[j][1][0]
                    y = wafer_info[j][1][1]
                    if float(x) == (float(wafer_x) + neighbour_tuple[0]) and float(y) == (float(wafer_y) + neighbour_tuple[1]):
                        test_seq[k] = 1
                        temp.append(wafer_info[j])
                        temp2[k] = wafer_info[j][0]
                        neighbours_map[str(wafer_info[i])] = temp
                        neighbours_map_clean[waferID] = temp2
        for t in range(len(test_seq)):
            if test_seq[t] == 0:
                neighbours_map_clean[waferID][t] = -1
    return neighbours_map, neighbours_map_clean

def ExtractTCcells(cells_tc_info, waferID, subdet):
    TC_cells_map = {}
    for i in range(cells_tc_info.shape[0]):
        temp = []
        if cells_tc_info[i][4] == 1:
            if cells_tc_info[i][3] == waferID and cells_tc_info[i][5] == subdet and cells_tc_info[i][6] == -1:
                temp.append(cells_tc_info[i][1].tolist()) 
                temp.append(cells_tc_info[i][2].tolist())
                TC_cells_map[cells_tc_info[i][0]] = temp
    return TC_cells_map

def invertTCmap(d):
    return dict( (v,k) for k in d for v in d[k][0] )

def invertWaferMap(d):
    return dict( (v,k) for k in d for v in d[k] )

def ExtractNeighborTCs(TC_cell_map_inv, wafer_cells_map_inv, cell_neighbor_map, relative_waferIDs, wafer_cells_map, wafer_types_map):
    TC_neighbor_map = {}
    TC_neighbor_map_final = {}
    TC_center = TC_neighbor = 0
    for key, cell_neighbors in cell_neighbor_map.iteritems():
        temp = []
        cellID = key[1:key.index(',')]
        TC_center = TC_cell_map_inv[int(cellID)]
        for n in cell_neighbors:
            neighborID = n[0]
            TC_neighbor = TC_cell_map_inv[neighborID]
            TC_neighbor_wafer = wafer_cells_map_inv[neighborID]
            tup = (TC_neighbor, relative_waferIDs[TC_neighbor_wafer])
            if TC_center != TC_neighbor and tup not in temp:
                temp.append(tup)
            if TC_center == TC_neighbor and tup not in temp and wafer_types_map[wafer_cells_map[int(cellID)]] != wafer_types_map[wafer_cells_map[neighborID]]:
                temp.append(tup)
        TC_neighbor_map.setdefault(TC_center,[]).append(temp)
    for k, v in TC_neighbor_map.iteritems():
        values_union = set().union(*v)
        TC_neighbor_map_final[k] = list(values_union)
    return TC_neighbor_map_final

def ExtractWaferRing(neighbours_wafer_map, wafer):
    key_list = neighbours_wafer_map.keys()
    central_wafer_ring = []
    for k in key_list:
        ID = k[1:k.index(",")]
        if ID == str(wafer):
            central_wafer_ring = neighbours_wafer_map[k]
    return central_wafer_ring

def ExtractCellIDsWaferRing(cells_wafer_info, central_wafer, central_wafer_ring, wafer_type, subdet):
    wafer_ring_cellsIDs = {}
    if wafer_type == -1:
        print "\n" + "Extract cellIDs from large central wafer ID " + str(central_wafer)
    else:
        print "\n" + "Extract cellIDs from small central wafer ID " + str(central_wafer)
    wafer_cellIDs = GetCellIDsFromWafer(cells_wafer_info, central_wafer, subdet)
    wafer_ring_cellsIDs[central_wafer] = wafer_cellIDs
    print "\n" + "Extract cellIDs from central wafer ring "
    for wafer in (central_wafer_ring):
        waferID = wafer[0]
        wafer_cellIDs = GetCellIDsFromWafer(cells_wafer_info, waferID, subdet)
        wafer_ring_cellsIDs[waferID] = wafer_cellIDs
    return wafer_ring_cellsIDs

def ExtractTCCellsWaferRing(cells_tc_info, central_wafer, central_wafer_ring, subdet):
    wafer_ring_TC_cellsIDs = {}
    print "\n" + "Extract TC->cellIDs map from central wafer " + str(central_wafer)
    wafer_TC_cellIDs = ExtractTCcells(cells_tc_info, central_wafer, subdet)
    wafer_ring_TC_cellsIDs[central_wafer] = wafer_TC_cellIDs
    print "\n" + "Extract TC->cellIDs map from central wafer ring "
    for wafer in (central_wafer_ring):
        waferID = wafer[0]
        wafer_TC_cellIDs = ExtractTCcells(cells_tc_info, waferID, subdet)
        wafer_ring_TC_cellsIDs[waferID] = wafer_TC_cellIDs
    return wafer_ring_TC_cellsIDs

def InvertMaps(wafer_ring_TC_cells_dict):                                                                               
    ## make dict cid:tcid ##                                                                                             
    TC_map_inv = []
    wafer_TC_cellIDs_inv = {}
    for wafID, TC_map in wafer_ring_TC_cells_dict.iteritems():
        TC_map_inv.append(invertTCmap(TC_map))
        for d in TC_map_inv:
            wafer_TC_cellIDs_inv.update(d)                                                                             
    ## make dict cid:waferid ##                                                                                             
    wafer_cells_map = {}
    wafer_cells_map_inv = {}
    for wafID, TC_map in wafer_ring_TC_cells_dict.iteritems():
        temp = []
        for tc, cell_list in TC_map.iteritems():
            temp.append(cell_list[0])
        flat_list = [item for sublist in temp for item in sublist]
        wafer_cells_map[wafID] = flat_list
    wafer_cells_map_inv = invertWaferMap(wafer_cells_map)
    return wafer_TC_cellIDs_inv, wafer_cells_map_inv

def ChangeWaferNeighborsToRelative(neighborIDs, waferID):
    relative_waferIDs = {}
    counter = 0
    relative_waferIDs[waferID] = counter
    for ID in neighborIDs:
        counter += 1
        relative_waferIDs[ID] = counter
    return relative_waferIDs

def UpdateTCMapWithConfiguration(neighbor_map, config):
    neighbor_map_update = {}
    for key in neighbor_map:
        new_k = '(' + config + ',' + str(key) + ')'
        neighbor_map_update[new_k] = neighbor_map[key]
    return neighbor_map_update

def WriteToFile(map_large, map_small):
    file_name = "output_files/wafer_map.txt"
    first_column = 3
    with open(file_name, "w") as myfile:
        for key in sorted(map_large):
            myfile.write(str(first_column) + " " + str(key) + " ")
            for value in map_large[key]:
                myfile.write(str(value) + " ")
            myfile.write("\n")
        if type(first_column) != str:
            first_column += 1
        for key in sorted(map_small):
            myfile.write(str(first_column) + " " + str(key) + " ")
            for value in map_small[key]:
                myfile.write(str(value) + " ")
            myfile.write("\n")

def DivideWafers(interesting_central_wafers):
    list_large = []
    list_small = []
    divided_wafers = []
    for tup in interesting_central_wafers:
        if tup[0][0] == '0':
            list_large.append(tup)
        else:
            list_small.append(tup)
    divided_wafers.append(list_large)
    divided_wafers.append(list_small)
    return divided_wafers

def ShortenListFH(wafersFH, wafersEE):
    for ee in wafersEE:
        for fh in wafersFH:
            if ee[0] == fh[0]:
                wafersFH.remove(fh)
    return wafersFH

def TCNeighborFinder(divided_wafers_one_type, neighbours_wafer_map, map_l1, w_type, d_cell, neighbors, subdet, file_name, wafer_types_map_l1):
    TC_map = {}
    for element in divided_wafers_one_type:
        waferID = element[1]
        central_wafer_ring = ExtractWaferRing(neighbours_wafer_map, waferID)
        print "\n" + "large wafer(ID=" + str(waferID) + ") neighbors"
        print central_wafer_ring
        wafer_ring_cellsIDs = ExtractCellIDsWaferRing(cells_wafer_info, waferID, central_wafer_ring, w_type, subdet)
        wafer_ring_cellIDs_mapping = GetMappingCellIDsWaferRing(wafer_ring_cellsIDs, d_cell, waferID)
        wafer_ring_TC_cellsIDs = ExtractTCCellsWaferRing(cells_tc_info, waferID, central_wafer_ring, subdet)
        wafer_TC_cellIDs_inv, wafer_cells_map_inv = InvertMaps(wafer_ring_TC_cellsIDs)
        neighbours_cells_wafer_ring = ExtractNeighbourCellsRing(wafer_ring_cellIDs_mapping, neighbors, waferID, \
                                                                    wafer_ring_cellsIDs, wafer_cells_map_inv, wafer_types_map_l1)
        relative_waferIDs = ChangeWaferNeighborsToRelative(map_l1[waferID], waferID)
        TC_neighbor_map = ExtractNeighborTCs(wafer_TC_cellIDs_inv, wafer_cells_map_inv, neighbours_cells_wafer_ring, relative_waferIDs, \
                                                 wafer_cells_map_inv, wafer_types_map_l1)
        TC_neighbor_map_update = UpdateTCMapWithConfiguration(TC_neighbor_map, element[0])
        print "\n" + "TC NEIGHBOR RESULT with configurations"
        pp.pprint(TC_neighbor_map_update)
        with open(file_name, "a") as myfile:
            for key in sorted(TC_neighbor_map_update, key=lambda k: int(k[k.find(',')+1:k.find(')')])):
                myfile.write(str(key) + " ")
                for value in TC_neighbor_map_update[key]:
                    myfile.write(str(value) + " ")
                myfile.write("\n")

#########################################################
pp = pprint.PrettyPrinter(indent=4)
d_cell_large = 0.562273
a_large = (2 * d_cell_large * m.sqrt(3))/3
d_cell_small = 0.41233
a_small = (2 * d_cell_small * m.sqrt(3))/3

print "\n" + "CELL INFO"
print "large cell d is " + str(d_cell_large) + " and cell a is " + str(a_large)
print "small cell d is " + str(d_cell_small) + " and cell a is " + str(a_small)

######################################################
## extract waferIDs and their (row, col) position -> wafer_info ##
wafer_infoEE = ExtractWaferInfo(cells_wafer_info, 3)
wafer_infoFH = ExtractWaferInfo(cells_wafer_info, 4)

print "\n" + "WAFER INFO"
print "Num of wafers in 1st and last layer (EE part) is " + str(len(wafer_infoEE))
print "Num of wafers in 1st and last layer (FH part) is " + str(len(wafer_infoFH))

## get neighbouring wafers that are different for EE and FH part ##
neighbours_wafer = [(0, 2), (-1, 1), (-1, -1), (0, -2), (1, -1), (1, 1)]
neighbours_wafer_mapEE, mapEE = ExtractNeighbourWafers(wafer_infoEE, neighbours_wafer)
neighbours_wafer_mapFH, mapFH = ExtractNeighbourWafers(wafer_infoFH, neighbours_wafer)

## write neighbouring wafers information to file ##
WriteToFile(mapEE, mapFH)

#####################################################
## extract wafer info only on Layer1 (for configurations codewords extraction) ##
## extract waferIDs and their (row, col) position -> wafer_info_l1 ##
## extract waferIDs and their type -> wafer_types_map_l1 ##
wafer_infoEE_l1, wafer_types_mapEE_l1 = ExtractWaferInfoLayer1(cells_wafer_info, 3)
wafer_infoFH_l1, wafer_types_mapFH_l1 = ExtractWaferInfoLayer1(cells_wafer_info, 4)

## get neighbors wafers only layer1 ##
neighbours_wafer_mapEE_l1, mapEE_l1 = ExtractNeighbourWafers(wafer_infoEE_l1, neighbours_wafer)
neighbours_wafer_mapFH_l1, mapFH_l1 = ExtractNeighbourWafers(wafer_infoFH_l1, neighbours_wafer) 

## extract interesting central wafers on the layer1 to extract wafer configurations ##
interesting_central_wafersEE = ExtractCentralWafers(mapEE_l1, wafer_types_mapEE_l1)
interesting_central_wafersFH = ExtractCentralWafers(mapFH_l1, wafer_types_mapFH_l1)

print "\n" + "INTERESTING CENTRAL WAFERS EE"
pp.pprint(interesting_central_wafersEE)
print "len" + str(len(interesting_central_wafersEE))

ShortenListFH(interesting_central_wafersFH, interesting_central_wafersEE)

print "\n" + "INTERESTING CENTRAL WAFERS FH"
pp.pprint(interesting_central_wafersFH)
print "len" + str(len(interesting_central_wafersFH))

## divide interesting central wafers in two subgroups based on their type (small or large) #
divided_wafers_EE = DivideWafers(interesting_central_wafersEE)
divided_wafers_FH = DivideWafers(interesting_central_wafersFH) 

print "\n" + "LARGE AND SMALL DIVISION OF WAFERS EE"
pp.pprint(divided_wafers_EE)
print "\n" + "LARGE AND SMALL DIVISION OF WAFERS FH"
pp.pprint(divided_wafers_FH)

#####################################################
import time
t1 = time.time()
neighbors = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1), (0, 0)]
file_name = "output_files/TC_map.txt"
f = open(file_name, "w")
## TC neighbor finder algo call for all interesting wafers EE and FH ## 
TCNeighborFinder(divided_wafers_EE[0], neighbours_wafer_mapEE_l1, mapEE_l1, -1, d_cell_large, neighbors, 3, file_name, wafer_types_mapEE_l1)
TCNeighborFinder(divided_wafers_FH[0], neighbours_wafer_mapFH_l1, mapFH_l1, -1, d_cell_large, neighbors, 4, file_name, wafer_types_mapFH_l1)
TCNeighborFinder(divided_wafers_EE[1], neighbours_wafer_mapEE_l1, mapEE_l1, 1, d_cell_small, neighbors, 3, file_name, wafer_types_mapEE_l1)
TCNeighborFinder(divided_wafers_FH[1], neighbours_wafer_mapFH_l1, mapFH_l1, 1, d_cell_small, neighbors, 4, file_name, wafer_types_mapFH_l1)
print "\n" + "execution time is " + str((time.time()-t1)/60) + " min"
f.close()
