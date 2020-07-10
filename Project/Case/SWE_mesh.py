import pandas as pd
import numpy as np
import networkx as nx
import random
import re
import os
import struct
import json
import h5py
from osgeo import gdal, ogr
from shapely.geometry import shape, Polygon, Point
import geopandas as gpd


def getRasterVal(rb_array, x, y, x0, y0, dx, dy):
    px = int((x-x0)/dx)
    py = int((y-y0)/dy)
    height, width = rb_array.shape
    if (py > height-1) or (px > width-1):
        print("Range Exceeds")

    return rb_array[min(height-1, py)][min(width-1, px)]


def getBoundarySeg(poly_addr):
    # closeness of Boundary must be checked before
    driver = ogr.GetDriverByName("ESRI Shapefile")
    polyline_file = driver.Open(poly_addr, 0)
    polyline_layer = polyline_file.GetLayer()

    pt_list = []
    polyline_list = []

    for pl in polyline_layer:
        pl_json = json.loads(pl.ExportToJson())
        shp_pl = shape(pl_json["geometry"])
        x, y = shp_pl.coords.xy
        coordinates = list(zip(x, y))
        # append for coordinates list
        for cr in coordinates:
            if cr not in pt_list:
                pt_list.append(cr)
        # append for polyline list
        polyline_list.append([pt_list.index(cr) for cr in coordinates])

    return pt_list, polyline_list


def getPolygonSeg(poly_addr):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    polygon_file = driver.Open(poly_addr, 0)
    polygon_layer = polygon_file.GetLayer()

    pt_list = []
    polygon_list = []

    for polygon in polygon_layer:
        pg_json = json.loads(polygon.ExportToJson())
        shp_pg = shape(pg_json["geometry"])
        if shp_pg.type == "Polygon":
            # only points of outer ring are used
            x, y = shp_pg.exterior.coords.xy
            coordinates = list(zip(x, y))
            # append for coordinates list
            for cr in coordinates:
                if cr not in pt_list:
                    pt_list.append(cr)
            # append for polygon list
            polygon_list.append([pt_list.index(cr) for cr in coordinates])
        else:
            for p in shp_pg:
                x, y = p.exterior.coords.xy
                coordinates = list(zip(x, y))
                for cr in coordinates:
                    if cr not in pt_list:
                        pt_list.append(cr)
                polygon_list.append([pt_list.index(cr) for cr in coordinates])

    return pt_list, polygon_list


def getBCaddr(sg_id):
    if sg_id == 1:
        return "/outer_H.txt"
    else:
        return "No"


def shp2poly(outer_addr, building_addr, save_addr, bd_marker=1):
    '''
    Parameters
    ----------
    outer_addr : the address of outer boundary information .shp file
    building_addr: the address of inner building information .shp file
    save_addr : the address to save Triangle's .poly file
    bd_marker: the total number of boundary markers

    BC configuration
    [1] Outer Boundary: Free outflow, 1
    [2] Building Boundary: no-slip wall, 2
    '''
    outer_pt_info, outer_pl_info = getBoundarySeg(outer_addr)
    building_pt_info, building_pg_info = getPolygonSeg(building_addr)

    outer_pt_num = len(outer_pt_info)
    outer_sg_num = sum([len(pl)-1 for pl in outer_pl_info])
    building_pt_num = len(building_pt_info)
    building_sg_num = sum([len(pg)-1 for pg in building_pg_info])

    # Writing .poly file
    # [1] create the output file
    poly_file = open(save_addr, mode="a")
    # [2] write summary for vertices
    poly_file.write("# Part of vertices\n")
    pt_num = building_pt_num + outer_pt_num
    poly_file.write(str(pt_num) + " 2 0 " + str(bd_marker) + "\n")
    # write vertex of buildings
    for i in range(1, building_pt_num + 1):
        poly_file.write(str(i) + " " + " ".join(map(str, building_pt_info[i-1])) + " 2\n")
    # write vertex of outer boundary
    for i in range(1, outer_pt_num + 1):
        poly_file.write(str(i+building_pt_num) + " " + " ".join(map(str, outer_pt_info[i-1])) + " 1\n")
    # [3] write summary for segments
    poly_file.write("# Part of segments\n")
    sg_num = outer_sg_num + building_sg_num
    poly_file.write(str(sg_num) + " " + str(bd_marker) + "\n")
    print("The number of segments is " + str(sg_num))
    # write segment of buildings
    counter = 1
    for pg in building_pg_info:
        pg_num_i = len(pg)
        for i in range(0, pg_num_i-1):
            poly_file.write(str(counter) + " " + str(pg[i]+1) + " " + str(pg[i+1]+1) + " 2\n")
            counter += 1
    # write segment of outer boundary
    for pl in outer_pl_info:
        pl_num_i = len(pl)
        for i in range(0, pl_num_i-1):
            poly_file.write(str(counter) + " " + str(pl[i]+1+building_pt_num) + " " + str(pl[i+1]+1+building_pt_num) + " 1\n")
            counter += 1
    # [4] write summary for holes
    hole_num = len(building_pg_info)
    poly_file.write("# Part of holes\n")
    poly_file.write(str(hole_num) + "\n")
    print("The number of hole is " + str(hole_num))
    # write hole i with an inner point for the polygon
    h_count = 1
    for pg in building_pg_info:
        pg_coord = [building_pt_info[p] for p in pg]
        polygon_i = Polygon(pg_coord)
        while True:
            test_sample = random.sample(pg_coord, 3)
            test_centroid = np.average(test_sample, axis=0)
            if polygon_i.contains(Point(test_centroid)):
                poly_file.write(str(h_count) + " " + " ".join(map(str, test_centroid)) + "\n")
                break
        h_count += 1


def Tri2IUHM(node_addr, num_node, edge_addr, num_edge, ele_addr, num_ele, voro_addr, DEM_addr, LC_addr, save_addr):
    # open mesh file
    node_df = pd.read_table(node_addr, names=["ID", "x", "y", "Attribute"], sep="\s+", skiprows=lambda x: x in [0, num_node+1])
    edge_df = pd.read_table(edge_addr, names=["ID", "From", "To", "BD_marker"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    ele_df = pd.read_table(ele_addr, names=["ID", "N1", "N2", "N3"], sep="\s+", skiprows=lambda x: x in [0, num_ele+1])
    voro_df = pd.read_table(voro_addr, names=["ID", "E1", "E2", "x", "y"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    # open geo-Raster file
    DEM_ras = gdal.Open(DEM_addr)
    DEM_gt = DEM_ras.GetGeoTransform()
    DEM_rb = DEM_ras.GetRasterBand(1)
    DEM_rb_array = DEM_rb.ReadAsArray()
    DEM_mean = np.mean(np.ma.masked_where(DEM_rb_array < 0.0, DEM_rb_array))
    DEM_rb_array = np.where(DEM_rb_array < 0, DEM_mean, DEM_rb_array)
    LC_ras = gdal.Open(LC_addr)
    LC_gt = LC_ras.GetGeoTransform()
    LC_rb = LC_ras.GetRasterBand(1)
    LC_rb_array = LC_rb.ReadAsArray()
    # Node information
    node_id = np.linspace(1, num_node, num_node)
    node_x = np.array(node_df["x"])
    node_y = np.array(node_df["y"])
    # ---get the elevation for the node
    node_ele = np.array([getRasterVal(DEM_rb_array, node_x[i], node_y[i], DEM_gt[0], DEM_gt[3], DEM_gt[1], DEM_gt[5])
                         for i in range(0, num_node)])
    with h5py.File(os.path.join(save_addr, 'node.h5'), 'w') as hf:
        hf.create_dataset("NodeID", data=node_id, dtype=np.int64)
        hf.create_dataset("NodeX", data=node_x, dtype=np.float64)
        hf.create_dataset("NodeY", data=node_y, dtype=np.float64)
        hf.create_dataset("NodeEle", data=node_ele, dtype=np.float64)
    # Edge information
    edge_id = np.linspace(1, num_edge, num_edge)
    edge_from = np.array(edge_df["From"])
    edge_to = np.array(edge_df["To"])
    # ---set the boundary flag for the edge
    edge_bd = np.array(edge_df["BD_marker"])
    # ---specify the address of the boundary condition file for the edge
    edge_BDaddr = np.array([getBCaddr(edge_bd[i]).encode() for i in range(0, num_edge)])
    with h5py.File(os.path.join(save_addr, 'edge.h5'), 'w') as hf:
        hf.create_dataset("EdgeID", data=edge_id, dtype=np.int64)
        hf.create_dataset("EdgeF", data=edge_from, dtype=np.int64)
        hf.create_dataset("EdgeT", data=edge_to, dtype=np.int64)
        hf.create_dataset("BDFlag", data=edge_bd, dtype=np.int64)
        hf.create_dataset("BCsource", data=edge_BDaddr, dtype=h5py.special_dtype(vlen=str))
    # cell information
    cell_id = np.linspace(1, num_ele, num_ele)
    cellEdge = [[] for i in range(0, num_ele)]
    voro_e1 = np.array(voro_df["E1"])
    voro_e2 = np.array(voro_df["E2"])
    for i in range(0, num_edge):
        e1 = int(voro_e1[i])
        e2 = int(voro_e2[i])
        cellEdge[e1-1].append(i+1)
        if e2 != -1:
            cellEdge[e2-1].append(i+1)
    cellEdge = np.array(cellEdge)
    N1 = np.array(ele_df["N1"])
    N2 = np.array(ele_df["N2"])
    N3 = np.array(ele_df["N3"])
    # ---get the Manning's coefficient for the cell
    c_cx = np.array([(node_x[N1[i]-1]+node_x[N2[i]-1]+node_x[N3[i]-1])/3 for i in range(0, num_ele)])
    c_cy = np.array([(node_y[N1[i]-1]+node_y[N2[i]-1]+node_y[N3[i]-1])/3 for i in range(0, num_ele)])
    lc_n_dict = {10: 0.035, 20: 0.05, 30: 0.03, 50: 0.05, 60: 0.001, 80: 0.012, 90: 0.02}
    cell_n = np.array([lc_n_dict[getRasterVal(LC_rb_array, c_cx[i], c_cy[i], LC_gt[0], LC_gt[3], LC_gt[1], LC_gt[5])]
                        for i in range(0, num_ele)])
    lc_rc_dict = {10: 0.6, 20: 0.2, 30: 0.2, 50: 0.3, 60: 1.0, 80: 1.0, 90: 0.7}
    cell_rc = np.array([lc_rc_dict[getRasterVal(LC_rb_array, c_cx[i], c_cy[i], LC_gt[0], LC_gt[3], LC_gt[1], LC_gt[5])]
                        for i in range(0, num_ele)])

    h_ini = np.array([0.0 for i in range(0, num_ele)])
    u_ini = np.array([0.0 for i in range(0, num_ele)])
    v_ini = np.array([0.0 for i in range(0, num_ele)])
    with h5py.File(os.path.join(save_addr, 'cell.h5'), 'w') as hf:
        hf.create_dataset("CellID", data=cell_id, dtype=np.int64)
        hf.create_dataset("CellE1", data=cellEdge[:, 0], dtype=np.int64)
        hf.create_dataset("CellE2", data=cellEdge[:, 1], dtype=np.int64)
        hf.create_dataset("CellE3", data=cellEdge[:, 2], dtype=np.int64)
        hf.create_dataset("Celln", data=cell_n, dtype=np.float64)
        hf.create_dataset("Cellrc", data=cell_rc, dtype=np.float64)
        hf.create_dataset("h", data=h_ini, dtype=np.float64)
        hf.create_dataset("U", data=u_ini, dtype=np.float64)
        hf.create_dataset("V", data=v_ini, dtype=np.float64)


def prepareRainfall(rain_addr, dt_rain, num_cell, save_addr):
    f = open(rain_addr)
    counter = 0
    for line in f.readlines():
        rain_array = np.ones(num_cell) * float(line)
        with h5py.File(os.path.join(save_addr, 'rain_20180723.h5'), 'a') as hf:
            hf.create_dataset(str(dt_rain*counter), data=rain_array, dtype=np.float64)
        counter += 1
        print(dt_rain*counter)


'''
shp2poly(outer_addr="IUHM/Case_THU/THU_data/campusBoundaryClean.shp", building_addr="IUHM/Case_THU/THU_data/THU_building.shp", save_addr="IUHM/Case_THU/Case_THU.poly", bd_marker=1)
Tri2IUHM("IUHM/Case_THU/Case_THU.1.node", 236778, "IUHM/Case_THU/Case_THU.1.edge", 688250, "IUHM/Case_THU/Case_THU.1.ele", 450996,
         "IUHM/Case_THU/Case_THU.1.v.edge", "IUHM/Case_THU/THU_data/campusDEM_2m1.tif", "IUHM/Case_THU/THU_data/THU_LandCover_UTM50N_extent.tif", "IUHM/Case_THU/Case_Input")
prepareRainfall("IUHM/Case_THU/20180723/20180723.txt", 60, 450996, "IUHM/Case_THU/20180723")
'''