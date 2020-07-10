import pandas as pd
import numpy as np
import re
import os
import h5py
from osgeo import gdal, ogr
from shapely.geometry import shape, Polygon, Point


def prepareRainfall_zero(dt_rain, length, num_cell, save_addr):
    counter = 0
    for i in range(0, length):
        rain_array = np.zeros(num_cell)
        with h5py.File(os.path.join(save_addr, 'rain_zero.h5'), 'a') as hf:
            hf.create_dataset(str(dt_rain*counter), data=rain_array, dtype=np.float64)
        counter += 1
        print(dt_rain*counter)


# ************************* Stocker Dam Break Test *************************
def Stocker2poly(save_addr):
    # Writing .poly file
    # [1] create the output file
    poly_file = open(save_addr, mode="a")
    # [2] write summary for vertices
    poly_file.write("# Part of vertices\n")
    poly_file.write("6 2 0 1\n")
    # write vertices
    poly_file.write("1 0.0 0.0 1\n")
    poly_file.write("2 500.0 0.0 1\n")
    poly_file.write("3 1000.0 0.0 1\n")
    poly_file.write("4 1000.0 100.0 1\n")
    poly_file.write("5 500.0 100.0 1\n")
    poly_file.write("6 0.0 100.0 1\n")
    # [3] write summary for segments
    poly_file.write("# Part of segments\n")
    poly_file.write("7 1\n")
    poly_file.write("1 1 2 2\n")
    poly_file.write("2 2 3 2\n")
    poly_file.write("3 1 5 0\n")
    poly_file.write("4 3 4 1\n")
    poly_file.write("5 4 5 2\n")
    poly_file.write("6 5 6 2\n")
    poly_file.write("7 6 1 3\n")
    # [4] write summary for holes
    poly_file.write("# Part of holes\n")
    poly_file.write("0\n")

    poly_file.close()


def get_hIni_Stocker(x, y):
    if x <= 500:
        return 5.0
    else:
        return 0.2


def Tri2IUHM_Stocker(node_addr, num_node, edge_addr, num_edge, ele_addr, num_ele, voro_addr, save_addr):
    # open mesh file
    node_df = pd.read_table(node_addr, names=["ID", "x", "y", "Attribute"], sep="\s+", skiprows=lambda x: x in [0, num_node+1])
    edge_df = pd.read_table(edge_addr, names=["ID", "From", "To", "BD_marker"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    ele_df = pd.read_table(ele_addr, names=["ID", "N1", "N2", "N3"], sep="\s+", skiprows=lambda x: x in [0, num_ele+1])
    voro_df = pd.read_table(voro_addr, names=["ID", "E1", "E2", "x", "y"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    # Node information
    node_id = np.linspace(1, num_node, num_node)
    node_x = np.array(node_df["x"])
    node_y = np.array(node_df["y"])
    # ---get the elevation for the node
    node_ele = np.zeros(num_node)
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
    edge_BDaddr = np.array(["NoBCs".encode() for i in range(0, num_edge)])
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
    cell_n = np.array([0.0 for i in range(0, num_ele)])
    cell_rc = np.array([1.0 for i in range(0, num_ele)])
    h_ini = np.array([get_hIni_Stocker(c_cx[i], c_cy[i]) for i in range(0, num_ele)])
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


# ************************* Still Water over a 2D Bump Test *************************
def Still2poly(save_addr):
    # Writing .poly file
    # [1] create the output file
    poly_file = open(save_addr, mode="a")
    # [2] write summary for vertices
    poly_file.write("# Part of vertices\n")
    poly_file.write("4 2 0 1\n")
    # write vertices
    poly_file.write("1 0.0 0.0 1\n")
    poly_file.write("2 1.0 0.0 1\n")
    poly_file.write("3 1.0 1.0 1\n")
    poly_file.write("4 0.0 1.0 1\n")
    # [3] write summary for segments
    poly_file.write("# Part of segments\n")
    poly_file.write("4 1\n")
    poly_file.write("1 1 2 2\n")
    poly_file.write("2 2 3 2\n")
    poly_file.write("3 3 4 0\n")
    poly_file.write("4 4 1 2\n")
    # [4] write summary for holes
    poly_file.write("# Part of holes\n")
    poly_file.write("0\n")

    poly_file.close()


def get_DEM_Still(x, y):
    return max(0.0, 0.25 - 5 * ((x-0.5)**2 + (y-0.5)**2))


def get_hIni_Still(z):
    if 0.1 > z[2]:
        return 0.1 - (z[0]+z[1]+z[2])/3.0
    elif 0.1 > z[1]:
        return (0.1*0.1 + 0.1*z[2] - 0.3*z[0] - z[2]*z[1] + z[0]*z[1] + z[0]*z[0]) / 3.0 / (z[2]-z[0])
    elif 0.1 > z[0]:
        return (0.1-z[0])**3 / 3.0 / (z[1]-z[0]) / (z[2]-z[0])
    else:
        return 0.0


def Tri2IUHM_Still(node_addr, num_node, edge_addr, num_edge, ele_addr, num_ele, voro_addr, save_addr):
    # open mesh file
    node_df = pd.read_table(node_addr, names=["ID", "x", "y", "Attribute"], sep="\s+", skiprows=lambda x: x in [0, num_node+1])
    edge_df = pd.read_table(edge_addr, names=["ID", "From", "To", "BD_marker"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    ele_df = pd.read_table(ele_addr, names=["ID", "N1", "N2", "N3"], sep="\s+", skiprows=lambda x: x in [0, num_ele+1])
    voro_df = pd.read_table(voro_addr, names=["ID", "E1", "E2", "x", "y"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    # Node information
    node_id = np.linspace(1, num_node, num_node)
    node_x = np.array(node_df["x"])
    node_y = np.array(node_df["y"])
    # ---get the elevation for the node
    node_ele = np.array([get_DEM_Still(node_x[i], node_y[i]) for i in range(0, num_node)])
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
    edge_BDaddr = np.array(["NoBCs".encode() for i in range(0, num_edge)])
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
    cell_n = np.array([0.0 for i in range(0, num_ele)])
    cell_rc = np.array([1.0 for i in range(0, num_ele)])
    h_ini = np.array([get_hIni_Still(sorted([node_ele[N1[i]-1], node_ele[N2[i]-1], node_ele[N3[i]-1]])) for i in range(0, num_ele)])
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


# ************************* Partial Dam-Break Problem Test *************************
def PDB2poly(save_addr):
    # Writing .poly file
    # [1] create the output file
    poly_file = open(save_addr, mode="a")
    # [2] write summary for vertices
    poly_file.write("# Part of vertices\n")
    poly_file.write("8 2 0 1\n")
    # write vertices
    poly_file.write("1 0.0 0.0 1\n")
    poly_file.write("2 1.0 0.0 1\n")
    poly_file.write("3 1.0 0.8 1\n")
    poly_file.write("4 4.0 0.0 1\n")
    poly_file.write("5 4.0 2.0 1\n")
    poly_file.write("6 1.0 2.0 1\n")
    poly_file.write("7 1.0 1.2 1\n")
    poly_file.write("8 0.0 2.0 1\n")
    # [3] write summary for segments
    poly_file.write("# Part of segments\n")
    poly_file.write("9 1\n")
    poly_file.write("1 1 2 2\n")
    poly_file.write("2 2 4 1\n")
    poly_file.write("3 2 3 6\n")
    poly_file.write("4 4 5 1\n")
    poly_file.write("5 5 6 1\n")
    poly_file.write("6 6 7 6\n")
    poly_file.write("7 6 8 2\n")
    poly_file.write("8 8 1 2\n")
    poly_file.write("9 7 3 0\n")
    # [4] write summary for holes
    poly_file.write("# Part of holes\n")
    poly_file.write("0\n")

    poly_file.close()


def get_hIni_PDB(x, y):
    if x <= 1.0:
        return 0.6
    else:
        return 0.0


def Tri2IUHM_PDB(node_addr, num_node, edge_addr, num_edge, ele_addr, num_ele, voro_addr, save_addr):
    # open mesh file
    node_df = pd.read_table(node_addr, names=["ID", "x", "y", "Attribute"], sep="\s+", skiprows=lambda x: x in [0, num_node+1])
    edge_df = pd.read_table(edge_addr, names=["ID", "From", "To", "BD_marker"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    ele_df = pd.read_table(ele_addr, names=["ID", "N1", "N2", "N3"], sep="\s+", skiprows=lambda x: x in [0, num_ele+1])
    voro_df = pd.read_table(voro_addr, names=["ID", "E1", "E2", "x", "y"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    # Node information
    node_id = np.linspace(1, num_node, num_node)
    node_x = np.array(node_df["x"])
    node_y = np.array(node_df["y"])
    # ---get the elevation for the node
    node_ele = np.array([0.0 for i in range(0, num_node)])
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
    edge_BDaddr = np.array(["NoBCs".encode() for i in range(0, num_edge)])
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
    cell_n = np.array([0.0 for i in range(0, num_ele)])
    cell_rc = np.array([1.0 for i in range(0, num_ele)])
    h_ini = np.array([get_hIni_PDB(c_cx[i], c_cy[i]) for i in range(0, num_ele)])
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


# ************************* Dam Break on a Channel with Three Humps *************************
def ThreeHumps2poly(save_addr):
    # Writing .poly file
    # [1] create the output file
    poly_file = open(save_addr, mode="a")
    # [2] write summary for vertices
    poly_file.write("# Part of vertices\n")
    poly_file.write("6 2 0 1\n")
    # write vertices
    poly_file.write("1 0.0 0.0 1\n")
    poly_file.write("2 16.0 0.0 1\n")
    poly_file.write("3 75.0 0.0 1\n")
    poly_file.write("4 75.0 30.0 1\n")
    poly_file.write("5 16.0 30.0 1\n")
    poly_file.write("6 0.0 30.0 1\n")
    # [3] write summary for segments
    poly_file.write("# Part of segments\n")
    poly_file.write("7 1\n")
    poly_file.write("1 1 2 2\n")
    poly_file.write("2 2 3 2\n")
    poly_file.write("3 2 5 0\n")
    poly_file.write("4 3 4 2\n")
    poly_file.write("5 4 5 2\n")
    poly_file.write("6 5 6 2\n")
    poly_file.write("7 6 1 2\n")
    # [4] write summary for holes
    poly_file.write("# Part of holes\n")
    poly_file.write("0\n")

    poly_file.close()


def get_hIni_ThreeHumps(x, y):
    if x <= 16.0:
        return 1.875
    else:
        return 0.0


def get_DEM_ThreeHumps(x, y):
    return max(0.0, 1.0 - 0.125*np.sqrt((x-30.0)**2+(y-6.0)**2), 1.0 - 0.125*np.sqrt((x-30.0)**2+(y-24.0)**2),
               3.0 - 0.3*np.sqrt((x-47.5)**2+(y-15.0)**2))


def Tri2IUHM_ThreeHumps(node_addr, num_node, edge_addr, num_edge, ele_addr, num_ele, voro_addr, save_addr):
    # open mesh file
    node_df = pd.read_table(node_addr, names=["ID", "x", "y", "Attribute"], sep="\s+", skiprows=lambda x: x in [0, num_node+1])
    edge_df = pd.read_table(edge_addr, names=["ID", "From", "To", "BD_marker"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    ele_df = pd.read_table(ele_addr, names=["ID", "N1", "N2", "N3"], sep="\s+", skiprows=lambda x: x in [0, num_ele+1])
    voro_df = pd.read_table(voro_addr, names=["ID", "E1", "E2", "x", "y"], sep="\s+", skiprows=lambda x: x in [0, num_edge+1])
    # Node information
    node_id = np.linspace(1, num_node, num_node)
    node_x = np.array(node_df["x"])
    node_y = np.array(node_df["y"])
    # ---get the elevation for the node
    node_ele = np.array([get_DEM_ThreeHumps(node_x[i], node_y[i]) for i in range(0, num_node)])
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
    edge_BDaddr = np.array(["NoBCs".encode() for i in range(0, num_edge)])
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
    cell_n = np.array([0.018 for i in range(0, num_ele)])
    cell_rc = np.array([1.0 for i in range(0, num_ele)])
    h_ini = np.array([get_hIni_ThreeHumps(c_cx[i], c_cy[i]) for i in range(0, num_ele)])
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


# [1] Test: Stocker Dam Break
'''
poly_stocker = os.path.join("Case_Stocker/Case_stocker.poly")
Stocker2poly(poly_stocker)

Tri2IUHM_Stocker("Case_Stocker/Case_stocker.1.node", 8177, "Case_Stocker/Case_stocker.1.edge", 23973, "Case_Stocker/Case_stocker.1.ele", 15797,
         "Case_Stocker/Case_stocker.1.v.edge", "Case_Stocker/Case_Input")

prepareRainfall_zero(1, 70, 15797, "Case_Stocker/Case_Input")
'''

# [2] Test: Still Water over a 2D Bump
'''
poly_still = os.path.join("Case_Still/Case_still.poly")
#Still2poly(poly_still)

Tri2IUHM_Still("Case_Still/Case_still.1.node", 3291, "Case_Still/Case_still.1.edge", 9662, "Case_Still/Case_still.1.ele", 6372,
         "Case_Still/Case_still.1.v.edge", "Case_Still/Case_Input")

prepareRainfall_zero(1, 400, 6372, "Case_Still/Case_Input")
'''

# [3] Test: Partial Dam-Break Problem
'''
poly_pdb = os.path.join("Case_PDB/Case_PDB.poly")
PDB2poly(poly_pdb)
Tri2IUHM_PDB("Case_PDB/Case_PDB.1.node", 32020, "Case_PDB/Case_PDB.1.edge", 95433, "Case_PDB/Case_PDB.1.ele", 63414,
         "Case_PDB/Case_PDB.1.v.edge", "Case_PDB/Case_Input")
prepareRainfall_zero(1, 15, 63310, "Case_PDB/Case_Input")
'''

# [4] Test: Dam Break on a Channel with Three Humps
'''
poly_ThreeHumps = os.path.join("Case_ThreeHumps/Case_ThreeHumps.poly")
Tri2IUHM_ThreeHumps("Case_ThreeHumps/Case_ThreeHumps.1.node", 7311, "Case_ThreeHumps/Case_ThreeHumps.1.edge", 21613, "Case_ThreeHumps/Case_ThreeHumps.1.ele", 14303,
         "Case_ThreeHumps/Case_ThreeHumps.1.v.edge", "Case_ThreeHumps/Case_Input")
ThreeHumps2poly(poly_ThreeHumps)
prepareRainfall_zero(1, 400, 14303, "Case_ThreeHumps/Case_Input")
'''



