import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt


reportNo = "7"
reportGrid = 300
var_dict = {"Rho": 0, "U": 1, "P": 2}
annoate_dict = {"Godunov": ["purple", "s", 2], "HLL": ["forestgreen", "^", 2], "HLLC": ["navy", "D", 3]}

fig, ax = plt.subplots(5, 3, figsize=(12, 15))

for root, dirs, files in os.walk("HW6_res"):
    if len(re.findall(r"(\d+)_HW", root))>0:
        hw_No = re.findall(r"\d+_HW(\d+)", root)[0]
        num_grid = int(re.findall(r"(\d+)_HW", root)[0])
        if (hw_No == reportNo) and (num_grid==reportGrid):
            for f in files:
                if len(re.findall(r"Case", f))>0:
                    resArray = np.loadtxt(os.path.join(root, f))
                    l_res = len(resArray)
                    var = re.findall(r"_([a-zA-Z]+)_", f)[0]
                    schemeName = re.findall(r"Case\d+([a-zA-Z]+)_", f)[0]
                    caseNo = int(re.findall(r"Case(\d+)[a-zA-Z]", f)[0])-1
                    x_loc = [0.0 + (i + 0.5) * 1.0/num_grid for i in range(0, l_res)]
                    if schemeName == "Exact":
                        ax[int(caseNo)][var_dict[var]].plot(x_loc, resArray, color="red", label=schemeName, zorder=1)
                    else:
                        ax[int(caseNo)][var_dict[var]].scatter(x_loc, resArray, s=2, c=annoate_dict[schemeName][0],
                                                               marker=annoate_dict[schemeName][1], label=schemeName, zorder=annoate_dict[schemeName][2])


col_label = ["$\\rho$", "U", "P"]
row_label = ["Case{}".format(i) for i in [1, 2, 3, 4, 5]]

for a, col in zip(ax[0], col_label):
    a.set_title(col, size='xx-large')

for a, row in zip(ax[:,0], row_label):
    a.set_ylabel(row, size='large', rotation=90)

for i in range(0, 5):
    for j in range(0, 3):
        ax[i, j].set_xticks(np.linspace(0, 1.0, 6))
        ax[i, j].legend(fontsize='medium', loc=1)

#plt.show()
plt.savefig("HW7_300.png", dpi=400)