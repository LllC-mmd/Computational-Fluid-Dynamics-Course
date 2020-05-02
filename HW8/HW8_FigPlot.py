import os
import re
import numpy as np
import matplotlib.pyplot as plt


case_dict = {100: 0, 200: 1, 500: 2, 1000: 3}
var_dict = {"Rho": 0, "U": 1, "P": 2}
annoate_dict = {"Minimod": ["forestgreen", "^", 3], "SuperBee": ["navy", "D", 2]}

fig, ax = plt.subplots(4, 3, figsize=(12, 12))

res_Table = {}

for root, dirs, files in os.walk("HW8_res"):
    if len(re.findall(r"\/(\d+)", root))>0:
        num_grid = int(re.findall(r"\/(\d+)", root)[0])
        for f in files:
            if len(re.findall(r".txt", f))>0:
                resArray = np.loadtxt(os.path.join(root, f))
                l_res = len(resArray)
                var = re.findall(r"_([a-zA-Z]+)_", f)[0]
                schemeName = re.findall(r"([a-zA-Z]+)_", f)[0]
                caseNo = case_dict[num_grid]
                x_loc = [-0.5 + (i + 0.5) * 1.0/num_grid for i in range(0, l_res)]
                res_Table[schemeName + str(num_grid)] = resArray
                if schemeName == "Exact":
                    ax[int(caseNo)][var_dict[var]].plot(x_loc, resArray, color="red", label=schemeName, zorder=1)

                else:
                    ax[int(caseNo)][var_dict[var]].scatter(x_loc, resArray, s=2, c=annoate_dict[schemeName][0],
                                                        marker=annoate_dict[schemeName][1], label=schemeName, zorder=annoate_dict[schemeName][2])

# Plot the results
col_label = ["$\\rho$", "U", "P"]
row_label = ["Case{}".format(i) for i in [1, 2, 3, 4]]

for a, col in zip(ax[0], col_label):
    a.set_title(col, size='xx-large')

for a, row in zip(ax[:,0], row_label):
    a.set_ylabel(row, size='large', rotation=90)

for i in range(0, 4):
    for j in range(0, 3):
        ax[i, j].set_xticks(np.linspace(-0.5, 0.5, 6))
        ax[i, j].legend(fontsize='medium', loc=1)

#plt.show()
plt.savefig("HW8_TVD.png", dpi=400)

# Calculate errors
for nx in [100, 200, 500, 1000]:
    print("L_inf-error of num_grid="+str(nx))
    print("SuperBee limiter: ", np.max(np.abs(res_Table["SuperBee"+str(nx)]-res_Table["Exact"+str(nx)])), "\t"
            "Minimod limiter: ", np.max(np.abs(res_Table["Minimod"+str(nx)]-res_Table["Exact"+str(nx)])))
    print("L_1-error of num_grid="+str(nx))
    print("SuperBee limiter: ", np.mean(np.abs(res_Table["SuperBee"+str(nx)]-res_Table["Exact"+str(nx)])), "\t"
            "Minimod limiter: ", np.mean(np.abs(res_Table["Minimod"+str(nx)]-res_Table["Exact"+str(nx)])))