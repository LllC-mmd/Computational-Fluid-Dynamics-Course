import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt


def shift(x, start, end, T):
    x1 = math.floor((x - start) / T)
    x2 = math.floor((x - end) / T)
    if x1 < x2:
        x_shift = x - x2 * T
    else:
        x_shift = x - x1 * T
    return x_shift


def case1_ini(x):
    x_shift = shift(x, -0.5, 0.5, 1)
    if (x_shift >= -0.5) and (x_shift < -0.25):
        return 0
    elif (x_shift >= -0.25) and (x_shift < 0.25):
        return 1
    else:
        return 0


def case2_ini(x):
    x_shift = shift(x, -0.5, 0.5, 1)
    return np.sin(4*np.pi*x_shift)


paraDict = {"1": "(0.01, 0.005)", "2": "(0.001, 0.0005)", "3": "(0.0001, 0.00005)", "4": "(0.00002, 0.00001)"}
locDict = {1: 0, 5: 1, 10: 2}
x_id = np.linspace(-0.5, 0.5, 100001)
t1_ref1 = np.array([case1_ini(i - 1) for i in x_id])
t5_ref1 = np.array([case1_ini(i - 5) for i in x_id])
t10_ref1 = np.array([case1_ini(i - 10) for i in x_id])
t1_ref2 = np.array([case2_ini(i - 1) for i in x_id])
t5_ref2 = np.array([case2_ini(i - 5) for i in x_id])
t10_ref2 = np.array([case2_ini(i - 10) for i in x_id])

fig, ax = plt.subplots(4, 3, sharex=True, figsize=(20, 20))

for i in range(0, 4):
    ax[i, 0].plot(x_id, t1_ref1, color="red", label="Exact")
    ax[i, 1].plot(x_id, t5_ref1, color="red", label="Exact")
    ax[i, 2].plot(x_id, t10_ref1, color="red", label="Exact")
    '''
    ax[i, 0].plot(x_id, t1_ref2, color="red", label="Exact")
    ax[i, 1].plot(x_id, t5_ref2, color="red", label="Exact")
    ax[i, 2].plot(x_id, t10_ref2, color="red", label="Exact")
    '''

for root, dirs, files in os.walk("HW4_res"):
    if (len(files) > 0) and (len(re.findall(r"/para(\d+)_", root))>0):
        for f in files:
            paraNo = int(re.findall(r"/para(\d+)_", root)[0])
            caseNo = int(re.findall(r"Case(\d+)", f)[0])
            schemeName = re.findall(r"_(\w+)_", f)[0]
            timeStamp = int(re.findall(r"(\d+).txt", f)[0])
            # calculate the Mean Absolute Error
            if (timeStamp in [1, 5, 10]) and (caseNo == 1):
            # if (timeStamp in [1, 5, 10]) and (caseNo == 2):
                resArray = np.loadtxt(os.path.join(root, f))
                l_res = len(resArray)
                x_loc = np.linspace(-0.5, 0.5, l_res)
                # trueArray = np.array([case1_ini(x - timeStamp) for x in x_loc])
                trueArray = np.array([case2_ini(x - timeStamp) for x in x_loc])
                mae = np.mean(np.abs(resArray-trueArray))
                if schemeName == "2rdUpWind":
                    ax[paraNo-1, locDict[timeStamp]].scatter(x_loc, resArray, s=0.5, c="darkorange", marker="x", label="2rdUpWind: %.2f %%"% (mae*100))
                else:
                    ax[paraNo-1, locDict[timeStamp]].scatter(x_loc, resArray, s=0.5, c="navy", marker="x", label="2rdCenter: %.2f %%"% (mae*100))

col_label = ["t={}".format(i) for i in [1, 5, 10]]
row_label = ["$\Delta x=0.01$\n$\Delta t=0.005$", "$\Delta x=0.001$\n$\Delta t=0.0005$",
             "$\Delta x=0.0001$\n$\Delta t=0.00005$", "$\Delta x=0.00002$\n$\Delta t=0.00001$"]

for a, col in zip(ax[0], col_label):
    a.set_title(col, size='xx-large')

for a, row in zip(ax[:,0], row_label):
    a.set_ylabel(row, rotation=90)

for i in range(0, 4):
    for j in range(0, 3):
        ax[i, j].set_xticks(np.linspace(-0.5, 0.5, 6))
        # Set for Case1
        ax[i, j].set_yticks(np.linspace(-0.5, 1.5, 11))
        # ax[i, j].set_yticks(np.linspace(-1.0, 1.0, 11))
        ax[i, j].legend(fontsize='medium', loc=1)

# plt.show()
plt.savefig("Case1_Res.png", dpi=400)
# plt.savefig("Case2_Res.png", dpi=400)
                

