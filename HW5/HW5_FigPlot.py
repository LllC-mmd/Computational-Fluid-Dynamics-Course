import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def shift(x, start, end, T):
    x1 = math.floor((x - start) / T)
    x2 = math.floor((x - end) / T)
    if x1 < x2:
        x_shift = x - x2 * T
    else:
        x_shift = x - x1 * T
    return x_shift


def case_ini(x, y, t=2):
    gamma = 1.4
    vs = 5.0
    x0 = 5.0
    y0 = 5.0
    x_shift = shift(x, 0, 10, 10)
    y_shift = shift(y, 0, 10, 10)
    r_square = (x_shift-(x0+1*t))**2 + (y_shift-(y0+1*t))**2
    res = np.power(1-(gamma-1)*vs*vs/(8*gamma*np.pi*np.pi)*np.exp(1-r_square), 1/(gamma-1))
    return res


paraDict = {"0": [0.02, 0.2], "1": [0.01, 0.1], "2": [0.005, 0.05], "3": [0.0025, 0.025],
            "4": [0.00125, 0.0125], "5": [0.000625, 0.00625]}

fig, ax = plt.subplots(5, 4, figsize=(15, 16))

for root, dirs, files in os.walk("HW5_res"):
    if (len(files) > 0) and (len(re.findall(r"/para(\d+)_", root))>0):
        for f in files:
            paraNo = int(re.findall(r"/para(\d+)_", root)[0])
            schemeName = re.findall(r"_(\w+)_", f)[0]
            timeStamp = float(re.findall(r"(\d+.\d+).txt", f)[0])
            # calculate the Mean Absolute Error
            if timeStamp == 2.0:
                resArray = np.loadtxt(os.path.join(root, f))
                l_res = len(resArray)
                x_loc = [0.0 + (i + 0.5) * paraDict[str(paraNo)][1] for i in range(0, l_res)]
                y_loc = [0.0 + (j + 0.5) * paraDict[str(paraNo)][1] for j in range(0, l_res)]
                rho_ref = np.reshape([case_ini(xi, yi) for xi in x_loc for yi in y_loc], (l_res, l_res))
                res_diff = resArray-rho_ref
                mae = np.mean(np.abs(res_diff))
                ini_error = np.max(np.abs(res_diff))
                # print(schemeName, "\t", paraNo, "\tLength: ", l_res, "\tMAE: ", mae, "\tL_ini: ", ini_error)
                para_x = str(paraDict[str(paraNo)][1])
                para_t = str(paraDict[str(paraNo)][0])
                if schemeName == "2rdJameson":
                    fmt = ticker.ScalarFormatter(useMathText=True)
                    fmt.set_powerlimits((0, 0))
                    ax00 = ax[paraNo, 0].imshow(X=resArray, origin='lower')
                    fig.colorbar(ax00, ax=ax[paraNo, 0], format=fmt)
                    ax[paraNo, 0].set_title("$\Delta x$=$\Delta y$="+para_x+", $\Delta t$="+para_t, fontsize='medium')
                    ax[paraNo, 0].set_xticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 0].set_xticklabels(np.linspace(0, 10, 6))
                    ax[paraNo, 0].set_yticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 0].set_yticklabels(np.linspace(0, 10, 6))
                    ax01 = ax[paraNo, 1].imshow(X=res_diff, origin='lower', cmap="coolwarm")
                    ax[paraNo, 1].set_xticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 1].set_xticklabels(np.linspace(0, 10, 6))
                    ax[paraNo, 1].set_yticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 1].set_yticklabels(np.linspace(0, 10, 6))
                    fig.colorbar(ax01, ax=ax[paraNo, 1], format=fmt)
                else:
                    fmt = ticker.ScalarFormatter(useMathText=True)
                    fmt.set_powerlimits((0, 0))
                    ax10 = ax[paraNo, 2].imshow(X=resArray, origin='lower')
                    fig.colorbar(ax10, ax=ax[paraNo, 2], format=fmt)
                    ax[paraNo, 2].set_title("$\Delta x$=$\Delta y$=" + para_x + ", $\Delta t$=" + para_t, fontsize='medium')
                    ax[paraNo, 2].set_xticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 2].set_xticklabels(np.linspace(0, 10, 6))
                    ax[paraNo, 2].set_yticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 2].set_yticklabels(np.linspace(0, 10, 6))
                    ax11 = ax[paraNo, 3].imshow(X=res_diff, origin='lower', cmap="coolwarm")
                    fig.colorbar(ax11, ax=ax[paraNo, 3], format=fmt)
                    ax[paraNo, 3].set_xticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 3].set_xticklabels(np.linspace(0, 10, 6))
                    ax[paraNo, 3].set_yticks(np.linspace(0, l_res, 6))
                    ax[paraNo, 3].set_yticklabels(np.linspace(0, 10, 6))

# plt.show()
plt.savefig("Res.png", dpi=400)
