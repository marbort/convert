import numpy as np
import argparse
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from matplotlib import font_manager
import json

# Specify the path to your custom font file
custom_font_path = "/home/marco/.fonts/f/Formular.ttf"

# Register the custom font
prop = font_manager.FontProperties(fname=custom_font_path)

# Use the custom font in your plot
plt.rcParams["font.family"] = prop.get_name()


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    print(array[idx])
    return array[idx]


def extract_data(input, split, min, fac):
    data = {}
    datalist = {}
    freesplit = []
    with open(input, "r") as ifile:
        lines = ifile.readlines()
    data["full"] = {
        "cv": [float(x.split()[0]) * fac for x in lines if "#" not in x],
        "free": [float(x.split()[1]) for x in lines if "#" not in x],
    }
    with open("{}_equil".format(input), "r") as ifile:
        lines = ifile.readlines()
    data["equil"] = {
        "cv": [float(x.split()[0]) * fac for x in lines if "#" not in x],
        "free": [float(x.split()[1]) for x in lines if "#" not in x],
    }
    for i in range(split):
        with open("{}_{}".format(input, i), "r") as ifile:
            lines = ifile.readlines()
        data["split_{}".format(i)] = {"cv": [], "free": []}
        for line in lines:
            if "#" not in line:
                # if 'inf' not in line:
                data["split_{}".format(i)]["cv"].append(float(line.split()[0]) * fac)
                data["split_{}".format(i)]["free"].append(float(line.split()[1]))
        if i != 0:
            freesplit.append(data["split_{}".format(i)]["free"])
    data["mean_free"] = [np.mean(np.ma.masked_invalid(x)) for x in zip(*freesplit)]
    if min:
        cv_min = find_nearest(data["split_0"]["cv"], min)
        data["mean_free_scaled"] = [
            x - data["mean_free"][data["split_0"]["cv"].index(cv_min)]
            for x in data["mean_free"]
        ]
        data["full"]["free_scaled"] = [
            x - data["full"]["free"][data["split_0"]["cv"].index(cv_min)]
            for x in data["full"]["free"]
        ]
        data["equil"]["free_scaled"] = [
            x - data["equil"]["free"][data["split_0"]["cv"].index(cv_min)]
            for x in data["equil"]["free"]
        ]
    data["std_free"] = [np.std(np.ma.masked_invalid(x)) for x in zip(*freesplit)]
    data["std_free_mean"] = [
        (np.std(np.ma.masked_invalid(x)) / np.sqrt(split - 1)) for x in zip(*freesplit)
    ]
    # print(data['mean_free_scaled'][data['split_0']['cv'].index(cv_min)])
    # print(data["split_1"]['free'][-1],data["split_2"]['free'][-1],data["split_3"]['free'][-1],data["split_4"]['free'][-1],data['std_free'][-1],data['std_free_mean'][-1])
    print(data["std_free"][12], data["std_free_mean"][12])
    for entry in data:
        if isinstance(data[entry], dict):
            datalist[entry] = data[entry]
            for elm in data[entry]:
                if isinstance(data[entry][elm], np.ndarray):
                    datalist[entry][elm] = data[entry][elm].tolist()
        if isinstance(data[entry], np.ndarray):
            datalist[entry] = data[entry].tolist()

    with open(f"{input}_data", "w") as rawfile:
        summary = {"summary" : {}}
        for i in data:
            try:
                with open("temp", "w") as rawfile2:
                    json.dump({i:data[i]},rawfile2,indent=2)
                summary["summary"].update({i:data[i]})
            except:
                newdata=[]
                for j,item in enumerate(data[i]):
                    newdata.append(item)
                    if np.ma.is_masked(data[i][j]):
                        newdata[j]=np.inf
                    else:
                        newdata[j]=data[i][j]
                summary["summary"].update({i:newdata})
            """
            # print(i)
            newdata = []
            for j, item in enumerate(data[i]):
                newdata.append(item)
                print(data[i])
                if np.ma.is_masked(data[i][j]):
                    newdata[j] = np.inf
                else:
                    newdata[j] = data[i][j]
            summary["summary"].append({i: data[i]})
            """
        json.dump(summary, rawfile, indent=2)
        os.remove("temp")

                # for j,item in data[i]:
            #        if type(item) is  MaskedConstant:
            #    print(i,data[i])
            # newdata=data[i].tolist()
            #    rawfile.write(json.dumps[newdata])

    # with open(f"{input}_data",'w') as ofile:
    #   json.dump(datalist,ofile,indent=4)
    return data




def plot_data(
    data,
    min,
    full,
    equil,
    input,
    kcal,
    x,
    y,
    labels,
    inverted,
    single_color,
    pres,
    yrange,
    title,
    wd=10,
    hg=10,
):
    fig = plt.figure(figsize=(wd, hg), dpi=300)
    font = {"weight": "normal", "size": 30}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=4, marker="o", markersize=6)
    mpl.rcParams["axes.linewidth"] = 3

    """
    for j,i in enumerate(data):
        plt.errorbar([x*fac[j] for x in data[i]['cv']],data[i]['free'],yerr=data[i]['err'],ecolor=erclr,marker='o',markersize='3', label=i)
    plt.xlabel('CV')
    plt.ylabel('Free Energy ($kJ\ mol^{-1}$)')
    plt.legend()
    plt.savefig('umbrella_fes.png',format='png')
    """
    kj_to_kcal = 0.239006
    erclr = [
        "1f77b4",
        "ff7f0e",
        "2ca02c",
        "d62728",
        "9467bd",
        "8c564b",
        "e377c2",
        "7f7f7f",
        "bcbd22",
        "17becf",
    ]
    erclr_rgba = [
        [
            int("".join(x[0:2]), 16) / 255,
            int("".join(x[2:4]), 16) / 255,
            int("".join(x[4:6]), 16) / 255,
            0.2,
        ]
        for x in erclr
    ]
    # pres_colors=["cb5649","eba938","32c47f","C179B9","17d9ff","4cb944"]
    pres_colors = ["c1272d", "0000a7", "eba938", "008176", "b3b3b3", "4cb944"]
    pres_colors_rgba = [
        [
            int("".join(x[0:2]), 16) / 255,
            int("".join(x[2:4]), 16) / 255,
            int("".join(x[4:6]), 16) / 255,
            1,
        ]
        for x in pres_colors
    ]
    for j, i in enumerate(data):
        if pres:
            clr = pres_colors_rgba[j]
        elif single_color:
            clr = single_color
        else:
            clr = "C{}".format(j)
        if equil:
            if kcal:
                if min:
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["equil"]["cv"]],
                        [x * kj_to_kcal for x in data[i]["equil"]["free_scaled"]],
                        yerr=[
                            x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba
                    plt.fill_between(
                        data[i]["equil"]["cv"],
                        [
                            x * kj_to_kcal
                            - [
                                x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free_scaled"])
                        ],
                        [
                            x * kj_to_kcal
                            + [
                                x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free_scaled"])
                        ],
                        linewidth=0,
                        alpha=0.1,
                        color=clr,
                    )
                    if full:
                        markers, caps, bars = plt.plot(
                            data[i]["full"]["cv"],
                            [x * kj_to_kcal for x in data[i]["full"]["free_scaled"]],
                            "--",
                            color=clr,
                        )
                else:
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["equil"]["cv"]],
                        [x * kj_to_kcal for x in data[i]["equil"]["free"]],
                        yerr=[
                            x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba,
                    plt.fill_between(
                        data[i]["equil"]["cv"],
                        [
                            x * kj_to_kcal
                            - [
                                x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free"])
                        ],
                        [
                            x * kj_to_kcal
                            + [
                                x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free"])
                        ],
                        linewidth=0,
                        alpha=0.1,
                        color=clr,
                    )
                    if full:
                        markers, caps, bars = plt.plot(
                            data[i]["full"]["cv"],
                            [x * kj_to_kcal for x in data[i]["full"]["free"]],
                            "--",
                            color=clr,
                        )
            else:
                if min:
                    # print(data[i]['equil']['free_scaled'])
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["equil"]["cv"]],
                        data[i]["equil"]["free_scaled"],
                        yerr=[
                            x * 1.96 if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba
                    plt.fill_between(
                        data[i]["equil"]["cv"],
                        [
                            x
                            - [
                                x * 1.96 if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free_scaled"])
                        ],
                        [
                            x
                            + [
                                x * 1.96 if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free_scaled"])
                        ],
                        linewidth=0,
                        alpha=0.1,
                        color=clr,
                    )
                    if full:
                        markers, caps, bars = plt.plot(
                            data[i]["full"]["cv"],
                            data[i]["full"]["free_scaled"],
                            "--",
                            color=clr,
                        )

                else:
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["equil"]["cv"]],
                        data[i]["equil"]["free"],
                        yerr=[
                            x * 1.96 if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba,
                    plt.fill_between(
                        data[i]["equil"]["cv"],
                        [
                            x
                            - [
                                x * 1.96 if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free"])
                        ],
                        [
                            x
                            + [
                                x * 1.96 if x is not np.ma.masked else 0
                                for x in data[i]["std_free_mean"]
                            ][k]
                            for k, x in enumerate(data[i]["equil"]["free"])
                        ],
                        linewidth=0,
                        alpha=0.1,
                        color=clr,
                    )
                    if full:
                        markers, caps, bars = plt.plot(
                            data[i]["full"]["cv"],
                            data[i]["full"]["free"],
                            "--",
                            color=clr,
                        )
        else:
            if kcal:
                if min:
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["split_0"]["cv"]],
                        [x * kj_to_kcal for x in data[i]["mean_free_scaled"]],
                        yerr=[
                            x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba
                    if full:
                        plt.plot(
                            data[i]["full"]["cv"],
                            [x * kj_to_kcal for x in data[i]["full"]["free_scaled"]],
                            "--",
                            color=clr,
                        )
                else:
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["split_0"]["cv"]],
                        [x * kj_to_kcal for x in data[i]["mean_free"]],
                        yerr=[
                            x * 1.96 * kj_to_kcal if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba,
                    if full:
                        markers, caps, bars = plt.plot(
                            data[i]["full"]["cv"],
                            [x * kj_to_kcal for x in data[i]["full"]["free"]],
                            "--",
                            color=clr,
                        )
            else:
                if min:
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["split_0"]["cv"]],
                        data[i]["mean_free_scaled"],
                        yerr=[
                            x * 1.96 if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba
                    if full:
                        markers, caps, bars = plt.plot(
                            data[i]["full"]["cv"],
                            data[i]["full"]["free_scaled"],
                            "--",
                            color=clr,
                        )
                else:
                    markers, caps, bars = plt.errorbar(
                        [x for x in data[i]["split_0"]["cv"]],
                        data[i]["mean_free"],
                        yerr=[
                            x * 1.96 if x is not np.ma.masked else 0
                            for x in data[i]["std_free_mean"]
                        ],
                        ecolor=erclr_rgba[j],
                        label="_".join([i.split("/")[x] for x in range(-6, -1)]),
                        color=clr,
                    )  # ecolor=erclr_rgba,
                    if full:
                        markers, caps, bars = plt.plot(
                            data[i]["full"]["cv"],
                            data[i]["full"]["free"],
                            "--",
                            color=clr,
                        )
        [bar.set_alpha(0.0) for bar in bars]
    # plt.plot([x for x in data[i]["split_0"]['cv']],np.zeros(len(data[i]["split_0"]['cv'])))
    plt.axhline(y=0, color="gray", linestyle="-", linewidth=1)
    plt.xlabel(x)
    plt.ylabel(y)
    # plt.xlim([3.4,4])
    # if kcal:
    #    Y=[x*kj_to_kcal for x in data[i]['equil']['free_scaled']
    #    MINY=min()

    plt.ylim(yrange[0], yrange[1])

    if labels:
        # leg=plt.legend(labels,loc='upper center')
        leg = plt.legend(labels)
        for lh in leg.legend_handles:
            lh.set_alpha(1)
    # else:
    #    plt.legend()
    if kcal:
        plt.ylabel("Free Energy ($kcal\ mol^{-1}$)")
    if inverted:
        plt.gca().invert_xaxis()
    plt.title(title)
    plt.tight_layout()
    plt.savefig("umbrella_fes_{}.png".format(os.path.basename(input[0])), format="png")


def plot_data_bootstrap(
    data,
    min,
    full,
    equil,
    input,
    kcal,
    x,
    y,
    labels,
    inverted,
    single_color,
    pres,
    yrange,
    title,
    fac,
    wd=12,
    hg=10,
):
    fig,ax = plt.subplots(figsize=(wd, hg), dpi=300,layout='constrained')
    font = {"weight": "normal", "size": 30}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=4, marker="o", markersize=6)
    mpl.rcParams["axes.linewidth"] = 3

    """
    for j,i in enumerate(data):
        plt.errorbar([x*fac[j] for x in data[i]['cv']],data[i]['free'],yerr=data[i]['err'],ecolor=erclr,marker='o',markersize='3', label=i)
    plt.xlabel('CV')
    plt.ylabel('Free Energy ($kJ\ mol^{-1}$)')
    plt.legend()
    plt.savefig('umbrella_fes.png',format='png')
    """
    kj_to_kcal = 0.239006
    erclr = [
        "1f77b4",
        "ff7f0e",
        "2ca02c",
        "d62728",
        "9467bd",
        "8c564b",
        "e377c2",
        "7f7f7f",
        "bcbd22",
        "17becf",
    ]
    erclr_rgba = [
        [
            int("".join(x[0:2]), 16) / 255,
            int("".join(x[2:4]), 16) / 255,
            int("".join(x[4:6]), 16) / 255,
            0.2,
        ]
        for x in erclr
    ]
    # pres_colors=["cb5649","eba938","32c47f","C179B9","17d9ff","4cb944"]
    pres_colors = ["c1272d", "0000a7", "eba938", "008176", "b3b3b3", "4cb944"]
    pres_colors_rgba = [
        [
            int("".join(x[0:2]), 16) / 255,
            int("".join(x[2:4]), 16) / 255,
            int("".join(x[4:6]), 16) / 255,
            1,
        ]
        for x in pres_colors
    ]
   
    for j, i in enumerate(data):
        if len(fac) == 1:
            factor = fac[0]    
        else:
            factor=fac[j]
        
        cv=np.array(data[i][0])*factor
        if pres:
            clr = pres_colors_rgba[j]
        elif single_color:
            clr = single_color
        else:
            clr = "C{}".format(j)
        
        if equil:
            if kcal:
                kj_to_kcal = 0.239006
            else:
                kj_to_kcal = 1
            if min:
                if len(min) == 1:
                    cv_min_idx=np.argmin(np.abs(np.array(data[i][0])-min[0]/factor))
                else:
                    cv_min_idx=np.argmin(np.abs(np.array(data[i][0])-min[j]/factor))
                free_min=data[i][1][cv_min_idx]
            else:
                free_min=0
        
        
        free=(np.array(data[i][1])-free_min)*kj_to_kcal
        yerrplus=(np.add(data[i][1],data[i][2])-free_min)*kj_to_kcal
        yerrminus=(np.subtract(data[i][1],data[i][2])-free_min)*kj_to_kcal 
        
        markers, caps, bars = plt.errorbar(
            cv,
            free,
            yerr=data[i][2],
            ecolor=erclr_rgba[j],
            label="_".join([i.split("/")[x] for x in range(-6, -1)]),
            color=clr,
        )  # ecolor=erclr_rgba
        plt.fill_between(
            cv,
            yerrplus,
            yerrminus,
            linewidth=0,
            alpha=0.1,
            color=clr,
        )
        [bar.set_alpha(0.0) for bar in bars]
    # plt.plot([x for x in data[i]["split_0"]['cv']],np.zeros(len(data[i]["split_0"]['cv'])))
    plt.axhline(y=0, color="gray", linestyle="-", linewidth=1)
    plt.xlabel(x)
    plt.ylabel(y)
    # plt.xlim([3.4,4])
    # if kcal:
    #    Y=[x*kj_to_kcal for x in data[i]['equil']['free_scaled']
    #    MINY=min()

    plt.ylim(yrange[0], yrange[1])

    if labels:
        # leg=plt.legend(labels,loc='upper center')
        leg = plt.legend(labels,loc='lower right',fontsize=24)
        for lh in leg.legend_handles:
            lh.set_alpha(1)
    # else:
    #    plt.legend()
    if kcal:
        plt.ylabel("Free Energy ($kcal\ mol^{-1}$)")
    if inverted:
        plt.gca().invert_xaxis()
    plt.title(title)
    plt.tight_layout()
    plt.savefig("umbrella_fes_{}_bootstrap.png".format(os.path.basename(input[0])), format="png")



def plot_splits(data,labels=None):
    fig,axes = plt.subplots(figsize=(14, 16), dpi=150, nrows=len(data), ncols=1)
    font = {"family": "sans", "weight": "normal", "size": 32}
    mpl.rc("font", **font)
    mpl.rcParams["axes.linewidth"] = 3
    for k,i in enumerate(data):
        for j in range(5):
            try:
                axes[k].plot(
                    data[i][f"split_{j}"]["cv"],
                    data[i][f"split_{j}"]["free"],
                    label=f"split_{j}",
                )
                axes[k].title.set_text(labels[k])
            except:
                axes.plot(
                data[i][f"split_{j}"]["cv"],
                data[i][f"split_{j}"]["free"],
                label=f"split_{j}",
                )
                axes.title.set_text(labels[k])
            #plt.xlabel("CV")
            #plt.ylabel("Free Energy ($kJ\ mol^{-1}$)")
            #plt.title(i)
    #plt.legend(ncols=5)
    plt.tight_layout()
    plt.savefig("umbrella_fes_splits.png", format="png")


def main():
    parser = argparse.ArgumentParser(description="Plot data")
    parser.add_argument("--input", dest="input", type=str, help="input data", nargs="+")
    parser.add_argument(
        "--fac", dest="fac", type=float, help="cv conversion factors", nargs="+"
    )
    parser.add_argument(
        "--split",
        dest="split",
        type=int,
        help="number of split for each umbrella",
        nargs="+",
    )
    parser.add_argument(
        "--min", dest="min", type=float, help="CV value for 0", nargs="+"
    )
    parser.add_argument(
        "--full",
        dest="full",
        action="store_true",
        help="plot also FES using whole trajectory",
    )
    parser.add_argument(
        "--equil",
        dest="equil",
        action="store_true",
        help="plot also FES using all but first block",
    )
    parser.add_argument(
        "--kcal", dest="kcal", action="store_true", help="plot FES in kcal"
    )
    parser.add_argument(
        "--xlabel", dest="xlabel", default="CV", type=str, help="x label"
    )
    parser.add_argument(
        "--ylabel",
        dest="ylabel",
        default="Free Energy ($kJ\ mol^{-1}$)",
        type=str,
        help="y label",
    )
    parser.add_argument(
        "--yrange", dest="yrange", default=[-5, 35], nargs=2, type=float, help="y range"
    )
    parser.add_argument(
        "--labels",
        dest="labels",
        default="Fes",
        type=str,
        nargs="+",
        help="Plot labels label",
    )
    parser.add_argument(
        "--inverted", dest="inverted", action="store_true", help="Inverts X axis"
    )
    parser.add_argument(
        "--pres", dest="pres", action="store_true", help="use presentation colors"
    )
    parser.add_argument("--color", dest="color", type=str, help="choose single color")
    parser.add_argument("--title",  type=str, help="Plot title",default="FES Profile")
    args = parser.parse_args()
    data = {}
    data_bootstrap = {}
    for i, input in enumerate(args.input):
        try:
            if len(args.split) == len(args.input) and len(args.fac) == len(args.input):
                data[input] = extract_data(input, args.split[i], args.min[i], args.fac[i])
            else:
                data[input] = extract_data(input, args.split[0], args.min[i], args.fac[0])
            
        except:
            print("Single minimum defined")
            data[input] = extract_data(input, args.split[i], args.min, args.fac[i])
            
        data_bootstrap[input] = np.loadtxt(f"{input}_equil_errors",unpack=True)

    
    print(data[input].keys())
    plot_data(
        data,
        args.min,
        args.full,
        args.equil,
        args.input,
        args.kcal,
        args.xlabel,
        args.ylabel,
        args.labels,
        args.inverted,
        args.color,
        args.pres,
        args.yrange,
        args.title,
    )
    
    plot_data_bootstrap(
        data_bootstrap,
        args.min,
        args.full,
        args.equil,
        args.input,
        args.kcal,
        args.xlabel,
        args.ylabel,
        args.labels,
        args.inverted,
        args.color,
        args.pres,
        args.yrange,
        args.title,
        args.fac
    )
    
    plot_splits(data,args.labels)

    # print(len(data['free']),len(data['err']))


if __name__ == "__main__":
    main()
