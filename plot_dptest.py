import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy


def plot_energy(input):
    file = []
    file = np.loadtxt(input, unpack=True)
    fig = plt.figure(figsize=(8, 8), dpi=150)
    plt.bar(
        [i for i in range(len(file[0]))],
        [x - file[1][i] for i, x in enumerate(file[0])],
    )
    plt.ylim([0, 1])
    plt.title("_".join(input.split("#")[9:]))
    plt.tight_layout()
    plt.savefig(input.replace(".out", "_hist.png"), format="png")
    plt.close()

    fig = plt.figure(figsize=(8, 8), dpi=150)
    plt.plot(file[0] - min(file[0]), file[1] - min(file[1]), "o", markersize=3)
    plt.plot(file[0] - min(file[0]), file[0] - min(file[0]), linewidth=2)
    plt.tight_layout()
    plt.savefig(input.replace(".out", ".png"), format="png")
    plt.close()


def plot_force(input):
    file = np.loadtxt(input, unpack=True)
    fig = plt.figure(figsize=(8, 8), dpi=150)
    corr = []
    for i in range(3):
        plt.subplot(3, 1, i + 1)
        plt.plot(file[i], file[3 + i], "o", markersize=3)
        plt.plot(file[3 + i], file[3 + i])
        corr.append(scipy.stats.pearsonr(file[i], file[3 + i]))
    fig.suptitle("_".join(input.split("#")[9:]))
    print(input, corr)
    for i in corr:
        if i[0] < 0.5:
            with open("bad_file.log", "a") as ofile:
                ofile.write(
                    "{}\n".format(input[1:].replace("#", "/").replace(".f.out", ""))
                )
                good = False
            break
        else:
            good = True
        if good:
            with open("good_files.log", "a") as ofile:
                ofile.write(
                    "{}\n".format(input[1:].replace("#", "/").replace(".f.out", ""))
                )
    plt.savefig(input.replace(".out", ".png"), format="png")
    plt.close()


files_e = glob.glob("*e.out")
files_f = glob.glob("*f.out")

for file in files_e:
    try:
        plot_energy(file)
    except:
        print(file)
        # exit()
for file in files_f:
    plot_force(file)