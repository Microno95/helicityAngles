import os, numpy, multiprocessing, matplotlib.pyplot as plt, sys, pandas as pd, math, shutil, scipy.stats as scstats

# Set some Pandas options
pd.set_option('html', False)
pd.set_option('max_columns', 30)
pd.set_option('max_rows', 20)

def ensureDir(dirPath):
    try:
        os.makedirs(dirPath)
    except OSError:
        if not os.path.isdir(dirPath):
            raise
            
def analyseSaveData(dataPath, savePath="./"):
    pdData = pd.read_csv(dataPath, sep=";", names=['Time', 'Current Deviation A', 'Current Deviation B', 'Derivative Deviation A', 'Derivative Deviation B',
                                                   'Helicity Angles', 'Twist Angles']).replace('nan', '0', regex=True).astype(str)
    # valCount = scstats.binned_statistic(eval(pdData['Current Deviation'][i]), eval(pdData['Current Deviation'][i]), statistic='count', bins=12, range=(-5, 5))
    i = 0
    print('Starting work on Frame {}'.format(i))
    relevantDataHelicity = [eval(i)[0] for i in pdData['Helicity Angles']]
    relevantDataTwist = [eval(i)[0] for i in pdData['Twist Angles']]
    time = [float(i) for i in pdData['Time']]
    fig1 = plt.figure(figsize=(20,20), tight_layout=True)
    ax1 = fig1.add_subplot(221)
    ax2 = fig1.add_subplot(223)
    ax3 = fig1.add_subplot(222)
    ax4 = fig1.add_subplot(224)
    n, bins, patches = ax1.hist(relevantDataHelicity, bins=200, histtype='step', normed=True, label="Helicity Spectra")
    mu = numpy.mean(relevantDataHelicity)
    sigma = numpy.std(relevantDataHelicity)
    y = scstats.norm.pdf(bins, mu, sigma)
    l = ax1.plot(bins, y, 'k--', linewidth=1.5, label="PDF of Helicity Spectra")
    n, bins, patches = ax2.hist(relevantDataTwist, bins=200, histtype='step', normed=True, label="Twist Spectra")
    mu = numpy.mean(relevantDataTwist)
    sigma = numpy.std(relevantDataTwist)
    y = scstats.norm.pdf(bins, mu, sigma)
    l = ax2.plot(bins, y, 'k--', linewidth=1.5, label="PDF of Twist Spectra")
    ax3.plot(time, relevantDataHelicity, label="Helicity Angle")
    ax4.plot(time, relevantDataTwist, label="Twist Angle")
    ax3.set_xscale('log')
    ax4.set_xscale('log')
    ax3.set_yscale('log')
    ax4.set_yscale('log')
    ax1.set_xlabel("Helicity Angle")
    ax2.set_xlabel("Twist Angle")
    ax3.set_xlabel("Time")
    ax4.set_xlabel("Time")
    ax3.set_ylabel("Helicity Angle")
    ax4.set_ylabel("Twist Angle")
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.tight_layout()
    fig1.savefig(savePath + 'Frame_{}.png'.format(i), dpi=150)
    plt.close(fig1)
    del fig1
    
    
if __name__ == "__main__":
    dataPath = "/home/perse/Documents/Helicity_Work/"
    for file in sorted(os.listdir(dataPath + "Variational Lyapunov/")):
        if file.lower().endswith(".txt"):
            print("Starting work on {}\n".format(file))
            savePath = dataPath + "/Analysis/" + file.lower().replace(".txt", "").replace(" ", "_").replace("-", "")  + "/"
            ensureDir(savePath)
            analyseSaveData(dataPath + "Variational Lyapunov/" + file, savePath=savePath)
            print("Completed work on {}".format(file))
            print(shutil.get_terminal_size((80, 20))[0]*"-")
        if '--single' in sys.argv:
            break
