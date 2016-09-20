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
                                                   'Helicity Angles A', 'Helicity Angles B', 'Twist Angles A', 'Twist Angles B']).replace('nan', '0', regex=True).astype(str)
    for i in range(len(pdData['Helicity Angles A']) - 5, len(pdData['Helicity Angles A'])):
        # valCount = scstats.binned_statistic(eval(pdData['Current Deviation'][i]), eval(pdData['Current Deviation'][i]), statistic='count', bins=12, range=(-5, 5))
        print('Starting work on Frame {}'.format(i))
        fig1 = plt.figure(figsize=(10,10), tight_layout=True)
        ax1 = fig1.add_subplot(121)
        ax = fig1.add_subplot(122)
        ax.hist([j for j in eval(pdData['Helicity Angles A'][i]) if j != 0], bins=250, normed=True, label="Positive initial Deviation")
        ax1.hist([j for j in eval(pdData['Helicity Angles B'][i]) if j != 0], bins=250, normed=True, label="Double initial Deviation")
        ax.legend()
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
