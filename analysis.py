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
    pdData = pd.read_csv(dataPath, sep=";", names=['Time', 'Current Deviation', 'Derivative Deviation',
                                                   'Helicity Angles', 'Twist Angles']).replace('nan', '0', regex=True).astype(str)
    for i in range(len(pdData['Current Deviation'])):
        valCount = scstats.binned_statistic(eval(pdData['Current Deviation'][i]), eval(pdData['Current Deviation'][i]), statistic='count', bins=12, range=(-5, 5))
        print('Starting work on Frame {}'.format(i))
        fig1 = plt.figure(figsize=(10,10))
        ax = fig1.add_subplot(111)
        ax.plot(valCount.bin_edges[0:12], valCount.statistic)
        fig1.savefig(savePath + 'Frame_{}.png'.format(i), dpi=600)
        plt.close(fig1)
        del fig1
    
    
if __name__ == "__main__":
    dataPath = "/home/perse/Documents/Helicity_Work/"
    for file in os.listdir(dataPath + "Variational Lyapunov/"):
        if file.lower().endswith(".txt"):
            print("Starting work on {}\n".format(file))
            savePath = dataPath + "/Analysis/" + file.lower().replace(".txt", "").replace(" ", "_").replace("-", "")  + "/"
            ensureDir(savePath)
            analyseSaveData(dataPath + "Variational Lyapunov/" + file, savePath=savePath)
            print("Completed work on {}".format(file))
            print(shutil.get_terminal_size((80, 20))[0]*"-")
        if '--single' in sys.argv:
            break
