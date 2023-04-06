import sys,re,os
import numpy as np
sys.path.append('/home/mathias/ftsrader')
from ftsreader import ftsreader

class cell_spectra:
    def __init__(self,csdir,nocsdir=''):
        self.csdir = csdir # directory containing cell spectra
        self.nocsdir=nocsdir # directory containg spectra without cell
        self.spectra=[]

        self.wvnrange={'HBR': [2400, 2540],
                       'N2O': [2150, 2250]}
    
    def load_cellspectra(self, snippet, date='last'):
        m = re.compile(snippet,re.I)
        files = filter(m.search,os.listdir(self.csdir))
        for ff in files:
            try:
                self.spectra.append(ftsreader.ftsreader(os.path.join(self.csdir,ff), verbose=False, getspc=True, getifg=False))
            except Exception as e:
                print(e)

            continue
            
    def average_cellspectra(self, direction='both',type='HBR'):
        nr = 0
        for s in self.spectra:
            if nr == 0:
                self.wvn = s.spcwvn[(s.spcwvn < self.wvnrange[type][1]) & (s.spcwvn > self.wvnrange[type][0])]
                self.average = s.spc[(s.spcwvn < self.wvnrange[type][1]) & (s.spcwvn > self.wvnrange[type][0])]
                nr = nr + 1
                continue
            self.average = self.average + s.spc[(s.spcwvn < self.wvnrange[type][1]) & (s.spcwvn > self.wvnrange[type][0])]
            nr = nr+1
        self.average = self.average/nr

    def norm_cellspectra(self,snippet,plot=False,norm=True):
        m = re.compile(snippet,re.I)
        files = filter(m.search,os.listdir(self.nocsdir))
        for ff in files:
            spc = ftsreader.ftsreader(os.path.join(self.nocsdir,ff), verbose=False, getspc=True, getifg=False)
        
        espec = np.interp(self.wvn,spc.spcwvn,spc.spc)

        self.average = self.average/espec

        if norm:
            self.average = self.average/np.max(self.average)

        if plot:
            import matplotlib.pyplot as plt
            f = plt.figure()
            plt.plot(self.wvn,self.average)
            f.show()
            input()
        
if __name__ == '__main__':

    cs = cell_spectra(sys.argv[1],sys.argv[2])
    cs.load_cellspectra(snippet='lhi04b')
    cs.average_cellspectra(type='HBR')
    cs.norm_cellspectra(snippet='l0i04b',plot=True,norm=True)

    np.savetxt(sys.argv[3], np.vstack((cs.wvn,cs.average)).T)
