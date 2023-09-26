import glob, os, ROOT
import ROOT
import os

def isValidRootFile(fname):
    if not os.path.exists(os.path.expandvars(fname)):
        return False
    try:
        f = ROOT.TFile(fname)
    except Exception as e:
        return False

    if not f:
        return False

    try:
        return not (f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered) or f.GetListOfKeys().IsEmpty())
    finally:
        f.Close()




def getfileset(DirList):
	fileset = {}
	for eraDir in DirList:
		for setName in glob.glob(os.path.join(eraDir,'*')):
			lst = []
			print(setName)
			for file in glob.glob(os.path.join(setName,'*/*/*/*/*.root')):
				#print(file)
				if os.path.isfile(file):
					if isValidRootFile(file):
						lst.append(file)
					else: 
						print(f"Excluding file {file} as invalid ROOT file")
			setName = setName.split('/')[-1]
			fileset[setName] = lst
	return fileset



def getFiles(folder_path):
    file_list = []
    for file_name in os.listdir(folder_path):
        if os.path.isfile(os.path.join(folder_path, file_name)):
            file_list.append(os.path.join(folder_path, file_name))
    return file_list

