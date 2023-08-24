import glob, os
def getfileset(eraDir):
	fileset = {}
	for setName in glob.glob(os.path.join(eraDir,'*/*')):
		if 'mu' in setName:
			continue
		lst = []
		for file in glob.glob(os.path.join(setName,'*/*/*/*/*.root')):
			print(file)
			if os.path.isfile(file):
				if 'mu' in file:
					pass
				lst.append(file)
		setName = setName.split('/')[-1]
		fileset[setName] = lst
	return fileset
rr = '/nfs/home/common/RUN2_UL/Tree_crab/SIXTEEN_preVFP'
ff = getfileset(rr)
print(ff.keys())		
	

