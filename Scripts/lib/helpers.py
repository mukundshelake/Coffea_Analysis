import glob, os
def getfileset(eraDir):
	fileset = {}
	for setName in glob.glob(os.path.join(eraDir,'*/*')):
		lst = []
		for file in glob.glob(os.path.join(setName,'*/*/*/*/*.root')):
			print(file)
			if os.path.isfile(file):
				lst.append(file)
		setName = setName.split('/')[-1]
		fileset[setName] = lst
	return fileset


def getFiles(folder_path):
    file_list = []
    for file_name in os.listdir(folder_path):
        if os.path.isfile(os.path.join(folder_path, file_name)):
            file_list.append(os.path.join(folder_path, file_name))
    return file_list

