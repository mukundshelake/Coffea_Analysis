import glob, os, ROOT
import ROOT
import os
import time
import paramiko
import secrets

problemFiles = []
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
			# print(setName)
			for file in glob.glob(os.path.join(setName,'*/*/*/*/*.root')):
				#print(file)
				if os.path.isfile(file):
					if isValidRootFile(file) and file not in problemFiles:
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


def getSamples(DirList):
	sampleFiles = {}
	for eraDir in DirList:
		for setName in glob.glob(os.path.join(eraDir,'*')):
			lst = []
			# print(setName)
			for file in glob.glob(os.path.join(setName,'*/*/*/*/tree_10.root')):
				#print(file)
				if os.path.isfile(file):
					if isValidRootFile(file):
						lst.append(file)
					else: 
						print(f"Excluding file {file} as invalid ROOT file")
			setName = setName.split('/')[-1]
			sampleFiles[setName] = lst
	return sampleFiles


def logScript(logFile, script_filename, coffeaFile, meta, runStatus, errorFile):
    # Get the name of the current script file
    with open(script_filename, 'r') as script_file:
        script_content = script_file.read()
    
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")

	# Header
    header = f"{timestamp}\t\t{script_filename}--->{coffeaFile}\n\n"
    status = ""
    if runStatus =="Success":
        status = f"Status: Successfull run; no error file was generated\n\n"
    elif runStatus == "Failure":
        status = f"Status: There was error; error log can be found at {errorFile}\n\n"
	

    # Create a formatted entry for the current script
    script_entry = f"Processor: {meta}\n" + "-"*80 + "-"*80 +f"\n{script_content}\n" + "-"*80 + "-"*80 +"\n"

    # Write the current script's entry with separator before the existing content
    with open(logFile, 'w') as output_file:
	    output_file.write(f"{header}{status}{script_entry}")

def scptoEOS(outputFile, eosLogDir):
    # Define the source file path (local file)
    filename = outputFile.split('/')[-1]
    # Define the destination file path on the remote server
    destination_file = os.path.join('/eos/user/m/mshelake/CoffeaIntro/Logs/', eosLogDir, filename)

    # SSH connection parameters
    hostname = secrets.hostname
    port = 22  # SSH port (typically 22)
    username = secrets.Username
    password = secrets.Password



    # Create an SSH client
    ssh_client = paramiko.SSHClient()

    try:
        # Automatically add the remote server's host key (this is insecure, use known_hosts in production)
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        # Connect to the remote server
        ssh_client.connect(hostname, port, username, password)

        # Create an SFTP client
        sftp = ssh_client.open_sftp()

        # Extract the directory path from the destination file
        destination_dir = "/".join(destination_file.split("/")[:-1])

        # Check if the destination directory exists, and create it if it doesn't
        try:
            sftp.stat(destination_dir)
        except FileNotFoundError:
            sftp.mkdir(destination_dir)

        # Copy the local file to the remote server
        sftp.put(outputFile, destination_file)

        # Close the SFTP client
        sftp.close()

        # Close the SSH connection
        ssh_client.close()

        print(f"File '{outputFile}' copied to '{hostname}:{destination_file}' successfully.")
    except Exception as e:
        print(f"Error: {e}")
		
def sayHello():
	print("Hello")