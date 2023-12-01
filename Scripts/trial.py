from lib.helpers import getSamples, scptoEOS
import paramiko, os
dirAddress = ['/nfs/home/common/RUN2_UL/Tree_crab/SEVENTEEN/MC']
sampleList = getSamples(dirAddress)
print(sampleList)

# SSH connection parameters
hostname = "lxplus.cern.ch"
port = 22  # SSH port (typically 22)
username = "mshelake"
password = "Mkshp400"
# Create an SSH client
ssh_client = paramiko.SSHClient()

for dset in sampleList:
    if len(sampleList[dset]) < 1:
        continue
    file = sampleList[dset][0]
    filename = f"2017_{dset}.root"
    destination_file = os.path.join('/eos/user/m/mshelake/CoffeaIntro/Files/', filename)
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
        sftp.put(file, destination_file)

        # Close the SFTP client
        sftp.close()

        # Close the SSH connection
        ssh_client.close()

        print(f"File '{file}' copied to '{hostname}:{destination_file}' successfully.")
    except Exception as e:
        print(f"Error: {e}")