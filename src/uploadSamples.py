import subprocess
import json
# List of files to be transferred
import paramiko
for era in ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']:
    with open(f'Datasets/sampleFiles_{era}.json', 'r') as json_file:
        dicti = json.load(json_file)
        for dataMC in dicti:
            for pr in dicti[dataMC]:
                datasetName = f'{era}_{pr}'
                for file in dicti[dataMC][pr]:
                    # print(f"Transferring {file}")
                    destination = f"/eos/user/m/mshelake/myAnanlysis/datasetSamples/{datasetName}.root"
                    try:
                        # Using pexpect to run the scp command and handle password prompt
                        scp_command = f"scp {file} {destination}"
                        print(scp_command)
                        hostname = "lxplus.cern.ch"
                        port = 22  # SSH port (typically 22)
                        username = "mshelake"
                        password = "Mkshp400"
                        # Create an SSH client
                        ssh_client = paramiko.SSHClient()
                        # Automatically add the remote server's host key (this is insecure, use known_hosts in production)
                        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

                        # Connect to the remote server
                        ssh_client.connect(hostname, port, username, password)

                        # Create an SFTP client
                        sftp = ssh_client.open_sftp()

                        # Copy the local file to the remote server
                        sftp.put(file, destination)

                        # Close the SFTP client
                        sftp.close()

                        # Close the SSH connection
                        ssh_client.close()

                        print(f"File '{file}' copied to '{hostname}:{destination}' successfully.")
                    except Exception as e:
                        print(f"Error: {e}")
