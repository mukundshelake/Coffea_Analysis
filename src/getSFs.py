import subprocess
from SFs.sfPaths import SFs

for sf in SFs:
    publicHash = SFs[sf]['publicHash']
    srcPath = SFs[sf]['srcPath']
    fileName = SFs[sf]['fileName']
    targetPath = SFs[sf]['targetPath']

    command = f"curl -JL https://cernbox.cern.ch/remote.php/dav/public-files/{publicHash}/{srcPath}/{fileName} -o {targetPath}/{fileName}"
    print(command)
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Read the output and error (if any)
    output, error = process.communicate()

    # Print the output
    print("Output:")
    print(output.decode())

    # Print the error (if any)
    if error:
        print("Error:")
        print(error.decode())