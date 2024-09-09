import subprocess

# Run subprocess
def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    while True:
        line = process.stdout.readline()
        if not line: break    
    return_code = process.wait()


