import subprocess

# Run subprocess
def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, text=True)
    while True:
        line = process.stdout.readline()
        if line:
            print(f"Output: {line.strip()}")
        else: 
            break    
    return_code = process.wait()


