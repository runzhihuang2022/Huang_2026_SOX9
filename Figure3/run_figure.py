from pathlib import Path
import subprocess,sys
if __name__=='__main__':
    subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/reproduce.py'),'spatial']+sys.argv[1:],check=True)
