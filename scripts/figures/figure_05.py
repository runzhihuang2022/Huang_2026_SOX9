"""Same locked eight-library culture analysis as Figure4, in a new output directory."""
from pathlib import Path
import subprocess,sys
if __name__=='__main__':
    subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'reproduce.py'),'culture']+sys.argv[1:],check=True)
