"""Cultured-Fib network/cached GO and library-level association replay."""
from pathlib import Path
import subprocess,sys
if __name__=='__main__':
    subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'reproduce.py'),'culture']+sys.argv[1:],check=True)
