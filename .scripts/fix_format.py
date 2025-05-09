import subprocess

cmd = ("ruff", "format", "../pycap/")
proc = subprocess.run(cmd)
assert proc.returncode == 0, f"Error running command: {' '.join(cmd)}"

cmd = ("ruff", "check", "../pycap/", "--fix")
proc = subprocess.run(cmd)
assert proc.returncode == 0, f"Error running command: {' '.join(cmd)}"
