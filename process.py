import os
if not os.path.exists('db'):
    exit(0)
os.system('rm -rf result/*')
filepaths = []
for path, dirs, files in os.walk('db'):
    for file in files:
        if file.endswith('.fasta'):
            filepath = os.path.join(path, file)
            filepath = filepath.replace('\\', '/')
            filepaths.append(filepath)


cmd_sequence = [100, 0, 1]
for t in cmd_sequence:
    for filepath in filepaths:
        cmd = 'python biopython.py {} {}'.format(filepath, t)
        # print(cmd)
        os.system(cmd)
