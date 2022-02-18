rm -rf result/*.csv

python biopython.py db/Alpha.fasta 0
python biopython.py db/B.1.427_B.1.429.fasta 0
python biopython.py db/BA.2.fasta 0
python biopython.py db/Beta.fasta 0
python biopython.py db/Delta.fasta 0
python biopython.py db/Gamma.fasta 0
python biopython.py db/Mu.fasta 0
python biopython.py db/Omicron.fasta 0
python biopython.py db/Vum.fasta 0

python biopython.py db/Alpha.fasta 1
python biopython.py db/B.1.427_B.1.429.fasta 1
python biopython.py db/BA.2.fasta 1
python biopython.py db/Beta.fasta 1
python biopython.py db/Delta.fasta 1
python biopython.py db/Gamma.fasta 1
python biopython.py db/Mu.fasta 1
python biopython.py db/Omicron.fasta 1
python biopython.py db/Vum.fasta 1
