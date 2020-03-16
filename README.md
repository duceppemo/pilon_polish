## pilon_polish
Polish assembly (fasta) with Illumina paired-end reads (fastq)

## Usage
```
usage: pilon_polish.py [-h] -i /input_folder/ -f /input_fastq/ -o
                       /output_folder/ [-t 48] [-m 459]

Polish genome assembly file (fasta) with Illumina paired-end reads (fastq)

optional arguments:
  -h, --help            show this help message and exit
  -i /input_folder/, --input /input_folder/
                        Folder that contains the genome assemblies to polish
  -f /input_fastq/, --fastq /input_fastq/
                        Folder that contains the matching fastq files. File
                        names before first "_" must beunique and matching
                        assembly names (also before first "_"
  -o /output_folder/, --output /output_folder/
                        Folder to save the polished assemblies
  -t 48, --threads 48   Number of CPU. Default is maximum CPU available(48)
  -m 459, --memory 459  Memory to use for Pilon (Java). Default is 85%
                        available memory (459GB)
```
