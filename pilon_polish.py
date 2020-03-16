#!/usr/local/env python3

from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from shutil import which
from glob import glob
import os
import pathlib
import subprocess


class SampleObject(object):
    pass


class PilonPolish(object):
    fasta_ext = ['.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz']  # TODO -> check if can use gzipped fasta
    fastq_ext = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']

    def __init__(self, args):
        self.input_fasta_folder = args.input
        self.input_fastq_folder = args.fastq
        self.output_folder = args.output
        self.cpu = args.threads
        self.mem = args.memory

        # data structure
        self.sample_dict = dict()

        # run script
        self.run()

    def run(self):
        # Get fasta and fastq file(s) info
        PilonPolish.get_fasta(self.input_fasta_folder, PilonPolish.fasta_ext, self.sample_dict)
        PilonPolish.get_fastq(self.input_fastq_folder, PilonPolish.fastq_ext, self.sample_dict)
        self.pilon_polish(self.sample_dict, self.cpu, self.mem, self.output_folder)  # TODO -> make parallel

    @staticmethod
    def get_fasta(folder, ext_list, sample_dict):
        file_list = list()
        for ext in ext_list:
            files = glob(folder + '/*' + ext)
            if files:
                file_list.extend(files)

        for f in file_list:
            sample_name = os.path.basename(f).split('_')[0].split('.')[0]
            sample_object = SampleObject()
            sample_object.fasta = f
            sample_object.name = sample_name
            sample_dict[sample_name] = sample_object

    @staticmethod
    def get_fastq(folder, ext_list, sample_dict):
        file_list = list()
        for ext in ext_list:
            files = glob(folder + '/*' + ext)
            if files:
                file_list.extend(files)

        for file in file_list:
            sample_name = os.path.basename(file).split('_')[0]  # It should have a "_" in the file name
            try:
                test = sample_dict[sample_name]
            except KeyError:
                raise Exception('Fastq and fasta files must be named the same.')
            try:
                if '_R1' in os.path.basename(file):
                    sample_dict[sample_name].fastq[0] = file
                elif '_R2' in os.path.basename(file):
                    sample_dict[sample_name].fastq[1] = file
            except AttributeError:
                setattr(sample_dict[sample_name], 'fastq', [None, None])
                if '_R1' in os.path.basename(file):
                    sample_dict[sample_name].fastq[0] = file
                elif '_R2' in os.path.basename(file):
                    sample_dict[sample_name].fastq[1] = file

    def check_files(self):
        pass

    @staticmethod
    def make_folder(folder):
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def pilon_polish(sample_dict, cpu, mem, output_folder):
        # do 5 rounds of polishing

        for name, sample_obj in sample_dict.items():
            for i in range(1, 5 + 1):
                if i == 1:  # the first polishing round is from the assembly
                    bam = PilonPolish.map_reads(sample_obj.fasta, sample_obj.fastq[0], sample_obj.fastq[1], cpu,
                                                output_folder)
                    PilonPolish.run_pilon(sample_obj.fasta, bam, cpu, mem, output_folder, i)
                else:  # use the output of the previous round as input
                    # Check if there were changes in the previous pilon round
                    previous_change_file = '{}/{}/{}_pilon{}.changes'.format(output_folder, name, name, i - 1)
                    previous_fasta_file = '{}/{}/{}_pilon{}.fasta'.format(output_folder, name, name, i - 1)
                    if os.stat(previous_change_file).st_size == 0:  # If no changes, stop here
                        break
                    else:
                        bam = PilonPolish.map_reads(previous_fasta_file, sample_obj.fastq[0], sample_obj.fastq[1], cpu,
                                                    output_folder)
                        PilonPolish.run_pilon(previous_fasta_file, bam, cpu, mem, output_folder, i)
                        os.remove(previous_fasta_file)  # delete previous fasta file

    @staticmethod
    def map_reads(fasta, r1, r2, cpu, output_folder):
        name = os.path.basename(fasta).split('_')[0].split('.')[0]
        # Create output directory
        PilonPolish.make_folder(output_folder + '/' + name)
        # Output bam file
        out_bam = '{}/{}/{}.bam'.format(output_folder, name, name)

        # Commands to run
        index_cmd = ['bwa', 'index', fasta]
        mapping_cmd = ['bwa', 'mem', '-t', str(cpu), '-M', fasta, r1, r2]
        samtools_view_cmd = ['samtools', 'view', '-@', str(cpu), '-bhF', '4', '-']
        samtools_sort_cmd = ['samtools', 'sort', '-@', str(cpu), '-']
        samtools_rmdup_cmd = ['samtools', 'rmdup', '-', out_bam]
        samtools_index_cmd = ['samtools', 'index', out_bam]

        # Index fasta for bwa
        subprocess.run(index_cmd)

        # Run bwa and samtools using pipes
        p1 = subprocess.Popen(mapping_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)
        p2.stdout.close()
        p4 = subprocess.Popen(samtools_rmdup_cmd, stdin=p3.stdout)
        p3.stdout.close()
        p4.communicate()

        # Index bam file
        subprocess.run(samtools_index_cmd)

        # Cleanup
        to_clean = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        for ext in to_clean:
            os.remove(fasta + ext)  # Index files

        return out_bam

    @staticmethod
    def run_pilon(fasta, bam, cpu, mem, output_folder, i):
        # Find where pilon binaries are installed
        pilon_path = which('pilon')
        pilon_path = '/'.join(pilon_path.split('/')[:-2])
        pilon_path += '/share/pilon-1.23-2/pilon-1.23.jar'  # this works if pilon was installed with conda

        # Output name considering which iteration we're at
        name = os.path.basename(fasta).split('_')[0].split('.')[0]
        output_name = '{}_pilon{}'.format(name, i)

        pilon_cmd = ['java', '-Xmx{}g'.format(mem), '-jar', pilon_path,
                     '--threads', str(cpu),
                     '--genome', fasta,
                     '--frags', bam,
                     '--outdir', '{}/{}'.format(output_folder, name),
                     '--output', output_name,
                     '--changes']
        p1 = subprocess.Popen(pilon_cmd)
        (stdout, stderr) = p1.communicate()

        # Cleanup files
        to_clean = ['.bam', '.bam.bai']
        for ext in to_clean:
            os.remove('{}/{}/{}{}'.format(output_folder, name, name, ext))  # Index files


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Polish genome assembly file (fasta) with Illumina paired-end reads (fastq)')
    parser.add_argument('-i', '--input', metavar='/input_folder/',
                        required=True,
                        help='Folder that contains the genome assemblies to polish')
    parser.add_argument('-f', '--fastq', metavar='/input_fastq/',
                        required=True,
                        help='Folder that contains the matching fastq files. File names before first "_" must be '
                             'unique and matching assembly names (also before first "_"')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        help='Folder to save the polished assemblies')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU. Default is maximum CPU available({})'.format(max_cpu))
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory to use for Pilon (Java). Default is 85%% available memory ({}GB)'.format(max_mem))

    # Get the arguments into an object
    arguments = parser.parse_args()
    PilonPolish(arguments)