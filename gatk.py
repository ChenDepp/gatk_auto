#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Aug  3 21:20:17 2020
@author: gchen

Changed on 10:23 2021 7 20
"""


#global variants 
Version = "0.2"

try:
    import logging 
    import os 
    import sys 
    import argparse 
    import textwrap
    import distutils.spawn
    import subprocess as sbp
    import multiprocessing as mp
    from abc import ABC, abstractmethod
except ImportError as e:
    raise ImportError(f"{e.name} package not found, please install it before running this program!")


#check python3 version 
def check_python_version(version:tuple):
    if sys.version_info[:2] < version:
        sys.stderr.write(f"the python version is less than {version}\n")
        sys.exit(1)

class DirectoryError(Exception):
    """
    when directory don't exist, raise this error!
    """

class FileError(Exception):
    """
    when file don't exist, raise this error!
    """

class RunSystemCmdError(sbp.CalledProcessError):
    """
    when system command running failed, raise this error!
    """

#检查文件是否存在
def check_file_exists(filenames):
    if isinstance(filenames, (tuple, list)):
        store_return = [] 
        for file in filenames:
            file = os.path.expanduser(file) #避免使用者使用~指代家目录时找不到路径
            if not os.path.isfile(file):
                raise FileError(f"{file} not is file or don't exists!")
            store_return.append(os.path.abspath(file))
        return store_return 
    elif isinstance(filenames, str):
        filenames = os.path.expanduser(filenames)
        if not os.path.isfile(filenames):
            raise FileError(f"{filenames} not is file or don't exist!")
        return os.path.abspath(filenames)
    else:
        raise TypeError('please check input!')

#检查文件夹是否存在
def check_dir_exists(dirnames):
    if isinstance(dirnames, (tuple, list)):
        store_return = []
        for dirname in dirnames:
            dirname = os.path.expanduser(dirname)
            if not os.path.isdir(dirname):
                raise DirectoryError(f"{dirname} not is directory or don't exist!")
            store_return.append(os.path.abspath(dirname))
        return store_return
    elif isinstance(dirnames, str):
        dirnames = os.path.expanduser(dirnames)
        if not os.path.isdir(dirnames):
            raise DirectoryError(f"{dirnames} not is directory or don't exist!")
        return os.path.abspath(dirnames)
    else:
        raise TypeError('please check input!')
    
#当输出目录不存在的时候，新建输出目录 
def check_outdir(dirname):
    dirname = os.path.expanduser(dirname)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return os.path.abspath(dirname)

#检测软件是否存在 
def check_software_exists(software):
    if isinstance(software, str):
        assert distutils.spawn.find_executable(software), f"{software} not exists in you environment!"
    elif isinstance(software, (list, tuple)):
        for i in software:
            assert distutils.spawn.find_executable(i), f"{i} not exists in you evironment!"
            
            
def prepare_argparser():
    
    description = textwrap.dedent("""
                    The purpose of this program is for calling snp used gatk!
                    --------------------------------------------------------------------""")
    
    epilog = "For command line option of each command, type: %(prog)s COMMAND -h"
    
    parser = argparse.ArgumentParser(description=description, 
                                  formatter_class=argparse.RawDescriptionHelpFormatter,
                                  epilog=epilog)
    parser.add_argument('-v', "--version", action="version", version="%(prog)s " + Version)
    
    parser.add_argument("-i", "--indir", type=check_dir_exists, required=True,
                                help="the directory where bam located", metavar="")
    parser.add_argument("-o", "--outdir", type=check_outdir, required=True,
                                help="the directory save the snp and indel calling result", metavar="")
    parser.add_argument("-m", "--mode", type=str, choices=("pc",'sc'), default="sc",
                                help="the calling mode, pc: population calling, sc: single calling [sc]", metavar="")
    gatk_split_chroms_method_group = parser.add_mutually_exclusive_group()
    gatk_split_chroms_method_group.add_argument("-sn", "--split-num", type=int,
                                help="the number of split genome", metavar="")
    gatk_split_chroms_method_group.add_argument("-ss", "--split-size", type=int,
                                help="the size for split genome [unit M]", metavar="")
    parser.add_argument('--split-mode', type=str, choices=("imprecise", "precise"), default="imprecise",
                                help="the mode of split genome [imprecise]", metavar="")
    parser.add_argument('-hcjo', '--haplotypecaller-java-options', type=str,
                                help="the java options for HaplotypeCaller, such as -Xmx100G limit the memory usage", metavar="")
    parser.add_argument('-cgjo', '--combinegvcf-java-options', type=str,
                                help="the java options for CombineGVCFs", metavar="")
    parser.add_argument('-ggjo', '--genotypegvcf-java-options', type=str,
                                help="the java options for GenotypeGVCFs", metavar="")
    parser.add_argument('-mvjo', '--mergevcf-java-options', type=str,
                                help="the java options for MergeVcfs", metavar="")
    parser.add_argument('--gvcf', action="store_true",
                                help="output gvcf file , default output is vcf file [False]")
    parser.add_argument("-s", "--suffix-name", nargs="+", type=str, default=["sort.fixmate.rmdup.bam"],
                                help="the suffix name of input file [default sort.fixmate.rmdup.bam]", metavar="")
    parser.add_argument("-p", "--process", type=int, default=5,
                                help="the porcess number for deal with multiple samples [5]", metavar="")
    parser.add_argument('-pbi', '--process-build-index', type=int, default=5,
                                help="the process of build bam index [5]", metavar="")
    parser.add_argument('-t', '--threads', type=int, default=5,
                                help="the threads for index and sort bam file [5]", metavar="")
    parser.add_argument("-r", "--ref", type=check_file_exists, required=True,
                                help="the reference genome", metavar="")
    parser.add_argument("-n", "--prefix-name", type=str, default="All", 
                                help="the prefix name of output vcf file [All]", metavar="")
    parser.add_argument('-l', '--log', type=str, default="gatk.log", 
                                help="the log file for record running information [gatk.log]", metavar="")

    return parser 



#single end files
class Single_end(object):
    def __init__(self, dirname, suffix_name=None):
        self.dirname = dirname
        if isinstance(suffix_name, list):
            suffix_name = tuple(suffix_name)
        self.suffix_name = suffix_name
        if self.suffix_name:
            self._files = [os.path.join(self.dirname, file) for file in os.listdir(dirname) if file.endswith(suffix_name) and os.path.isfile(os.path.join(self.dirname, file))]
        else:
            self._files = [os.path.join(self.dirname, file) for file in os.listdir(dirname) if os.path.isfile(os.path.join(self.dirname, file))]
        if not self._files:
            raise ValueError(f"the directory {dirname} is empty or there is no file with the special suffix-name")
    def sort_by_size(self):
        self._files.sort(key=lambda file:os.stat(file).st_size)
        return self
    def __getitem__(self, index):
        return self._files[index]
    def __len__(self):
        return len(self._files)
    def __repr__(self):
        return f'{self.__class__.__name__}{(self.dirname, self.suffix_name)}'
    def __str__(self):
        return os.path.join(self.dirname, f"*{self.suffix_name if self.suffix_name else ''}")
    def __iter__(self):
        return iter(self._files)
    
    
#基类，统一调用接口
class Base_Class(ABC):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if "log" in self.__dict__:
            self.set_logger(self.log)
    ###concatenate command
    @staticmethod
    def concat_cmd(main_cmd, preset_part, *link_part, **kwargs):
        """
        main_cmd:一般指定运行的程序名
        preset_part:中间预设的参数部分
        link_part:后续链接的参数部分
        kwargs:一些参数
        """
        for key, value in kwargs.items():
            if value:
                if len(key) == 1:
                    main_cmd += f" -{key}"
                else:
                    main_cmd += f" --{key}"
                if not isinstance(value, bool):
                    main_cmd += f" {value}"
        link_part_cmd = ""
        for cmd in link_part:
            if cmd:
                link_part_cmd += f" && {cmd}"
        return main_cmd + ' ' + preset_part + link_part_cmd
    @staticmethod
    def RunSystemCmd(cmd, wkdir='.'): #调用系统命令
        try:
            sbp.run(cmd, shell=True, check=True, cwd=wkdir)
        except sbp.CalledProcessError as e:
            raise RunSystemCmdError(f"{e.args} running failed")
    def set_logger(self, filename, level=logging.DEBUG):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(level=level) 
        filehandler = logging.FileHandler(filename, mode='w', encoding='utf-8')
        stdouthandler = logging.StreamHandler()
        filehandler_formater = logging.Formatter("%(asctime)s:%(name)s:%(levelname)s:%(message)s", 
                                datefmt='%Y/%m/%d %I:%M:%S %p')
        filehandler.setFormatter(filehandler_formater)
        stdouthandler.setFormatter(filehandler_formater)
        self.logger.addHandler(filehandler)
        self.logger.addHandler(stdouthandler)
    @abstractmethod
    def __call__(self):
        """
        for same interface!
        """
        raise NotImplementedError("this is abstractmethod of base class, please implement the method in childern class!")
    


###GATK
class GATK(Base_Class):
    build_bam_index_cmd = "samtools index -%s -@ {threads} {file}"
    build_reference_index_cmd = "samtools faidx {reference}"
    build_reference_dict_cmd = "gatk CreateSequenceDictionary -R {reference}"
    gatk_pc_variant_call_cmd = "gatk %sHaplotypeCaller -R {reference} -I {file} -O {prefix_name}.vcf"
    gatk_pc_variant_call_by_chrom_cmd = "gatk %sHaplotypeCaller -R {reference} -I {file} -L {chroms} -O {outfile}"
    gatk_sc_variant_call_cmd = "gatk %sHaplotypeCaller --ERC GVCF -R {reference} -I {file} -O {prefix_name}.gvcf"
    gatk_sc_variant_call_by_chrom_cmd = "gatk %sHaplotypeCaller --ERC GVCF -R {reference} -I {file} -L {chroms} -O {outfile}"
    gatk_sc_combine_gvcf_cmd = "gatk %sCombineGVCFs -R {reference} -V {file} -O {prefix_name}.gvcf"
    gatk_sc_genotype_gvcf_cmd = "gatk %sGenotypeGVCFs -R {reference} -V {prefix_name}.gvcf -O {prefix_name}.vcf && rm {prefix_name}.gvcf {prefix_name}.gvcf.idx"
    gatk_merge_vcf_cmd = "gatk %sMergeVcfs -I {file} -O {prefix_name}.vcf"
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        #set java options 
        java_options_map_relation = {"gatk_pc_variant_call_cmd": self.haplotypecaller_java_options, "gatk_pc_variant_call_by_chrom_cmd": self.haplotypecaller_java_options,
                                        "gatk_sc_variant_call_cmd": self.haplotypecaller_java_options, "gatk_sc_variant_call_by_chrom_cmd": self.haplotypecaller_java_options,
                                            "gatk_sc_combine_gvcf_cmd": self.combinegvcf_java_options, "gatk_sc_genotype_gvcf_cmd": self.genotypegvcf_java_options,
                                                "gatk_merge_vcf_cmd": self.mergevcf_java_options}
        for var_name in java_options_map_relation:
            if java_options_map_relation[var_name]:
                self.__dict__[var_name] = GATK.__dict__[var_name] % (f"--java-options '{java_options_map_relation[var_name]}' ")
            else:
                self.__dict__[var_name] = GATK.__dict__[var_name] % ("")    
        #running each chromosome
        self.split_running = self.split_num or self.split_size
        if self.gvcf and self.split_running:
            sys.stderr.write(f"can't split chromosome calling in --gvcf mode\n")
            sys.exit(1)
    def remove_files(self, files, *args, outdir=None):   #args set suffix name of file (need to be remove)
        for file in files:
            if isinstance(file, (list, tuple, set)):
                self.remove_files(file, *args, outdir=outdir)
            else:
                if not outdir is None:
                    file = os.path.join(outdir, file)
                if os.path.isfile(file):
                    os.remove(file)
                for suffix_name in args:
                    file = f"{file}.{suffix_name}"
                    if os.path.isfile(file):
                        os.remove(file)
    def check_index_format(self): #check index format based on fasta file or fai index file
        thershold_size = 2 ** 29
        self.index_format = "bai"
        with open(f"{self.ref}.fai") as f:
            for line in f:
                chrom_len = int(line.split()[1])
                if chrom_len > thershold_size:
                    self.index_format = "csi"
                    return 
    def build_bam_index(self): #when bam file with index file, skip build bam index file
        bam_files = []
        for file in Single_end(self.indir, self.suffix_name):
            if not os.path.isfile(f"{file}.bai") and not os.path.isfile(f"{file}.csi"):
                bam_files.append(file)
        if bam_files: #not empty
            self.check_index_format()
            self.build_bam_index_cmd = self.build_bam_index_cmd % ("-b" if self.index_format == "bai" else "-c")
            with mp.Pool(self.process_build_index) as P:
                for file in bam_files:
                    try:
                        self.logger.info(f"start build bam index for {file}")
                        P.apply_async(self.RunSystemCmd, args=(self.build_bam_index_cmd.format(**vars(), threads=self.threads),))
                    except sbp.CalledProcessError as e:
                        self.logger.error(f"{e.cmd} running failed")
                P.close()
                P.join() 
    @property
    def check_ref_fai(self):
        if os.path.isfile(f"{self.ref}.fai"):
            return True 
        return False 
    @property
    def check_ref_dict(self):
        if os.path.isfile(f"{self.ref.rsplit('.', 1)[0]}.dict"):
            return True 
        return False 
    def build_ref_fai(self):
        if not self.check_ref_fai:
            self.logger.info(f"start build reference genome fai index for {os.path.basename(self.ref)}")
            self.RunSystemCmd(self.build_reference_index_cmd.format(reference=self.ref))
    def build_ref_dict(self):
        if not self.check_ref_dict:
            self.logger.info(f"start build reference genome dict for {os.path.basename(self.ref)}")
            self.RunSystemCmd(self.build_reference_dict_cmd.format(reference=self.ref))
    def genome_split(self):
        fai_file = f"{self.ref}.fai"
        chrom_info = []
        with open(fai_file) as f:
            for line in f:
                chrom_name, chrom_size = line.split()[:2]
                chrom_info.append((chrom_name, int(chrom_size)))
        if self.split_num: #split by num
            total_size = sum(x[1] for x in chrom_info)
            split_size = total_size // self.split_num
        else:       #split by size
            split_size = self.split_size * 10 ** 6
        #imprecise split mode
        if self.split_mode == "imprecise":
            chrom_name_lists, record_size = [], 0
            for chrom_name, chrom_size in chrom_info:
                chrom_name_lists.append(f"'{chrom_name}'")
                record_size += chrom_size
                if record_size >= split_size:
                    yield chrom_name_lists
                    chrom_name_lists, record_size = [], 0
            if chrom_name_lists:
                yield chrom_name_lists
        else: #precise split mode
            chroms_info = chrom_info[:]
            #record region len // #store region info
            remain_size, store_region_info = split_size, []
            while chroms_info:
                #get one region
                region = chroms_info[0]
                if store_region_info and remain_size == 0:
                    yield store_region_info
                    remain_size, store_region_info = split_size, []
                if len(region) == 2:
                    if region[1] >= remain_size:
                        store_region_info.append(f"'{region[0]}':{1}-{remain_size}")
                        if region[1] != remain_size:
                            chroms_info[0] = (region[0], remain_size + 1, region[1])
                        else:
                            chroms_info.pop(0)
                        remain_size = 0
                    else:
                        store_region_info.append(f"'{region[0]}':{1}-{region[1]}")
                        chroms_info.pop(0)
                        remain_size -= region[1]
                else:
                    region_size = region[2] - region[1] + 1
                    if region_size >= remain_size:
                        region_end = region[1] + remain_size
                        store_region_info.append(f"'{region[0]}':{region[1]}-{region_end}")
                        if region_size != remain_size:
                            chroms_info[0] = (region[0], region_end + 1, region[2])
                        else:
                            chroms_info.pop(0)
                        remain_size = 0
                    else:
                        store_region_info.append(f"'{region[0]}':{region[1]}-{region[2]}")
                        chroms_info.pop(0)
                        remain_size -= region_size
            yield store_region_info
    @staticmethod
    def check_gvcf_idx_timestamp(gvcf_file):
        #check gvcf_idx file timestamp, gvcf idx file creation date should be after gvcf file creation date
        if os.path.isfile(gvcf_file) and os.path.isfile(f"{gvcf_file}.idx"):
            if os.path.getctime(gvcf_file) > os.path.getctime(f"{gvcf_file}.idx"):
                sys.stderr.write(f"gvcf idx file {gvcf_file}.idx creation date should be after gvcf file {gvcf_file} creation date\n")
                sys.exit(1)
    def __call__(self):
        self.build_ref_fai()
        self.build_ref_dict()
        self.build_bam_index()
        files = Single_end(self.indir, self.suffix_name)       
        if not self.split_running:
            if self.mode == "pc":
                self.logger.info("start variant detect with population calling strategy")
                file = " -I ".join(files)
                if self.gvcf:
                    self.RunSystemCmd(self.gatk_sc_variant_call_cmd.format(**vars(), reference=self.ref, prefix_name=self.prefix_name), self.outdir)
                else:
                    self.RunSystemCmd(self.gatk_pc_variant_call_cmd.format(**vars(), reference=self.ref, prefix_name=self.prefix_name), self.outdir)
            else:
                if len(files) == 1:
                    sys.stderr.write(f"there only one sample, please don't use sc calling mode!\n")
                    sys.exit(1)
                self.logger.info("start variant detect with single calling strategy")
                with mp.Pool(self.process) as P:
                    gvcf_files = []
                    for file in files:
                        file_prefix_name = os.path.basename(file).rsplit('.', 1)[0]
                        if not os.path.isfile(os.path.join(self.outdir, f"{file_prefix_name}.gvcf.idx")):
                            P.apply_async(self.RunSystemCmd, args=(self.gatk_sc_variant_call_cmd.format(**vars(), prefix_name=file_prefix_name, reference=self.ref), self.outdir))
                        else:
                            self.check_gvcf_idx_timestamp(os.path.join(self.outdir, f"{file_prefix_name}.gvcf")) #check gvcf_idx timestamp
                        gvcf_files.append(f"{file_prefix_name}.gvcf")
                    P.close()
                    P.join() 
                self.logger.info("start combine gvcf")
                file = " -V ".join(gvcf_files)
                self.RunSystemCmd(self.gatk_sc_combine_gvcf_cmd.format(**vars(), reference=self.ref, prefix_name=self.prefix_name), self.outdir)
                self.logger.info("combine gvcf successfully")
                #remove each sample gvcf 
                self.remove_files(gvcf_files, "idx", outdir=self.outdir)
                #start genotypeing gvcf file
                if not self.gvcf:
                    self.logger.info("start genotypeing gvcf file")
                    self.RunSystemCmd(self.gatk_sc_genotype_gvcf_cmd.format(reference=self.ref, prefix_name=self.prefix_name), self.outdir)
        else:
            if self.mode == "pc":
                file = " -I ".join(files)
                vcf_files = []
                with mp.Pool(self.process) as P:
                    for i, chroms_list in enumerate(self.genome_split(), 1):
                        chroms = " -L ".join(chroms_list)
                        outfile = f"{self.prefix_name}.batch{i}.vcf"
                        #skip some batch
                        if os.path.isfile(f"{outfile}.idx"):
                            sys.stdout.write(f"it will skip produce the gvcf file {outfile}, but note that the splitting related parameters must be the same as last time\n")
                            self.check_gvcf_idx_timestamp(os.path.join(self.outdir, f"{outfile}")) #check gvcf_idx timstamp
                        else:
                            P.apply_async(self.RunSystemCmd, args=(self.gatk_pc_variant_call_by_chrom_cmd.format(**vars(), reference=self.ref), self.outdir))
                        vcf_files.append(outfile)
                    P.close()
                    P.join()
                self.logger.info("start merge vcf file")
                file = " -I ".join(vcf_files)
                self.RunSystemCmd(self.gatk_merge_vcf_cmd.format(**vars(), prefix_name=self.prefix_name), self.outdir)
                self.logger.info("merge vcf file successfully")
                #remove vcf of each batch
                self.remove_files(vcf_files, "idx", outdir=self.outdir)
            else:
                #for get gvcf of each chromsome for all samples
                if len(files) == 1:
                    sys.stderr.write(f"there only one sample, please don't use sc calling mode!\n")
                    sys.exit(1)
                with mp.Pool(self.process) as P:
                    gvcf_files = {}
                    for i, chroms_list in enumerate(self.genome_split(), 1):
                        chroms = " -L ".join(chroms_list)
                        for file in files:
                            prefix_name = os.path.basename(file).rsplit('.', 1)[0] 
                            outfile = f"{prefix_name}.batch{i}.gvcf"
                            #skip some batch 
                            if os.path.isfile(f"{outfile}.idx"):
                                sys.stdout.write(f"it will skip produce the gvcf file {outfile}, but note that the splitting related parameters must be the same as last time\n")
                                self.check_gvcf_idx_timestamp(os.path.join(self.outdir, f"{outfile}"))
                            else:
                                P.apply_async(self.RunSystemCmd, args=(self.gatk_sc_variant_call_by_chrom_cmd.format(**vars(), reference=self.ref), self.outdir))
                            gvcf_files.setdefault(i, []).append(outfile)
                    P.close()
                    P.join()
                #combine gvcf of each chromsome for all samples
                self.logger.info(f"start combine and genotypeing gvcf")
                with mp.Pool(self.process) as P:
                    prefix_names = []
                    for i, gvcf_file in gvcf_files.items():
                        file = " -V ".join(gvcf_file)
                        prefix_name = f"{self.prefix_name}.batch{i}"
                        #setting skip code
                        if os.path.isfile(f"{prefix_name}.vcf.idx"):
                            sys.stdout.write(f"it will skip produce the vcf file {prefix_name}.vcf, but note that the splitting related parameters must be the same as last time\n")
                            self.check_gvcf_idx_timestamp(os.path.join(self.outdir, f"{prefix_name}.vcf"))
                        elif os.path.isfile(f"{prefix_name}.gvcf.idx"):
                            sys.stdout.write(f"it will skip produce the gvcf file {prefix_name}.gvcf, but note that the splitting related parameters must be the same as last time\n")
                            self.check_gvcf_idx_timestamp(os.path.join(self.outdir, f"{prefix_name}.gvcf"))
                            P.apply_async(self.RunSystemCmd, args=(self.gatk_sc_genotype_gvcf_cmd.format(**vars(), reference=self.ref), self.outdir))
                        else:
                            P.apply_async(self.RunSystemCmd, args=(" && ".join([self.gatk_sc_combine_gvcf_cmd, self.gatk_sc_genotype_gvcf_cmd]).format(**vars(), reference=self.ref), self.outdir))
                        prefix_names.append(prefix_name)
                    P.close()
                    P.join()
                self.logger.info("combine and genotypeing gvcf successfully")
                #remove gvcf
                self.remove_files(list(gvcf_files.values()), "idx", outdir=self.outdir)
                #merge vcf
                self.logger.info("start merge vcf file")
                vcf_files_chrom = [f"{prefix_name}.vcf" for prefix_name in prefix_names]
                file = " -I ".join(vcf_files_chrom)
                self.RunSystemCmd(self.gatk_merge_vcf_cmd.format(file=file, prefix_name=self.prefix_name), self.outdir)
                self.logger.info("merge vcf file successfully")
                #remove vcf
                self.remove_files(vcf_files_chrom, "idx", outdir=self.outdir)
        self.logger.info("variant calling successfully")
    
    
def main():
    parser = prepare_argparser()
    args = parser.parse_args()
    sys.stdout.write('Notice:if you want set the java-options, reference the write style, such as -hcjo "\'-Xmx100G\'"\n')
    check_software_exists(["gatk", "samtools"])
    GATK(**args.__dict__)()
    
if __name__ == "__main__":
    check_python_version((3, 8))
    main()

