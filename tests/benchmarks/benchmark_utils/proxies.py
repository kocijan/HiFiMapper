from pathlib import Path
import os
import time
import random
import re
from typing import List, Dict
import hashlib
import inspect


WORK_DIRECTORY = Path(os.path.dirname(__file__))
EXE_DIRECTORY = WORK_DIRECTORY.parents[2].joinpath("build/bin")
RESULTS_DIRECTORY = WORK_DIRECTORY.parents[0].joinpath("results")
DATA_DIRECTORY = WORK_DIRECTORY.parents[0].joinpath("data")


os.makedirs(RESULTS_DIRECTORY, exist_ok=True)
os.makedirs(DATA_DIRECTORY, exist_ok=True)


def hash_str(s):
    return int(hashlib.sha1(s.encode("utf-8")).hexdigest(), 16)


def hash_file(full_name):
    content_hash_str = os.popen(f"git hash-object {full_name}").read()
    return hash((hash_str(full_name), hash_str(content_hash_str)))


def hash_folder(folder):
    all_hashes = []
    for dirName, _, fileList in os.walk(folder):
        for fname in fileList:
            full_name = f"{dirName}/{fname}"
            all_hashes.append(hash_file(full_name))
    return hash(tuple(all_hashes))


def hash_fixed(*argv):
    return hash(argv) % 10000000000000000


class ProxyObject:
    def __init__(self, _path):
        self.path = str(_path)

    def __str__(self):
        return str(self.path)

    def __hash__(self):
        return hash_str(self.path)


class Reference(ProxyObject):
    @staticmethod
    def generate(
        reference_size: int,
        seed: int = None,
        probability: float = 0.01,
        name: str = "ARTIFICIAL_REF",
        instructions: Dict[int, int] = None,
    ):
        if seed is None:
            seed = random.randint(0, 1000)

        reference_generator = EXE_DIRECTORY.joinpath("Refgen")

        if instructions is None:
            instructions = {}

        r = len(instructions.keys())
        add_seeds = tuple([i for item in instructions.items() for i in item])
        cons_hash = hash_fixed(
            hash_str("Reference"),
            hash_str(str(reference_generator)),
            reference_size,
            seed,
            probability,
            hash_str(name)
        )

        ref_name = str(cons_hash) + ".fasta"
        path = DATA_DIRECTORY / ref_name

        if not path.exists():
            path.touch()
            instructions_str = str(add_seeds)[1:-1].replace(",", "")
            os.system(
                f"{reference_generator} -r {r} -s {seed} -p {probability} -n {name} \
                {reference_size} {instructions_str} > {path}"
            )

        return Reference(path)


class Reads(ProxyObject):
    @staticmethod
    def generate(
        reference_path: str,
        read_length: int = 25000,
        num_reads: int = 10000,
        complement: float = 0.5,
        error_probability: float = 0.01,
        seed: int = None,
    ):
        if seed is None:
            seed = random.randint(0, 1000)

        reads_generator = EXE_DIRECTORY.joinpath("Readsgen")

        cons_hash = hash_fixed(
            hash_str("Reads"),
            hash_str(str(reference_path)),
            hash_str(str(reads_generator)),
            read_length,
            num_reads,
            complement,
            error_probability,
            seed,
        )

        reads_name = str(cons_hash) + ".fastq"
        path = DATA_DIRECTORY / reads_name
        if not path.exists():
            path.touch()
            command = f"{reads_generator} -l {read_length} -n {num_reads} -c {complement} -p {error_probability} -s {seed} {reference_path} > {path}"
            os.system(command)

        return Reads(path)

    @staticmethod
    def simulate(
        reference_path: str,
        subreads_path: str,
        alignment_path: str,
        read_length: int = 25000,
        num_reads: int = 10000,
        complement: float = 0.5,
        seed: int = None,
    ):
        if seed is None:
            seed = random.randint(0, 1000)

        reads_generator = EXE_DIRECTORY.joinpath("Readsgen")

        cons_hash = hash_fixed(
            hash_str("Reads"),
            hash_str(reference_path),
            hash_str(subreads_path),
            hash_str(alignment_path),
            hash_str(str(reads_generator)),
            read_length,
            num_reads,
            complement,
            seed,
        )

        reads_name = str(cons_hash) + ".fastq"
        path = DATA_DIRECTORY / reads_name
        if not path.exists():
            path.touch()
            os.system(
                f"{reads_generator} -m -l {read_length} -n {num_reads} -c {complement} -s {seed} \
                {reference_path} {subreads_path} {alignment_path} > {path}"
            )

        return Reads(path)


class ResultsHiFiMapper(ProxyObject):
    @staticmethod
    def create(
        reference: str,
        reads: str,
        threads: int = 8,
        sample_length=75,
        sample_count=20,
        seed: int = None,
        min_match:float = 0.9,
        quality: int = 90,
        frequency: int = 10,
        discard: bool = False,
        extended_search: bool = False,
        sequential: bool = False,
        lcp_information: bool = False,
        fix_unnmaped: bool = False,
        minimal_anchors: int = 3,
        bandwidth: int = 10,
        gap: int = 10000,
        secondary_alignements: int = 5,
        secondary_to_primary_ratio: float = 0.8,
        lcp_search_size: int = 100,
        DEBUG=False,
        PROFILE=False,
    ):
        if seed is None:
            seed = random.randint(0, 1000)
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)

        abbrevations = ['t','l','c', 's', 'm', 'q', 'f', 'd', 'e', 'a', 'i', 'u', 'n', 'b', 'g', 'N', 'p','x']
        comm = "HiFiMapper"
        for i, value in enumerate(ResultsHiFiMapper.create.__defaults__):
            if value != values[args[i+2]]:
                if type(value) is bool:
                    comm += f" -{abbrevations[i]}"
                else:
                    comm += f" -{abbrevations[i]} {values[args[i+2]]}"

        mapper = EXE_DIRECTORY.joinpath("mapper")

        # TODO TODO TODO 
        mapper = "/home/mjkocijan/hifimapper/HiFiMapper2/build/bin/mapper"

        cons_hash = hash_fixed(
            hash_str("Results"),
            hash_str(str(mapper)),
            hash_str(str(reference)),
            hash_str(str(reads)),
                    threads,
                    sample_length,
                    sample_count,
                    seed,
                    min_match,
                    quality,
                    frequency,
                    discard,
                    extended_search,
                    sequential,
                    lcp_information,
                    fix_unnmaped,
                    minimal_anchors,
                    bandwidth,
                    gap,
                    secondary_alignements,
                    secondary_to_primary_ratio,
                    lcp_search_size
        )
        
        # reads = "/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.fixed.limited.fastq"
        # reads = "/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.fixed.fastq"
        results_name = str(cons_hash) + ".paf"
        descriptor_name = str(cons_hash) + "_stderr"
        path = RESULTS_DIRECTORY.joinpath(results_name)
        descriptor_path = RESULTS_DIRECTORY.joinpath(descriptor_name)

        command = f"{mapper} -t {threads}  -l {sample_length} -c {sample_count} -s {seed} -m {min_match} -q {quality} -f {frequency} -n {minimal_anchors} -b {bandwidth} -g {gap} -N {secondary_alignements} -p {secondary_to_primary_ratio} -x {lcp_search_size} "
        if discard:
            command += "-d "
        if extended_search:
            command += "-e "
        if sequential:
            command += "-a "
        if lcp_information:
            command += "-i "
        if fix_unnmaped:
            command += "-u "
        command += f"{reference} "
        if reads is not None:
            command += str(reads)
        print (f"{command} > {path} 2> {descriptor_path}")
        #if not path.exists(): # TODO
        os.system(f"{command} > {path} 2> {descriptor_path}")
        try:
            if DEBUG and PROFILE:
                os.system(
                    f"valgrind \
                            --trace-children=yes \
                            --vgdb-error=0 \
                            --tool=callgrind \
                            {command} > /dev/null\
                            & vgdb --wait=3 --pid=$! --port=8021"
                )
            if not DEBUG and PROFILE:
                os.system(
                    f"valgrind \
                            --trace-children=yes \
                            --tool=callgrind \
                            {command} > /dev/null"
                )
            if DEBUG and not PROFILE:
                os.system(
                    f"gdbserver \
                            :8021 \
                            {command} > /dev/null"
                )
        except Exception as e:
            print(e)

        results_hifi = ResultsHiFiMapper(path)
        results_hifi.descriptor_path = descriptor_path
        return results_hifi, comm

    def get_time(self):
        with open(self.descriptor_path, "r") as f:
            l = re.findall(r"[0-9]*[.][0-9]*", f.readlines()[-1])
            return l[0] if len(l) > 0 else 0
    
    def get_mapping_time(self):
        with open(self.descriptor_path, "r") as f:
            l = re.findall(r"[0-9]*[.][0-9]*", f.readlines()[-3])
            return l[0] if len(l) > 0 else 0
    
    def get_unmapped(self):
        with open(self.descriptor_path, "r") as f:
            l = re.findall(r"[0-9][0-9]*", f.readlines()[-2])
            return l[0] if len(l) > 0 else 0


class ResultsMinimap2(ProxyObject):
    @staticmethod
    def create(reference: str, reads: str, threads: int = 8, ssecondary_alignments = 5, min_ratio = 0.8):
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)

        abbrevations = ['t','N','p']
        comm = "minimap2 -x map-hifi"
        for i, value in enumerate(ResultsMinimap2.create.__defaults__):
            if value != values[args[i+2]]:
                comm += f" -{abbrevations[i]} {values[args[i+2]]}"

        mapper = EXE_DIRECTORY.parents[0].joinpath("minimap2").joinpath("minimap2")

        cons_hash = hash_fixed(
            hash_str("Results"), hash_str(str(reference)), hash_str(str(reads)), threads, ssecondary_alignments, min_ratio
        )

        results_name = str(cons_hash) + "_minimap2.paf"
        descriptor_name = str(cons_hash) + "_minimap2_stderr"
        path = RESULTS_DIRECTORY.joinpath(results_name)
        descriptor_path = RESULTS_DIRECTORY.joinpath(descriptor_name)

        if reads is None:
            reads = reference
        
        # reads = "/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.fixed.limited.fastq"
        # reads = "/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.fixed.fastq"
        
        command = f"{mapper} -x map-hifi -t {threads} -N {ssecondary_alignments} -p {min_ratio} {reference} {reads} > {path} 2> {descriptor_path}"
        print (command)
        #if not path.exists(): #TODO
        os.system(command)
        # else:
        #     print(f"# {mapper} -x map-hifi -k 5 -f {format(filter, '.10f')} -t {threads} -N {ssecondary_alignments} -p {min_ratio} {reference} {reads} > {path} 2> {descriptor_path}")

        results_minimap2 = ResultsMinimap2(path)
        results_minimap2.descriptor_path = descriptor_path
        return results_minimap2, comm

    def get_time(self):
        with open(self.descriptor_path, "r") as f:
            l = re.findall(r"[0-9]*[.][0-9]*", f.readlines()[-1])
            return l[0] if len(l) > 0 else 0


class ResultsWinnowmap(ProxyObject):
    @staticmethod
    def create(reference: str, reads: str, threads: int = 8, min_ratio = 0.8):
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)

        abbrevations = ['t','p']
        comm = "Winnowmap -x map-pb"
        for i, value in enumerate(ResultsWinnowmap.create.__defaults__):
            if value != values[args[i+2]]:
                comm += f" -{abbrevations[i]} {values[args[i+2]]}"

        mapper_bin_directory = (
            EXE_DIRECTORY.parents[0].joinpath("Winnowmap").joinpath("bin")
        )

        cons_hash = hash_fixed(
            hash_str("Results"), hash_str(str(reference)), hash_str(str(reads)), threads
        )

        results_name = str(cons_hash) + "_winnowmap.paf"
        descriptor_name = str(cons_hash) + "_winnowmap_stderr"
        path = RESULTS_DIRECTORY.joinpath(results_name)
        descriptor_path = RESULTS_DIRECTORY.joinpath(descriptor_name)

        if reads is None:
            reads = reference

        command = f"{mapper_bin_directory}/meryl count k=15 output \
        {mapper_bin_directory}/merylDB {reference} 2> /dev/null"
        #if not path.exists():
        os.system(command)
        print (command)
        command = f"{mapper_bin_directory}/meryl print greater-than distinct=0.9998 \
        {mapper_bin_directory}/merylDB > {mapper_bin_directory}/repetitive_k15.txt 2> /dev/null"
        #if not path.exists():
        os.system(command)
        print (command)
        command = f"{mapper_bin_directory}/winnowmap -W {mapper_bin_directory}/repetitive_k15.txt \
            -x map-pb -t {threads} {reference} {reads} > {path} 2> {descriptor_path}"
        #if not path.exists():
        # TODO
        os.system(command)
        print (command)

        results_winnowmap = ResultsWinnowmap(path)
        results_winnowmap.descriptor_path = descriptor_path
        return results_winnowmap, comm

    def get_time(self):
        with open(self.descriptor_path, "r") as f:
            l = re.findall(r"[0-9]*[.][0-9]*", f.readlines()[-1])
            return l[0] if len(l) > 0 else 0
