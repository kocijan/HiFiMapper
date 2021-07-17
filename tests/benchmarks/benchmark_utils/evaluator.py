from enum import Enum
import pyfastx
import matplotlib.pyplot as plt


class Strand(Enum):
    POSITIVE = "+"
    NEGATIVE = "-"


class PafLine:
    def __init__(self, line: str):
        self._content = line.split("\t")

    @property  
    def query_name(self):
        return str(self._content[0])
    
    @property
    def query_length(self):
        return int(self._content[1])
    
    @property
    def query_start(self):
        return int(self._content[2])
    
    @property
    def query_end(self):
        return int(self._content[3])

    @property
    def strand(self):
        return Strand.POSITIVE if str(self._content[4]) == "+" else Strand.NEGATIVE
    
    @property  
    def target_name(self):
        return str(self._content[5])
    
    @property
    def target_length(self):
        return int(self._content[6])
    
    @property
    def target_start(self):
        return int(self._content[7])
    
    @property
    def target_end(self):
        return int(self._content[8])

    @property
    def residue_matches(self):
        return int(self._content[9])
    
    @property
    def alignment_block_length(self):
        return int(self._content[10])
    
    @property
    def mapping_quality(self):
        return int(self._content[11])


class PafFile:
    def __init__(self, _path):
        self._path = _path
        self.lines_by_query = dict()
        self._parse_lines()
    
    def _parse_lines(self):
        with open(self._path, "r") as f:
            lines = f.readlines()
            for line in lines:
                pline = PafLine(line)
                if pline.query_name not in self.lines_by_query:
                    self.lines_by_query[pline.query_name] = []
                self.lines_by_query[pline.query_name].append(pline)


def _jaccard_sim(x1, y1, x2, y2):
    unx, uny = min(x1, x2), max(y1, y2)
    prx, pry = max(x1, x2), min(y1, y2)
    if pry < prx:
        return 0
    return float(pry - prx + 1) / (uny - unx + 1)


def _is_correct(start, end, true_start, true_end, percantage=0.1):
    if _jaccard_sim(start, end, true_start, true_end) >= percantage:
        return True
    return False


class PafEvaluator:
    def __init__(self, _path, reads, filter_mappings = False):
        self.paf_file = PafFile(_path)
        self._process(reads, filter_mappings)

    def jaccard_similarity(self):
        uk_qual = 0
        for qname in self.paf_file.lines_by_query:
            sp = qname.split("_")
            position = int(sp[-1].lstrip("0") + "0") // 10
            strain = Strand.POSITIVE if sp[7] == "+" else Strand.NEGATIVE
            best_result = 0
            for pline in self.paf_file.lines_by_query[qname]:
                if pline.strand != strain:
                    continue
                sim = _jaccard_sim(pline.target_start, pline.target_end, position, position + pline.query_length - 1)
                if sim > best_result:
                    best_result = sim
            uk_qual += best_result
        return uk_qual / len(self.paf_file.lines_by_query)
    
    def correctly_mapped(self):
        correct_total = 0
        total = 0
        for qname in self.paf_file.lines_by_query:
            sp = qname.split("_")
            position = int(sp[-1].lstrip("0") + "0") // 10
            strain = Strand.POSITIVE if sp[-3] == "+" else Strand.NEGATIVE
            for pline in self.paf_file.lines_by_query[qname]:
                total += 1
                if pline.strand != strain:
                    continue
                if _is_correct(pline.target_start, pline.target_end, position, position + pline.query_length - 1):
                    correct_total += 1
        if total != 0:
            return correct_total / total
        return 0
    
    def _process(self, reads: str, filter_mappings: bool):
        self.correct = 0
        self.wrong = 0
        self.unmapped = 0
        self.total = 0

        for qname in self.paf_file.lines_by_query:
            sp = qname.split("_")
            position = int(sp[-1].lstrip("0") + "0") // 10
            strain = Strand.POSITIVE if sp[-3] == "+" else Strand.NEGATIVE
            if filter_mappings:
                best = None
                for pline in self.paf_file.lines_by_query[qname]:
                    if best is None or pline.mapping_quality > best.mapping_quality:
                        best = pline
                self.paf_file.lines_by_query[qname] = [best]

            c = False
            for pline in self.paf_file.lines_by_query[qname]:
                self.total += 1
                if pline.strand != strain:
                    self.wrong += 1
                    continue
                if _is_correct(pline.target_start, pline.target_end, position, position + pline.query_length - 1):
                    c = True
                else:
                    self.wrong += 1
            if c:
                self.correct += 1

        fq = pyfastx.Fastx(reads)
        for name,seq,qual,comment in fq:
            if name not in self.paf_file.lines_by_query:
                self.unmapped += 1
    
    def draw_corrcet_coverage(self, reference_size: int, filename: str, ylim: int = None):
        coverage_list = [0 for i in range (reference_size)]

        for qname in self.paf_file.lines_by_query:
            sp = qname.split("_")
            position = int(sp[-1].lstrip("0") + "0") // 10
            strain = Strand.POSITIVE if sp[-3] == "+" else Strand.NEGATIVE

            for pline in self.paf_file.lines_by_query[qname]:
                if pline.strand != strain:
                    continue
                if _is_correct(pline.target_start, pline.target_end, position, position + pline.query_length - 1):
                    coverage_list[pline.target_start] += 1
                    coverage_list[pline.target_end if pline.target_end < reference_size else reference_size - 1] -= 1
        
        res = []
        suma = 0
        for i in range (len(coverage_list)):
            suma += coverage_list[i]
            res.append(suma)
            res.append(0)

        x = []
        for i in range(len(res)//2):
            x.append(i)
            x.append(i + 0.5)
        
        fig = plt.figure()
        if ylim is not None:
            plt.ylim(top=ylim)

        plt.plot(x, res)
        fig.savefig(filename)

    def draw_wrong_coverage(self, reference_size: int, filename: str, ylim: int = None):
        coverage_list = [0 for i in range (reference_size)]

        for qname in self.paf_file.lines_by_query:
            sp = qname.split("_")
            position = int(sp[-1].lstrip("0") + "0") // 10
            strain = Strand.POSITIVE if sp[-3] == "+" else Strand.NEGATIVE

            for pline in self.paf_file.lines_by_query[qname]:
                if pline.strand != strain:
                    coverage_list[pline.target_start] += 1
                    coverage_list[pline.target_end if pline.target_end < reference_size else reference_size - 1] -= 1
                    continue
                if not _is_correct(pline.target_start, pline.target_end, position, position + pline.query_length - 1):
                    coverage_list[pline.target_start] += 1
                    coverage_list[pline.target_end if pline.target_end < reference_size else reference_size - 1] -= 1
        
        res = []
        suma = 0
        for i in range (len(coverage_list)):
            suma += coverage_list[i]
            res.append(suma)
            res.append(0)

        x = []
        for i in range(len(res)//2):
            x.append(i)
            x.append(i + 0.5)
        
        fig = plt.figure()
        if ylim is not None:
            plt.ylim(top=ylim)

        plt.plot(x, res, color = 'red')
        fig.savefig(filename)


class RealEvaluator:
    def __init__(self, _path, reads, filter_mappings: bool = False):
        self.paf_file = PafFile(_path)
        self._process(reads, filter_mappings)
    
    def _process(self, reads, filter_mappings: bool):
        self.unmapped = 0
        self.total = 0
        fq = pyfastx.Fastx(reads)
        for name,seq,qual,comment in fq:
            self.total += 1
            if name not in self.paf_file.lines_by_query:
                self.unmapped += 1
        if filter_mappings:
            for qname in self.paf_file.lines_by_query:
                best = None
                for pline in self.paf_file.lines_by_query[qname]:
                    if best is None or pline.mapping_quality > best.mapping_quality:
                        best = pline
                self.paf_file.lines_by_query[qname] = [best]

    def draw_coverage(self, reference_size, filename, ylim: int = None):
        coverage_list = [0 for i in range (reference_size)]

        for qname in self.paf_file.lines_by_query:
            if len(self.paf_file.lines_by_query[qname]) > 1 :
                print(qname)
            for pline in self.paf_file.lines_by_query[qname]:
                try:
                    coverage_list[pline.target_start] += 1
                    coverage_list[pline.target_end if pline.target_end < reference_size else reference_size - 1] -= 1
                except Exception as e:
                    print(pline.target_start, pline.target_end)
                    print(e)
        
        res = []
        suma = 0
        for i in range (len(coverage_list)):
            suma += coverage_list[i]
            res.append(suma)
            res.append(0)

        x = []
        for i in range(len(res)//2):
            x.append(i)
            x.append(i + 0.5)
        
        fig = plt.figure()
        if ylim is not None:
            plt.ylim(top=ylim)
        plt.plot(x, res)
        fig.savefig(filename)
