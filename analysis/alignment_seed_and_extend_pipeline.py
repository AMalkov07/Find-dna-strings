from collections import defaultdict
from typing import List, Dict, Tuple

from utils.data_structures import SeedExtendSettings, SeedExtendCluster

class SeedAndExtend:
    def __init__(self, mutagenic_zone: str, pattern: str, seed_extend_settings: SeedExtendSettings):
        self.mutagenic_zone = mutagenic_zone
        self.pattern = pattern
        self.seed_extend_settings = seed_extend_settings

    def _build_kmer_index(self) -> Dict[str, List[int]]:
        mutagenic_zone = self.mutagenic_zone
        n_mutagenic_zone = len(mutagenic_zone)
        k = self.seed_extend_settings.k
        idx = defaultdict(list)

        for i in range(n_mutagenic_zone - k + 1):
            kmer = mutagenic_zone[i:i+k]
            idx[kmer].append(i)
        return idx

    def _find_seed_hits(self, kmer_indexes: Dict[str, List[int]]) -> List[Tuple[int, int]]:
        hits: List[Tuple[int, int]] = []
        mutagenic_zone = self.mutagenic_zone
        k = self.seed_extend_settings.k
        
        for qpos in range(len(mutagenic_zone) - k + 1):
            kmer = mutagenic_zone[qpos:qpos + k]
            for rpos in kmer_indexes.get(kmer, ()):
                hits.append((qpos, rpos))
        return hits

    def _cluster_hits_by_offset(self, hits: List[Tuple[int, int]]) -> List[SeedExtendCluster]:
        offset_tolerance = self.seed_extend_settings.offset_tolerance
        min_seeds = self.seed_extend_settings.min_seed

        offsets = [r - q for (q, r) in hits]
        # group offsets by rounding/nearby bucket
        # Greedy clustering by offset value
        offsets_sorted = sorted(set(offsets))
        clusters = []
        for off in offsets_sorted:
            placed = False
            for cl in clusters:
                #cluster_mean = statistics.median(cl["offsets"]) #<<<<<<<<note: figure out if the mean is actually the metric that should be used to measure distance from
                if abs(cl["offsets"][0] - off) <= offset_tolerance: 
                    cl['offsets'].append(off)
                    placed = True
                    break
            if not placed:
                clusters.append({'offsets': [off]})
        # For each cluster, gather supporting hits
        clusters_out: List[SeedExtendCluster] = []
        for cl in clusters:
            off_vals = set(cl['offsets'])
            supporting = [(q, r) for (q, r) in hits if (r - q) in off_vals]
            if len(supporting) < min_seeds: #min_seeds determines the minimum number of kmer's that have to be close enough to eachother for us to consider the cluster. To my understanding, if the kmer size is 4 and min_seed is 2, then a single substring match of 5 bp would qualify as a valid cluster
                continue
            q_positions = [q for q, _ in supporting]
            r_positions = [r for _, r in supporting]
            clusters_out.append(SeedExtendCluster(
                offset=int(sum(off_vals) / len(off_vals)),
                hits=supporting,
                qmin=min(q_positions),
                qmax=max(q_positions),
                rmin=min(r_positions),
                rmax=max(r_positions),
                n_seeds=len(supporting)
            ))
        return clusters_out

    def _extend_cluster(self, cluster: SeedExtendCluster):
        pass


    def execute(self):
        kmer_indexes: Dict[str, List[int]] = self._build_kmer_index()
        hits: List[Tuple[int, int]] = self._find_seed_hits(kmer_indexes)
        if not hits:
            return []

        clusters: List[SeedExtendCluster] = self._cluster_hits_by_offset(hits)

        if not clusters:
            return []

        results = []
        all_used_regions = set()
        for cl in clusters:
            alns, used_regions = self._extend_cluster(cl)

