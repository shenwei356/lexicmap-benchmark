#!/usr/bin/env python

# by Antonio Camargo

import re
import sys


class Record:
    def __init__(self, fields: list[str]) -> None:
        self.qname: str = fields[0]
        self.tname: str = fields[2]
        self.pos: int = int(fields[3])
        self.cigar: str = fields[5]
        self.query_seq: str = fields[9]
        self.nm = self._get_nm_tag(fields)
        self.query_length = self._get_query_length()

    def _get_nm_tag(self, fields: list[str]) -> int:
        """Extracts the NM tag (edit distance) from SAM optional fields."""
        for field in fields[11:]:
            if field.startswith("NM:i:"):
                return int(field.split(":")[-1])
        return 0

    def _parse_cigar(self, cigar: str) -> tuple[int, int, int, int, int]:
        """
        Parses the CIGAR string and returns:
        - total_M: sum of M operations (alignment matches/mismatches)
        - total_I: sum of I operations (insertions in query)
        - total_D: sum of D operations (deletions from reference)
        - total_S: sum of S operations (soft clips)
        - total_H: sum of H operations (hard clips)
        - clip_start: number of clipped bases (soft or hard) at beginning
        """
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        total_M = total_I = total_D = total_S = total_H = 0
        clip_start = 0

        if ops:
            first_num, first_op = ops[0]
            if first_op in ("S", "H"):
                clip_start = int(first_num)

        for num, op in ops:
            num = int(num)
            if op == "M" or op == "=":
                total_M += num
            elif op == "I":
                total_I += num
            elif op == "D":
                total_D += num
            elif op == "S":
                total_S += num
            elif op == "H":
                total_H += num

        return total_M, total_I, total_D, total_S, total_H, clip_start

    def _get_query_length(self) -> int:
        """Determines the query length from SAM fields."""
        h_clips = re.findall(r"(\d+)H", self.cigar)
        hard_clipped = sum(int(h) for h in h_clips)
        return len(self.query_seq) + hard_clipped

    def process_alignment(self) -> tuple[str, str, int, int, int, int, int, int, float]:
        """Processes a SAM alignment record and returns relevant information."""
        if self.nm is None:
            return None

        total_M, total_I, total_D, total_S, total_H, clip_start = self._parse_cigar(self.cigar)
        query_start = clip_start + 1
        query_end = query_start + total_M + total_I - 1
        target_end = self.pos + total_M + total_D - 1
        alignment_length = total_M + total_I + total_D
        identity = (
            (alignment_length - self.nm) / alignment_length
            if alignment_length > 0
            else 0
        )

        return (
            self.qname,
            self.tname,
            self.query_length,
            query_start,
            query_end,
            self.pos,
            target_end,
            alignment_length,
            identity,
        )


def main() -> None:
    print(
        "query\ttarget\tquery_length\tquery_start\tquery_end\ttarget_start\t"
        "target_end\talignment_length\tidentity"
    )
    for line in sys.stdin:
        line = line.rstrip()

        if not line or line.startswith("@"):
            continue

        fields = line.split("\t")

        if len(fields) < 12:
            continue

        (
            qname,
            tname,
            query_length,
            query_start,
            query_end,
            target_start,
            target_end,
            alignment_length,
            identity,
        ) = Record(fields).process_alignment()

        print(
            f"{qname}\t{tname}\t{query_length}\t{query_start}\t{query_end}\t"
            f"{target_start}\t{target_end}\t{alignment_length}\t{identity:.4f}"
        )


if __name__ == "__main__":
    main()
