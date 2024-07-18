import primer3 as p3
from types import SimpleNamespace
from General.constants import *
import multiprocessing
import argparse
from types import SimpleNamespace

max_tm=13.48

class PairFinder:

    def __init__(self, sequence_nt1, sequence_nt2=None):
        self.sequence_nt1 = sequence_nt1
        self.sequence_nt2 = sequence_nt2 if sequence_nt2 else sequence_nt1

    def forbidden_pairs(self, primer, p_end, sequence_nt, args):

        memo = {}

        def search(start, end):
            results = []

            # Base case
            if (end - start) <= 2 * (args.primer_lmax):

                tm = p3.bindings.calc_heterodimer(sequence_nt[start:end], primer, mv_conc=50.0, dv_conc=1.5,
                                                  dntp_conc=0.6,
                                                  dna_conc=50.0, temp_c=37.0, max_loop=30).tm
                if tm >= max_tm:
                    return handle_special_case(start, end)
                else:
                    return []

            mid = (start + end) // 2
            left_end = mid + (args.primer_lmax - 1)
            right_start = max(mid - (args.primer_lmax - 1), 0)

            first_tm = p3.bindings.calc_heterodimer(sequence_nt[start:left_end], primer, mv_conc=50.0, dv_conc=1.5,
                                                    dntp_conc=0.6,
                                                    dna_conc=50.0, temp_c=37.0, max_loop=30).tm

            second_tm = p3.bindings.calc_heterodimer(sequence_nt[right_start:end], primer, mv_conc=50.0, dv_conc=1.5,
                                                     dntp_conc=0.6,
                                                     dna_conc=50.0, temp_c=37.0, max_loop=30).tm

            # Check and recurse on the left segment
            if first_tm >= max_tm:
                results.extend(search(start, left_end))

            # Check and recurse on the right segment
            if second_tm >= max_tm:
                results.extend(search(right_start, end))

            return results

        def handle_special_case(start, end):
            if (start, end) in memo:
                return memo[(start, end)]

            # Avoid unnecessary recursion if the segment is already at the minimum length
            if end - start <= args.primer_lmax:
                tm = p3.bindings.calc_heterodimer(sequence_nt[start:end], primer, mv_conc=50.0, dv_conc=1.5,
                                                  dntp_conc=0.6,
                                                  dna_conc=50.0, temp_c=37.0, max_loop=30).tm
                if tm >= max_tm:
                    memo[(start, end)] = [(start, end)]
                    return [(start, end)]
                else:
                    memo[(start, end)] = []
                    return []

            # Recursive case: try shortening the segment from both ends
            shorten_left = handle_special_case(start + 1, end)
            shorten_right = handle_special_case(start, end - 1)

            result = list(set(shorten_left + shorten_right))

            memo[(start, end)] = result
            return result

        off_targets = []

        # search second half of sequence (with max allowed overlap)
        off_targets += search(max(p_end - args.allowed_overlap, 0), len(sequence_nt))

        return off_targets

    def find_all_pairs(self, args):

        forbidden_pairs = []
        # iterate over every start position
        for p_start in range(len(self.sequence_nt1) - args.primer_lmax + 1):
            p_end = p_start + args.primer_lmax
            primer = self.sequence_nt1[p_start:p_end]
            # single sequence logic
            if self.sequence_nt1 == self.sequence_nt2:
                off_targets = self.forbidden_pairs(primer, p_end, self.sequence_nt1, args)
            else:  # multiple sequence logic
                off_targets = self.forbidden_pairs(primer, 0, self.sequence_nt2, args)

            for target in off_targets:
                forbidden_pairs += [((p_start, p_end), target)]

        forbidden_pairs = set(forbidden_pairs)

        # process each forbidden pair pmax-mer in parallel to find forbidden pairs in p_lmin <= length <= p_lmax
        sub_pairs = self.sub_pairs_parallel(forbidden_pairs,args)

        return sub_pairs

    def process_pair(self, pair, args):
        # process each forbidden pair of pmax-mers to find forbidden pairs in p_lmin <= length <= p_lmax
        sub_forbidden_pairs = []

        primer1_start, primer1_end = pair[0]
        primer2_start, primer2_end = pair[1]

        primer1 = self.sequence_nt1[primer1_start:primer1_end]
        primer2 = self.sequence_nt2[primer2_start:primer2_end]

        mutreg_start = 0

        for length1 in range(args.primer_lmin, args.primer_lmax + 1):
            for length2 in range(args.primer_lmin, args.primer_lmax + 1):
                for start1 in range(args.primer_lmax - length1 + 1):
                    for start2 in range(args.primer_lmax - length2 + 1):
                        sub1 = primer1[start1:start1 + length1]
                        sub2 = primer2[start2:start2 + length2]
                        tm = p3.bindings.calc_heterodimer(sub1, sub2, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6,
                                                          dna_conc=50.0, temp_c=37.0, max_loop=30).tm
                        if tm >= max_tm:
                            sub_primer1= (primer1_start + start1 - mutreg_start, primer1_start +  start1 + length1 - mutreg_start)
                            sub_primer2 = (primer2_start + start2 - mutreg_start, primer2_start +start2 + length2 - mutreg_start)
                            sub_forbidden_pairs.append((sub_primer1,sub_primer2))


        return sub_forbidden_pairs

    def sub_pairs_parallel(self, forbidden_pairs,args):

        # Create a list of arguments to pass to process_pair
        tasks = [(pair, args) for pair in forbidden_pairs]

        # Process all forbidden pairs pmax_mers in parallel
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        result = pool.starmap(self.process_pair, tasks)
        pool.close()
        pool.join()

        # Flatten the list of lists
        sub_forbidden_pairs = set([item for sublist in result for item in sublist])

        return sub_forbidden_pairs




def subsequences(sequence, primer_lmin, primer_lmax):
    ls = []
    for j in range(primer_lmin, primer_lmax + 1):
        for i in range(len(sequence) - j + 1):
            start = i
            stop = i + j
            ls.append((sequence[start:stop], start, stop))
    return ls


def naive_brute_force(sequence1, sequence2, primer_lmin, primer_lmax, max_tm):
    primers1 = subsequences(sequence1, primer_lmin, primer_lmax)
    primers2 = subsequences(sequence2, primer_lmin, primer_lmax)

    high_risk_pairs = []

    for seq1, start1, stop1 in primers1:
        for seq2, start2, stop2 in primers2:
            tm = p3.bindings.calc_heterodimer(seq1, seq2, mv_conc=50.0, dv_conc=1.5,
                                              dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30).tm
            if tm >= max_tm:
                high_risk_pairs.append(((start1, stop1), (start2, stop2)))

    high_risk_pairs = set(high_risk_pairs)
    return high_risk_pairs

def main():

    parser = argparse.ArgumentParser(description='Process some sequences.')
    parser.add_argument('--primer_lmin', type=int, default=18, help='Minimum primer length')
    parser.add_argument('--primer_lmax', type=int, default=30, help='Maximum primer length')
    parser.add_argument('--allowed_overlap', type=int, default=6,
                        help='Maximum allowed overlap between selected primers')

    args = parser.parse_args()

    sequence1 = 'ATGGCGCTCCTGCTGTGCTTCGTGCTCCTGTGCGGAGTAGTGGATTTCGCCAGAA'
    sequence2 = 'ATGACATCGCGGAGATGGTTTCACCCAAATATCACTGGTGTGGAGGCAGAAAACC'

    naive_pairs = naive_brute_force(sequence1, sequence2, args.primer_lmin, args.primer_lmax, max_tm)

    print("Naive: ",naive_pairs)

    pair_finder= PairFinder(sequence1,sequence2)

    algorithm_pairs = pair_finder.find_all_pairs(args)


    print("Our algorithm: ",algorithm_pairs)

    if algorithm_pairs==naive_pairs:
        print("Equal!")



if __name__ == '__main__':
    main()