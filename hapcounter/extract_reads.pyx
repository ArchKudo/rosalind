import pysam
import polars as pl
cimport cython

# CIGAR op codes
cdef int OP_M = 0      # M
cdef int OP_I = 1      # I
cdef int OP_D = 2      # D
cdef int OP_N = 3      # N (ref skip)
cdef int OP_S = 4      # S
cdef int OP_H = 5      # H
cdef int OP_P = 6      # P
cdef int OP_EQ = 7     # =
cdef int OP_X = 8      # X

@cython.cfunc
@cython.inline
cdef int _qpos_for_refpos(object cigartuples, int ref_start, int target_ref_pos) noexcept:
    """
    Return query index for target_ref_pos, or:
      -1 if target falls in deletion/refskip
      -2 if not covered by alignment
    """
    cdef int ref_pos = ref_start
    cdef int qpos = 0
    cdef int op, length

    if cigartuples is None:
        return -2

    for ct in cigartuples:
        op = ct[0]
        length = ct[1]

        if op == OP_M or op == OP_EQ or op == OP_X:
            if target_ref_pos >= ref_pos and target_ref_pos < ref_pos + length:
                return qpos + (target_ref_pos - ref_pos)
            ref_pos += length
            qpos += length
        elif op == OP_I:
            qpos += length
        elif op == OP_D or op == OP_N:
            if target_ref_pos >= ref_pos and target_ref_pos < ref_pos + length:
                return -1  # deletion or ref skip at this position
            ref_pos += length
        elif op == OP_S:
            qpos += length
        elif op == OP_H or op == OP_P:
            pass
        else:
            pass

    return -2

@cython.cfunc
@cython.inline
cdef int _ref_end_from_cigar(object cigartuples, int ref_start) noexcept:
    """Compute reference_end (one past last ref base) from CIGAR."""
    if cigartuples is None:
        return ref_start
    cdef int ref_pos = ref_start
    cdef int op, length
    for ct in cigartuples:
        op = ct[0]
        length = ct[1]
        if op == OP_M or op == OP_EQ or op == OP_X or op == OP_D or op == OP_N:
            ref_pos += length
        else:
            # I/S/H/P do not consume reference
            pass
    return ref_pos

@cython.cfunc
@cython.inline
cdef int _find_first_ge_list(list arr, long x) noexcept:
    """Binary search on a Python list of ints: first index i where arr[i] >= x; returns len(arr) if none."""
    cdef Py_ssize_t lo = 0
    cdef Py_ssize_t hi = len(arr)
    cdef Py_ssize_t mid
    cdef long v
    while lo < hi:
        mid = (lo + hi) >> 1
        v = <long>arr[mid]
        if v < x:
            lo = mid + 1
        else:
            hi = mid
    return <int>lo


# Disable bounds check and negative indexing for lists
@cython.boundscheck(False)
@cython.wraparound(False)
def extract_read_info(str bam_path, list variant_records, int window_bp=50000):
    """
    window_bp: max span (bp) for grouping nearby positions into a single fetch.
    """
    
    
    cdef object bam = pysam.AlignmentFile(bam_path)

    # Group variants by (chrom, pos)
    cdef dict posmap = {}   # {chrom: {pos: [(ref, alt), ...]}}
    cdef object chrom_key
    cdef int pos_i
    cdef str ref_s, alt_s

    for chrom_key, pos_i, ref_s, alt_s in variant_records:
        if chrom_key not in posmap:
            posmap[chrom_key] = {}
        if pos_i not in posmap[chrom_key]:
            posmap[chrom_key][pos_i] = []
        posmap[chrom_key][pos_i].append((ref_s, alt_s))

    # Column buffers
    cdef list chroms = []
    cdef list poss = []
    cdef list refs = []
    cdef list alts = []
    cdef list haps = []
    cdef list bases = []

    # Declarations for loop variables
    cdef object aln
    cdef int qpos, hp, target0
    cdef str base
    cdef object cigartuples
    cdef int ref_start, ref_end
    cdef object seq_obj
    cdef bytes seq_bytes
    cdef unsigned char b
    cdef list pos_list_py
    cdef list win_positions_py
    cdef dict chrom_posmap
    cdef list var_list
    cdef int i, j, k, n, pos_len
    cdef int win_start_pos, win_end_pos
    cdef int fetch_start0, fetch_end0, seq_len
    cdef long read_ref_start_pos1, read_ref_end_pos1, p1

    # Iterate chromosomes
    for chrom_key in posmap:
        chrom_posmap = posmap[chrom_key]
        pos_list_py = sorted(chrom_posmap.keys())
        n = len(pos_list_py)
        i = 0
        while i < n:
            win_start_pos = <int>pos_list_py[i]
            win_end_pos = win_start_pos
            j = i + 1
            while j < n and (<long>pos_list_py[j]) - win_start_pos <= window_bp:
                win_end_pos = <int>pos_list_py[j]
                j += 1

            win_positions_py = pos_list_py[i:j]

            # Fetch reads once for the whole window (0-based half-open)
            fetch_start0 = win_start_pos - 1
            fetch_end0 = win_end_pos  # half-open end

            for aln in bam.fetch(chrom_key, fetch_start0, fetch_end0):
                # Fast filters
                if getattr(aln, "is_unmapped", False) \
                   or getattr(aln, "is_secondary", False) \
                   or getattr(aln, "is_supplementary", False) \
                   or getattr(aln, "is_qcfail", False) \
                   or getattr(aln, "is_duplicate", False):
                    continue

                # Try to get HP; skip if absent
                try:
                    hp = aln.get_tag("HP")
                except Exception:
                    continue

                # Get essentials once
                try:
                    cigartuples = aln.cigartuples
                    ref_start = aln.reference_start
                except Exception:
                    cigartuples = None
                    ref_start = getattr(aln, "reference_start", -1)

                if ref_start < 0:
                    continue

                ref_end = _ref_end_from_cigar(cigartuples, ref_start)

                # Window prune: only positions overlapping read ref span
                read_ref_start_pos1 = ref_start + 1
                read_ref_end_pos1 = ref_end  # 1-based inclusive end

                pos_len = len(win_positions_py)
                k = _find_first_ge_list(win_positions_py, read_ref_start_pos1)
                if k >= pos_len:
                    continue

                # Prepare sequence once per read
                seq_obj = aln.query_sequence
                if seq_obj is None:
                    continue
                if isinstance(seq_obj, bytes):
                    seq_bytes = <bytes>seq_obj
                else:
                    try:
                        seq_bytes = (<str>seq_obj).encode("ascii")
                    except Exception:
                        try:
                            seq_bytes = bytes(seq_obj)
                        except Exception:
                            continue
                seq_len = len(seq_bytes)

                # Iterate positions in window that are within read span
                while k < pos_len:
                    p1 = <long>win_positions_py[k]
                    if p1 > read_ref_end_pos1:
                        break

                    target0 = <int>(p1 - 1)

                    if cigartuples is not None:
                        qpos = _qpos_for_refpos(cigartuples, ref_start, target0)
                        if qpos >= 0 and qpos < seq_len:
                            b = seq_bytes[qpos]
                            base = chr(b)
                            var_list = chrom_posmap[<int>p1]
                            for ref_s, alt_s in var_list:
                                chroms.append(chrom_key)
                                poss.append(<int>p1)
                                refs.append(ref_s)
                                alts.append(alt_s)
                                haps.append(hp)
                                bases.append(base)
                    else:
                        # Very slow fallback for mocks lacking cigartuples
                        try:
                            ref_positions = aln.get_reference_positions(full_length=True)
                            if 0 <= target0 < len(ref_positions):
                                qpos = -1
                                for ii in range(len(ref_positions)):
                                    if ref_positions[ii] == target0:
                                        qpos = ii
                                        break
                                if qpos >= 0 and qpos < seq_len:
                                    b = seq_bytes[qpos]
                                    base = chr(b)
                                    var_list = chrom_posmap[<int>p1]
                                    for ref_s, alt_s in var_list:
                                        chroms.append(chrom_key)
                                        poss.append(<int>p1)
                                        refs.append(ref_s)
                                        alts.append(alt_s)
                                        haps.append(hp)
                                        bases.append(base)
                        except Exception:
                            pass

                    k += 1

            # Advance to next window
            i = j

    try:
        bam.close()
    except Exception:
        pass

    if chroms:
        return pl.DataFrame(
            {
                "chrom": chroms,
                "pos": poss,
                "ref": refs,
                "alt": alts,
                "hap": haps,
                "base": bases,
            }
        )
    else:
        return pl.DataFrame(schema=["chrom", "pos", "ref", "alt", "hap", "base"])
