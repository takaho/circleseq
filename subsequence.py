import os, sys
import sqlite3

"""Subsequence retriever sscripts from genome.
Sequences were cached into SQLite database and reused if the same sequence was requested.
"""

__cacheofcache = {}

def get_subsequence(filename:str, seqname:str, start:int, stop:int=-1, step:int=4096) -> str:
    global __cacheofcache
    if filename in __cacheofcache:
        cache = __cacheofcache[filename]
    else:
        cache = __setup_dictionary(filename, step)
        __cacheofcache[filename] = cache
    if seqname is None and len(cache) > 0:
        seqname = list(cache.keys())[0]
    d = cache.get(seqname, {})
    # if the given chromosome is missing, return empty sequence
    if d is None or 'length' not in d or 'checkpoints' not in d:
        return ''
    length = d['length']
    checkpoints = d['checkpoints']
    if start < 0:
        start = 0
    if start >= length: # return empty sequence when the given position is out of- ange .
        return ''
    if stop < 0:
        stop = min(length, start + 100)
    left = 0
    right = len(checkpoints)
    while left < right:
        i = (left + right) // 2
        spos, fpos = checkpoints[i]
        if start < spos:
            right = i
        elif spos + step < start:
            right = i
        else:
            break
    center = (left + right) // 2
    spos, fpos = checkpoints[center]
    offset = start - spos
    required = stop - spos
    seq = ''
    with open(filename) as fi:
        fi.seek(fpos) # skip to cached position
        buffer = ''
        while 1:
            line = fi.readline()
            if line == '': break
            buffer += line.strip()
            if len(buffer) >= required:
                break
        return buffer[offset:offset + (stop - start)]

def __setup_dictionary(filename:str, step:int=65536)->dict:
    import time, json
    mtime:int = int(os.path.getmtime(filename))
    current:int = int(time.time())
    parent:str = os.path.dirname(filename)
    fn_cache:str = os.path.join(parent, '.subseq.db')
    # print(fn_cache)
    # exit()
    with sqlite3.Connection(fn_cache) as cnx:
        cur = cnx.cursor()
        cur.execute('create table if not exists fasta_cache (filename, num_seqs int, lastmodified int, step int, cache blob)')
        cur.execute('select filename, num_seqs, lastmodified, cache from fasta_cache where filename=? and step=?', 
            (os.path.basename(filename), step))
        r = cur.fetchone()
        if r is not None and r[2] >= mtime:
            data:dict = json.loads(r[3].decode('utf8'))
        else:
            cur.execute('delete from fasta_cache where filename=?', (os.path.basename(filename), ))
            data:dict = {}
            num_nuc:int = 0
            name:str = None
            checkpoints:list = []
            with open(filename) as fi:
                while 1:
                    fpos = fi.tell()
                    line = fi.readline()
                    if line == '': 
                        if name is not None and len(checkpoints) > 0 and num_nuc > 0:
                            checkpoints.append((num_nuc, fpos))
                            data[name] = {'length':num_nuc, 'checkpoints':checkpoints}
                        break
                    if line.startswith('>'):
                        if name is not None and len(checkpoints) > 0 and num_nuc > 0:
                            checkpoints.append((num_nuc, fpos))
                            data[name] = {'length':num_nuc, 'checkpoints':checkpoints}
                            sys.stderr.write('{}\t{}\t{}\n'.format(name, num_nuc, checkpoints[-1]))
                        name = line[1:].strip()
                        num_nuc = 0
                        checkpoints = [(0, fpos + len(line))]
                        next_step = step
                    else:
                        if num_nuc >= next_step:
                            checkpoints.append((num_nuc, fpos))
                            next_step += step
                        num_nuc += len(line.strip())
            cur.execute('insert into fasta_cache values(?, ?, ?, ?, ?)', 
                (os.path.basename(filename), len(data), current, step, json.dumps(data).encode('utf8')))
            cnx.commit()
    return data

def gen_random():
    import numpy as np
    for i in range(10):
        sys.stdout.write('>random_{}'.format(i))
        span = np.random.randint(2000, 10000)
        for j in range(span):
            if (j % 60 == 0): sys.stdout.write('\n')
            sys.stdout.write('ACGT'[np.random.randint(0, 4)])
        sys.stdout.write('\n')

if __name__ == '__main__':
    filename:str = sys.argv[1]
    pos:int = int(sys.argv[2]) if len(sys.argv) > 2 else 72933868
    data = __setup_dictionary(filename)
    for chrom in data.keys():
        print('{}\t{}\t{}'.format(chrom, pos, get_subsequence(filename, chrom, 72933868 -21, 72933868 + 20)))
