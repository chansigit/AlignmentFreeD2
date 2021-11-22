def build_kmer(sequence, k):
    kmers = []
    n_kmers = len(sequence) - k + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + k]
        kmers.append(kmer)
    return kmers

def enumerateKmer(k):
    assert k<=13 # k=13 requires 16s
    import itertools
    x = ['A','G','C','T']
    return ["".join(p) for p in itertools.product(x, repeat=k)]


def get_letter_prob(kmerList):
    longtext = "".join(kmerList)
    nBase = len(longtext)
    nA= (longtext.count("A")+longtext.count("a") )/nBase
    nT= (longtext.count("T")+longtext.count("t") )/nBase
    nG= (longtext.count("G")+longtext.count("g") )/nBase
    nC= (longtext.count("C")+longtext.count("c") )/nBase
    nU= (longtext.count("U")+longtext.count("u") )/nBase
    nN= (longtext.count("N")+longtext.count("n") )/nBase
    assert (nT==0 and nU!=0) or (nT!=0 and nU==0)
    if nT ==0:
        nA += nN/4
        nU += nN/4
        nG += nN/4
        nC += nN/4
        return {'A':nA,'U':nU,'G':nG,'C':nC}
    if nU == 0:
        nA += nN/4
        nT += nN/4
        nG += nN/4
        nC += nN/4
        return {'A':nA,'T':nT,'G':nG,'C':nC}

def compute_prob(w, letterProb):
    prob=1.0
    for s in w:
        prob = prob*letterProb[s]
    return prob

def enumerate_kmer_prob(kmerList, letterProb):
    kmerProb = dict()
    for kmer in kmerList:
        kmerProb[kmer] = compute_prob(kmer, letterProb)
    return kmerProb


def D2a(seqX, seqY, k=5, kmerEnumerated=None, kmerProb=None):
    from collections import Counter
    import math
    kmerCntX = Counter(build_kmer(seqX, k))
    kmerCntY = Counter(build_kmer(seqY, k))
    
    if kmerEnumerated is None:
        pass
    if kmerProb is None:
        # calculate kmerProb from seqX and seqY
        pass
    
    result=0
    adj_n =len(seqX)-k+1
    adj_m =len(seqY)-k+1
    for kmer in kmerEnumerated:
        prob = kmerProb[kmer]
        adj_kmerCntX = kmerCntX[kmer] - adj_n * prob
        adj_kmerCntY = kmerCntY[kmer] - adj_m * prob
        numerator    = adj_kmerCntX*adj_kmerCntY
        denominator  = math.sqrt(adj_n*adj_m) * prob
        quotient     = numerator/denominator if denominator!=0 else 0
        result+= quotient
    return result

def D2S(seqX, seqY, k=5, kmerEnumerated=None, kmerProb=None):
    from collections import Counter
    import math
    kmerCntX = Counter(build_kmer(seqX, k))
    kmerCntY = Counter(build_kmer(seqY, k))
    
    if kmerEnumerated is None:
        pass
    if kmerProb is None:
        # calculate kmerProb from seqX and seqY
        pass
    
    result=0
    adj_n =len(seqX)-k+1
    adj_m =len(seqY)-k+1
    for kmer in kmerEnumerated:
        prob = kmerProb[kmer]
        adj_kmerCntX = kmerCntX[kmer] - adj_n * prob
        adj_kmerCntY = kmerCntY[kmer] - adj_m * prob
        numerator    = adj_kmerCntX*adj_kmerCntY
        denominator  = math.sqrt(adj_kmerCntX*adj_kmerCntX + adj_kmerCntY*adj_kmerCntY) 
        quotient     = numerator/denominator if denominator!=0 else 0
        result+= quotient
    return result
  
  
  
#similarities = []
#for r in tqdm(read_unmapped):
#    sim= D2a(r.query_sequence, ITR1, k=5, kmerEnumerated= kmerEnumerated, kmerProb=kmerProb)
#    similarities.append(sim)
    
from joblib import Parallel, delayed
similarities_d2a= Parallel(n_jobs=10)(
    delayed(D2a)(
        r.query_sequence, ITR1, 
        k=5, kmerEnumerated= kmerEnumerated, kmerProb=kmerProb
    ) for r in read_unmapped
)

import matplotlib.pyplot as plt
_ = plt.hist(similarities_d2a,bins=100, log=False)


#similarities2 = []
#for r in tqdm(read_unmapped):
#    sim= D2S(r.query_sequence, ITR1, k=5, kmerProb=prob_model)
#    similarities2.append(sim)
    
    
from joblib import Parallel, delayed
similarities_d2S= Parallel(n_jobs=10)(
    delayed(D2S)(
        r.query_sequence, ITR1, 
        k=5, kmerEnumerated= kmerEnumerated, kmerProb=kmerProb
    ) for r in read_unmapped
)

import matplotlib.pyplot as plt
_ = plt.hist(similarities_d2S,bins=100, log=False)
