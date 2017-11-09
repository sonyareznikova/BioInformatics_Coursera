import random
import urllib2

def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = [0]*k
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def CountWithPseudocounts(Motifs):
    count = {}  # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = [1] * k
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)  # number of rows
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs)
    for key in profile.keys():
        for j in range(k):
            profile[key][j] = float(profile[key][j]) / (t+4)
    return profile

def Profile(Motifs):
    t = len(Motifs) #number of rows
    k = len(Motifs[0])
    profile = Count(Motifs)
    for key in profile.keys():
        for j in range(k):
            p = profile[key][j]
            profile[key][j] = float(profile[key][j]) / t
    return profile

def ConsensusWithPseudocounts(Motifs):
    profile = ProfileWithPseudocounts(Motifs)
    consensus = ""
    for j in range(len(Motifs[0])):
        m = 0.0
        maxKey = ""
        for key in profile.keys():
            if profile[key][j] > m:
                m = profile[key][j]
                maxKey = key
        consensus += maxKey
    return consensus

def Consensus(Motifs):
    profile = Profile(Motifs)
    consensus = ""
    for j in range(len(Motifs[0])):
        m = 0.0
        maxKey = ""
        for key in profile.keys():
            if profile[key][j] > m:
                m = profile[key][j]
                maxKey = key
        consensus += maxKey
    return consensus

def ScoreWithPseudocounts(Motifs):
    score = 0
    consensus = ConsensusWithPseudocounts(Motifs)
    for j in range(len(Motifs[0])):
        for i in range(len(Motifs)):
            if Motifs[i][j] != consensus[j]: score += 1
    return score

def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    for j in range(len(Motifs[0])):
        for i in range(len(Motifs)):
            if Motifs[i][j] != consensus[j]: score += 1
    return score

def Pr(Text, Profile):
    firstLetter = Text[0]
    final_pr = 1
    for j in range(len(Text)):
        for key in Profile.keys():
            if key == Text[j]: final_pr *= Profile[key][j]
    return final_pr

def ProfileMostProbablePattern(Text, k, Profile):
    maxProbability = -1
    final_kMer = ""
    for i in range(len(Text)-k+1):
        someProbability = Pr(Text[i:i+k], Profile)
        if someProbability > maxProbability:
            maxProbability = someProbability
            final_kMer = Text[i:i+k]
    return final_kMer

def GreedyMotifSearch(Dna, k, t):
    bestMotifs = []
    for i in range(0, t):
        bestMotifs.append(Dna[i][0:k])
    for i in range(len(Dna[0])-k+1):
        motifs = []
        motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(motifs[0:j])
            motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(motifs) < Score(bestMotifs): bestMotifs = motifs
    return bestMotifs

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    bestMotifs = []
    for i in range(0, t):
        bestMotifs.append(Dna[i][0:k])
    for i in range(len(Dna[0]) - k + 1):
        motifs = []
        motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(motifs[0:j])
            motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if ScoreWithPseudocounts(motifs) < ScoreWithPseudocounts(bestMotifs): bestMotifs = motifs
    return bestMotifs

def Motifs(Profile, Dna):
    k = len(Profile['A'])
    bestMotifs = []
    for i in range(len(Dna)):
        bestMotifs.append(ProfileMostProbablePattern(Dna[i], k, Profile))
    return bestMotifs

def RandomMotifs(Dna, k, t):
    randomMotifs = []
    for i in range(t):
        startNumber = random.randint(1, len(Dna[i])-k)
        randomMotifs.append(Dna[i][startNumber:startNumber+k])
    return randomMotifs

def RandomizedMotifSearch(Dna, k, t):
    bestMotifs = RandomMotifs(Dna, k, t)
    while True:
        Profile = ProfileWithPseudocounts(bestMotifs)
        newMotifs = Motifs(Profile, Dna)
        if ScoreWithPseudocounts(newMotifs) < ScoreWithPseudocounts(bestMotifs):
            bestMotifs = newMotifs
        else:
            return bestMotifs

def RandomizedMotifSearchV1(Dna, motifs):
    bestMotifs = motifs
    Profile = ProfileWithPseudocounts(bestMotifs)
    newMotifs = Motifs(Profile, Dna)
    return newMotifs

def Normalize(Probabilities):
    sum = 0
    for value in Probabilities.values():
        sum += value
    for key in Probabilities.keys():
        Probabilities[key] = Probabilities[key]/sum
    return Probabilities

def NormalizeList(items):
    number = sum(items)
    normalized = list(map(lambda x: x/number, items))
    print normalized


def WeightedDie(Probabilities):
    kMer = ''
    someNumber = random.uniform(0,1)
    sum = 0
    for key in Probabilities.keys():
        if someNumber > sum and someNumber < sum + Probabilities[key]:
            kMer = key
        sum += Probabilities[key]
    return kMer

def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    motifs = RandomMotifs(Dna, k, t)
    bestMotifs = motifs
    for j in range(1,N):
        i = random.randint(1,t-1)
        motifs.remove(motifs[i])
        profile = ProfileWithPseudocounts(motifs)
        someMotif = ProfileGeneratedString(Dna[i], profile, k)
        motifs.insert(i, someMotif)
        if ScoreWithPseudocounts(motifs) < ScoreWithPseudocounts(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

dna = ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"]
#dnaRaw = urllib2.urlopen("http://bioinformaticsalgorithms.com/data/challengedatasets/DosR.txt")
#dna = []
#for line in dnaRaw:
#    dna.append(line[:len(line)-2])
NormalizeList([0.45, 0.63, 0.09, 0.27, 0.36])

motifs = ['TGA', 'GTT', 'GAA', 'TGT']
bm = RandomizedMotifSearchV1(dna, motifs)
print bm

