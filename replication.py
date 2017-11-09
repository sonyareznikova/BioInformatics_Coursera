import urllib2
import matplotlib.pyplot as plt

def FrequentWords(Text, k):
    FrequentPatterns = []
    count = CountDict(Text, k)
    m = max(count.values())
    for c in count:
        if count[c] == m:
            FrequentPatterns.append(Text[c:c+k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates

def FrequentWordsOver(Text, k, n):
    FrequentPatterns = {}
    count = CountDict(Text, k)
    for c in count:
        if count[c] >= n:
            FrequentPatterns[Text[c:c+k]] = count[c];
    FrequentPatternsNoDuplicates = remove_dups_over(FrequentPatterns)
    return FrequentPatterns

def remove_dups_over(items):
    no_dups = {}
    for key,value in items.items():
        if key not in no_dups.keys():
            no_dups[key] = value
    return no_dups

def remove_duplicates(Items):
    ItemsNoDuplicates = [] # output variable
    for i in Items:
        if i not in ItemsNoDuplicates:
            ItemsNoDuplicates.append(i)
    return ItemsNoDuplicates

def CountDict(Text, k):
    Count = {} # output variable
    for i in range(len(Text)-k+1):
        c = PatternCount(Text[i:i+k], Text)
        Count[i] = c
    return Count

def PatternCount(Pattern, Text):
    count = 0 # output variable
    for i in range(len(Text)-1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(len(Text)-1):
        if Text[i:i+len(Pattern)] == Pattern or hammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count += 1
    return count

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def ReverseComplement(Pattern):
    rev = reverse(Pattern)
    revComp = complement(rev)
    return revComp

def reverse(string):
    list = []
    for s in string:
        list.append(s)
    list.reverse()
    return "".join(list)

def complement(Nucleotide):
    comp = ''  # output variable
    for n in Nucleotide:
        if n == 'A': comp += 'T'
        elif n == 'T': comp += 'A'
        elif n == 'C': comp += 'G'
        elif n == 'G': comp += 'C'
    return comp

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range(len(Genome)-1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-1):
        if Text[i:i+len(Pattern)] == Pattern or hammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

def hammingDistance(p, q):
    hd = 0
    for i in range(len(p)):
        if p[i] != q[i]: hd += 1
    if len(p) < len(q): hd += abs(len(p) - len(q))
    return hd

def Skew(Genome):
    skew = {} #initializing the dictionary
    skew[0] = 0
    for i in range(1, len(Genome)+1):
        if Genome[i-1] == "G": skew[i] = skew[i-1]+1
        elif Genome[i-1] == "C": skew[i] = skew[i-1]-1
        else: skew[i] = skew[i-1]
    return skew

def MinimumSkew(Genome):
    positions = [] # output variable
    sk = Skew(Genome)
    m = min(sk.values())
    for key,value in sk.items():
        if value == m: positions.append(key)
    return positions

#returns a dictionary with 9-mers as keys and the sum of their and their reverse complements count as values - a list of potential dnaA boxes
def findDNAaBox(ori):
    initList = FrequentWordsOver(ori, 9, 3)
    for i in initList:
        rev = ReverseComplement(i)
        count = PatternCount(rev, ori)
        initList[i] += count
    return initList



# data = urllib2.urlopen("http://bioinformaticsalgorithms.com/data/realdatasets/Replication/v_cholerae_oric.txt")
# data = urllib2.urlopen("http://bioinformaticsalgorithms.com/data/realdatasets/Replication/t_petrophila_oriC.txt")
# data = urllib2.urlopen("http://bioinformaticsalgorithms.com/data/realdatasets/Replication/E_coli.txt")

# f = open('e.coli.txt')

# plt.plot(sk.items())
# plt.ylabel('E.coli skew')
# plt.show()


