#Esta función devuelve la cantidad de cada nucleótido
def count_nt(seq):
    print('Cant A: ' + str(seq.upper().count('A')) + ', Cant C: ' + str(seq.upper().count('C')) + ', Cant G: ' + str(seq.upper().count('G')) + ', Cant T: ' + str(seq.upper().count('T')))

    
#Esta función devuelve la cantidad del/de los aminoácido/s indicado/s en la secuencia de proteína (ingresar los aminoácidos deseados como una lista)
def porc_aa(seq, aa_list):
    total = 0
    
    for aa in aa_list:
        aa = aa.upper()
        print('La proteína es {:.2f}% {}'.format(seq.upper().count(aa.upper()), aa.upper()))
        total += seq.upper().count(aa.upper())

    print('La proteína es {:.2f}% {}'.format(total, ', '.join(aa_list)))


#Transcripción ADN a ARN
def transcripcion(seq):
    arn = seq.upper().replace('T', 'U')
    print('Secuencia de ARN: ' + seq.upper().replace('T', 'U'))
    return arn


#Traducción ARN a proteína
def traduccion(seq):
    seq = seq.upper()
    code = {'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T', 'AAC':'N', 'AAU':'N', 
            'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R', 'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P', 'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R',  
            'CGG':'R', 'CGU':'R', 'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 
            'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G', 'UCA':'S', 'UCC':'S', 
            'UCG':'S', 'UCU':'S', 'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L', 'UAC':'Y', 'UAU':'Y', 'UGC':'C', 'UGU':'C', 
            'UGG':'W'}

    stop = ['UAA', 'UAG', 'UGA']

    protein = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in stop:
            break
        else:
            protein += code[codon]

    print('Secuencia de proteína: ' + protein)  
    return protein


#Cálculo del peso molecular de una proteína
def pm_prot(seq):
    pm_dic = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841, 'G': 57.02146, 'H': 137.05891,
              'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858,
              'R': 156.10111, 'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333 }

    pm = 0
    for aa in seq.upper():
        pm += pm_dic[aa]
    
    print('El peso molecular de la proteína ingresada es ' + str(pm))
    return pm


#Esta función devuelve la hebra reversa complementaria
def rev_compl(seq):
    seq_rc = seq.upper()[::-1]
               
    for nucleotide in seq_rc:    
        if nucleotide == 'A':
            seq_rc = seq_rc.replace('A', 'T')
        elif nucleotide == 'T':
            seq_rc = seq_rc.replace('T', 'A')
        elif nucleotide == 'C':
            seq_rc = seq_rc.replace('C', 'G')
        elif nucleotide == 'G':
            seq_rc = seq_rc.replace('G', 'C')
             
    print('Hebra reversa complementaria: ' + seq_rc)
    return seq_rc

#Esta función devuelve el número de mutaciones puntuales entre dos secuencias del mismo largo
def mut_punt(seq1, seq2):
    if len(seq1) != len(seq2):
        print('Las secuencias no presentan el mismo largo.')
    else:
        q = 0
        n = 0
        for nucleotide in seq1:
            if n < len(seq1):
                if seq1.upper()[n] != seq2.upper()[n]:
                    q += 1
                    n += 1
                else:
                    n += 1

        print('Entre las dos secuencias hay ' + str(q) + ' mutaciones puntuales') 
        return q


#Calcular contenido GC
def cont_gc(seq):
    contenidogc = (seq.count('G') + seq.count('g') + seq.count('C') + seq.count('c')) *100 / len(seq)
    print('El contenido de GC es del ' + str(contenidogc) + '%')        
    return contenidogc


#Esta función devuelve el patrón más frecuente de largo k ingresado
def pat_freq(seq, k):
    freq = {}
    n = len(seq)
    for i in range(n-k+1):
        patron = seq[i:i+k]
        freq[patron] = 0
        for i in range(n-k+1):
            if seq[i:i+k] == patron:
                freq[patron] += 1

    values = list(freq.values())
    max_freq = max(values)
    most_freq = []

    for key in freq:
        if freq[key] == max_freq:
            most_freq.append(key)
        
    print ('El patrón más frecuente es: ' + str(most_freq))
    return most_freq


#Esta función devuelve las posiciones donde se encuentra un patrón considerando la cantidad de diferencias indicada como dif
def pos_patron(seq, pat, dif):
    seq = seq.upper()
    pat = pat.upper()

    positions = [] 

    for i in range(len(seq)-len(pat)+1):
        q = 0
        subtext = seq[i:i+len(pat)]
        for j in range(len(pat)):
            if subtext[j] != pat[j]:
                q += 1
        if q <= dif:
            positions.append(i)
   
    print ('El patrón se encuentra en las posiciones: ' + str(positions) + ' considerando como máximo ' + str(dif) + ' diferencias')
    return positions


#Esta función devuelve la cantidad de veces que aparece un patrón considerando la cantidad de diferencias indicada como dif
def cant_patron(seq, pat, dif):
    seq = seq.upper()
    pat = pat.upper()

    cant = 0 

    for i in range(len(seq)-len(pat)+1):
        q = 0
        subtext = seq[i:i+len(pat)]
        for j in range(len(pat)):
            if subtext[j] != pat[j]:
                q += 1
        if q <= dif:
            cant += 1
   
    print ('El patrón se encuentra ' + str(cant) + ' veces considerando como máximo ' + str(dif) + ' diferencias')
    return cant


#Esta función devuelve la secuencia consenso a partir de una lista de secuencias cargadas
def seq_consenso(s):
    count = {}
    k = len(s[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(s)
    for i in range(t):
        for j in range(k):
            symbol = s[i].upper()[j]
            count[symbol][j] += 1 

    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    print('La secuencia consenso es ' + str(consensus))
    return consensus


#Esta función devuelve una matriz de probabilidad para cada nucleótido a partir de una lista de secuencias cargadas
def matriz_prob(s):
    count = {}
    k = len(s[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(s)
    for i in range(t):
        for j in range(k):
            symbol = s[i].upper()[j]
            count[symbol][j] += 1 

    profile = {}
    for symbol in "ACGT":
        profile[symbol] = []
        for element in count[symbol]:
            freq = round(element / t, 2)
            profile[symbol].append(freq)

    print('La matriz de probabilidad es ' + str(profile))
    return profile
