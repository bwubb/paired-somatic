

def get_principleTX():
    pTX={}
    with open('/home/bwubb/resources/preferred_transcripts.20170505.txt','r') as file:
        for line in file:
            try:
                k,v=line.rstrip().split('\t')
                pTX[k]=v
            except IndexError:
                pass
    return pTX

def fix_blanks(x):
    for i,j in enumerate(x):
        if j=='.':
            x[i]=''
    return x

def _exonic_genes(fields,eg,aa,pTX):
    #eg=row_dict['Gene_refGene'].split(',')
    #aa=row_dict['AAChange_refGene'].split(',')
    g,a,n,e,t=([] for x in range(5))
    for i in eg:
        V=None
        for j in aa:
            if pTX[i] in j:
                V=j.split(':')
        if V:
            assert len(V)>3,'AssertionError: _exonic_genes - gene_info = {0}'.format(V)#Nonframeshift substitutions do not have amino annotation
            try:
                a.append(V[4])
            except IndexError:
                a.append('')
            n.append(V[3])
            e.append(V[2][4:])
            t.append(V[1])
            g.append(V[0])
    if len(g)>0:
        for k,v in zip(['Gene','Transcript','Exon','NTChange','AAChange'],[g,t,e,n,a]):
            fields[k]+=[','.join(v)]
        fields['GenomicRegion']+=['exonic']
    return fields

def _splicing_genes(fields,sp,gd,pTX):
        g,n,e,t=([] for x in range(4))
        for i in sp:
            V=None
            for j in gd:
                if pTX[i] in j:
                    V=j.split(':')
            if V:
                assert len(V)==3,'AssertionError: _splicing_genes - gene_info = {0}'.format(V)
                n.append(V[2])
                e.append(V[1])
                t.append(V[0])
                g.append(i)
        if len(g)>0:
            for k,v in zip(['Gene','Transcript','Exon','NTChange'],[g,t,e,n]):
                fields[k]+=[','.join(v)]
            fields['GenomicRegion']+=['splicing']
        return fields

def _UTR_genes(self,fields,ug,gd,utr,pTX):
    g,n,t=([] for x in range(3))
    for i in ug:
        V=None
        for j in gd:
            if pTX[i] in j:
                V=j.split(':')
        if V:
            assert len(V)==2,'AssertionError: _UTR_genes - gene_info = {0}'.format(V)
            n.append(V[1])
            t.append(V[0])
            g.append(i)
    if len(g)>0:
        for k,v in zip(['Gene','Transcript','NTChange'],[g,t,n]):
            fields[k]+=[','.join(v)]
        fields['GenomicRegion']+=[utr]
    return fields

def _intergenic_genes(self,fields,ig,gd,pTX):
    d,g,t=([] for x in range(3))
    for n,i in enumerate(ig):
        if pTX[i]:
            g.append(i)
            d.append(':'.join([i,gd[n]]))
            t.append(pTX[i])
    if len(g)>0:
        for k,v in zip(['IntergenicDist','Gene','Transcript'],[d,g,t]):
            fields[k]+=[','.join(v)]
        fields['GenomicRegion']+=['intergenic']
    return fields

def _other_genes(self,fields,og,func,pTX):
    g,t=([] for x in range(2))
    for i in og:
        if pTX.get(i,False):#####Look at this notation! Otherwise you can get a key error from shit like ASB16-AS1;ASB16-AS1
            g.append(i)
            t.append(pTX[i])
    if len(g)>0:
        for k,v in zip(['Gene','Transcript'],[g,t]):
            fields[k]+=[','.join(v)]
        fields['GenomicRegion']+=[func]
    return fields

