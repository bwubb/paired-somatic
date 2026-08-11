#Author: Brad Wubbenhorst

#Check known BRCA1/2 germline SNVs/indels against DeepVariant+VEP report CSVs.
#Samples in cnvkit_exclude_normals.list OR structural/? HGVS -> skip_cnv (not in SNV/indel pass).
#Sample ids ending in -T1, -T2, or -DZ2 are mapped to -GL for report lookup.
#Gene==WT: flag any BRCA1/2 Variant.LoF_level==1.
#Gene==BRCA1/BRCA2 with empty HGVS: report LoF1 in that gene if present.
#Gene==BRCA1/BRCA2 with HGVS: fuzzy HGVSc/HGVSp match (dup↔ins, ±1 shifts).
#
#Known TSV (tab): Sample  Gene  HGVSc  HGVSp
#  python check_known_germline_pvs.py -k known_pvs.tsv -e cnvkit_exclude_normals.list -d data/work -o known_pv_check.tsv

import argparse
import csv
import os
import re
import sys
from collections import Counter

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('-k','--known',required=True,help='TSV with Sample, Gene, HGVSc, HGVSp')
    p.add_argument('-e','--exclude',default='cnvkit_exclude_normals.list',help='Large del/dup sample list')
    p.add_argument('-d','--workdir',default='data/work',help='Work directory root')
    p.add_argument('-p','--report-path',default='{sample}/deepvariant/germline.norm.vep.report.csv',help='Report path under workdir')
    p.add_argument('-o','--output',default='known_pv_check.tsv',help='Summary TSV')
    p.add_argument('--wt-lof-out',default='known_pv_wt_lof1.tsv',help='Unexpected WT LoF1 detail TSV')
    p.add_argument('--rescue-out',default='known_pv_rescue.tsv',help='Cross-sample rescue TSV for unresolved variants')
    p.add_argument('--final-out',default='known_pv_final_summary.tsv',help='Final human-readable summary TSV')
    argv=p.parse_args()
    return vars(argv)

def map_sample(sample):
    sample=sample.strip()
    for suffix in ('-T1','-T2','-DZ2'):
        if sample.endswith(suffix):
            return sample[:-len(suffix)]+'-GL'
    return sample

def load_exclude(fp):
    out=set()
    if not os.path.exists(fp):
        print(f"WARNING: exclude list missing ({fp}); continuing with none",file=sys.stderr)
        return out
    with open(fp,'r') as fh:
        for line in fh:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            out.add(map_sample(line))
    return out

def load_known(fp):
    rows=[]
    with open(fp,'r',newline='') as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError(f"empty known file: {fp}")
        fields={f.lower():f for f in reader.fieldnames}
        if 'sample' not in fields or 'gene' not in fields:
            raise ValueError(f"known TSV needs Sample and Gene (found: {reader.fieldnames})")
        sample_k=fields['sample']
        gene_k=fields['gene']
        hgvsc_k=fields.get('hgvsc')
        hgvsp_k=fields.get('hgvsp')
        for raw in reader:
            sample=map_sample(raw[sample_k] or '')
            if not sample:
                continue
            rows.append({
                'sample':sample,
                'gene':(raw[gene_k] or '').strip().upper(),
                'hgvsc':(raw[hgvsc_k] or '').strip() if hgvsc_k else '',
                'hgvsp':(raw[hgvsp_k] or '').strip() if hgvsp_k else ''
            })
    return rows

def strip_q(s):
    return s.strip().strip('"')

def normalize_c(hgvs):
    if not hgvs or hgvs in {'.','NA','na','None'}:
        return ''
    s=strip_q(hgvs)
    if ':' in s and s.split(':',1)[1].startswith('c.'):
        s=s.split(':',1)[1]
    return s.replace(' ','')

def normalize_p(hgvs):
    if not hgvs or hgvs in {'.','NA','na','None'}:
        return ''
    s=strip_q(hgvs)
    idx=s.find('p.')
    if idx>=0:
        s=s[idx:]
    s=s.replace(' ','')
    s=re.sub(r'p\.\(([^)]+)\)',r'p.\1',s)
    return s.replace('*','Ter')

def is_structural(hgvsc):
    #Exon-scale / uncertain-boundary events are not DeepVariant SNV/indel targets.
    s=normalize_c(hgvsc)
    if not s:
        return False
    if '?' in s:
        return True
    if ';' in s:
        return True
    if re.search(r'c\.\d+[+-]\d+_\d+[+-]\d+(delins|del|dup)',s,re.I):
        return True
    return False

def cds_tokens(hgvsc):
    s=normalize_c(hgvsc)
    out={'raw':s,'kind':'','start':None,'end':None,'ins':''}
    if not s.startswith('c.'):
        return out
    body=s[2:]
    m=re.match(r'(\d+)(?:_(\d+))?(delins|del|dup|ins)([A-Z]+)?',body,re.I)
    if m:
        out['start']=int(m.group(1))
        out['end']=int(m.group(2) or m.group(1))
        out['kind']=m.group(3).lower()
        out['ins']=(m.group(4) or '').upper()
        return out
    m=re.match(r'(\d+)([ACGT])>([ACGT])',body,re.I)
    if m:
        out['start']=int(m.group(1))
        out['end']=out['start']
        out['kind']='snv'
        out['ins']=f"{m.group(2).upper()}>{m.group(3).upper()}"
    return out

def protein_tokens(hgvsp):
    s=normalize_p(hgvsp)
    out={'raw':s,'pos':None,'fs':False}
    if not s.startswith('p.'):
        return out
    body=s[2:]
    m=re.match(r'[A-Za-z]+(\d+)',body)
    if m:
        out['pos']=int(m.group(1))
    out['fs']='fs' in body.lower()
    return out

def hgvsc_equivalent(a,b):
    na,nb=normalize_c(a),normalize_c(b)
    if not na or not nb:
        return None
    if na==nb:
        return 'exact_c'
    ta,tb=cds_tokens(na),cds_tokens(nb)
    if not ta['kind'] or not tb['kind']:
        return None
    kinds={ta['kind'],tb['kind']}
    if kinds<= {'dup','ins'} and ta['start'] and tb['start']:
        if abs(ta['start']-tb['start'])<=1 and abs((ta['end'] or 0)-(tb['end'] or 0))<=1:
            if ta['ins'] and tb['ins'] and len(ta['ins'])==len(tb['ins']):
                return 'equiv_c_dup_ins'
            if not ta['ins'] or not tb['ins']:
                return 'equiv_c_dup_ins'
    if ta['kind']==tb['kind'] and ta['kind'] in {'del','delins','ins','dup'}:
        if ta['start'] and tb['start'] and abs(ta['start']-tb['start'])<=1:
            if abs((ta['end'] or ta['start'])-(tb['end'] or tb['start']))<=1:
                return 'equiv_c_indel_shift'
    if ta['kind']=='snv' and tb['kind']=='snv' and ta['start']==tb['start'] and ta['ins']==tb['ins']:
        return 'exact_c'
    return None

def hgvsp_equivalent(a,b):
    na,nb=normalize_p(a),normalize_p(b)
    if not na or not nb:
        return None
    if na==nb:
        return 'exact_p'
    pa,pb=protein_tokens(na),protein_tokens(nb)
    if pa['pos'] and pb['pos'] and pa['fs'] and pb['fs'] and abs(pa['pos']-pb['pos'])<=1:
        return 'equiv_p_fs_shift'
    return None

def load_brca_rows(csv_path):
    if not os.path.exists(csv_path):
        return []
    rows=[]
    with open(csv_path,'r',newline='') as fh:
        reader=csv.DictReader(fh)
        for raw in reader:
            gene=strip_q(raw.get('Gene',''))
            if gene not in {'BRCA1','BRCA2'}:
                continue
            rows.append({
                'gene':gene,
                'hgvsc':strip_q(raw.get('HGVSc','')),
                'hgvsp':strip_q(raw.get('HGVSp','')),
                'lof':strip_q(raw.get('Variant.LoF_level','')),
                'consequence':strip_q(raw.get('Variant.Consequence','')),
                'clinvar_sig':strip_q(raw.get('ClinVar.SIG','')),
                'filter':strip_q(raw.get('FILTER','')),
                'chr':strip_q(raw.get('Chr','')),
                'start':strip_q(raw.get('Start','')),
                'ref':strip_q(raw.get('REF','')),
                'alt':strip_q(raw.get('ALT','')),
                'zyg':strip_q(raw.get('Sample.Zyg','')),
                'autogvp':strip_q(raw.get('AutoGVP',''))
            })
    return rows

def lof1_brca(rows,gene=None):
    out=[r for r in rows if r['lof']=='1']
    if gene is not None:
        out=[r for r in out if r['gene']==gene]
    return out

def match_in_rows(known,candidates):
    tier_rank={
        'exact_c':1,
        'exact_p':2,
        'equiv_c_dup_ins':3,
        'equiv_c_indel_shift':4,
        'equiv_p_fs_shift':5
    }
    best=None
    best_tier=None
    for r in candidates:
        tiers=[]
        tc=hgvsc_equivalent(known['hgvsc'],r['hgvsc'])
        tp=hgvsp_equivalent(known['hgvsp'],r['hgvsp'])
        if tc:
            tiers.append(tc)
        if tp:
            tiers.append(tp)
        if not tiers:
            continue
        tier=sorted(tiers,key=lambda t:tier_rank.get(t,99))[0]
        if best is None or tier_rank[tier]<tier_rank[best_tier]:
            best=r
            best_tier=tier
    return best,best_tier

def best_match(known,rows):
    gene=known['gene']
    candidates=[r for r in rows if r['gene']==gene]
    best,best_tier=match_in_rows(known,candidates)
    if best is not None:
        status='found' if best_tier in {'exact_c','exact_p'} else 'likely_equivalent'
        return status,best,best_tier

    #Known table sometimes has the wrong BRCA gene; try the other one.
    other='BRCA2' if gene=='BRCA1' else 'BRCA1'
    other_rows=[r for r in rows if r['gene']==other]
    best,best_tier=match_in_rows(known,other_rows)
    if best is not None:
        status='found_other_gene' if best_tier in {'exact_c','exact_p'} else 'likely_equivalent_other_gene'
        return status,best,f"{best_tier}|known_gene={gene}|hit_gene={other}"

    if not candidates:
        return 'not_found',None,'no BRCA rows for gene'
    lof1=lof1_brca(candidates)
    if len(lof1)==1:
        return 'review_lof1',lof1[0],'single_lof1_in_gene'
    if len(lof1)>1:
        return 'review_lof1',lof1[0],f"{len(lof1)}_lof1_in_gene"
    return 'not_found',None,'no HGVS or LoF1 match'

def fmt_hit(row):
    if not row:
        return ''
    return (
        f"{row['gene']}|{row['hgvsc']}|{row['hgvsp']}|LoF={row['lof']}|"
        f"{row['chr']}:{row['start']} {row['ref']}>{row['alt']}|{row['zyg']}"
    )

def fmt_hits(rows):
    if not rows:
        return ''
    return ' ; '.join(fmt_hit(r) for r in rows)

def append_wt_lof(wt_lof_rows,sample,l1):
    for r in l1:
        wt_lof_rows.append({
            'sample':sample,
            'gene':r['gene'],
            'hgvsc':r['hgvsc'],
            'hgvsp':r['hgvsp'],
            'lof':r['lof'],
            'consequence':r['consequence'],
            'clinvar_sig':r['clinvar_sig'],
            'filter':r['filter'],
            'variant':f"{r['chr']}:{r['start']}:{r['ref']}>{r['alt']}",
            'zyg':r['zyg'],
            'autogvp':r['autogvp']
        })

def all_report_samples(workdir):
    out=[]
    if not os.path.exists(workdir):
        return out
    for sample in sorted(os.listdir(workdir)):
        rep=os.path.join(workdir,sample,'deepvariant','germline.norm.vep.report.csv')
        if os.path.exists(rep):
            out.append((sample,rep))
    return out

def build_report_cache(workdir,exclude):
    cache={}
    for sample,rep in all_report_samples(workdir):
        cache[sample]={
            'report_path':rep,
            'excluded':sample in exclude,
            'rows':load_brca_rows(rep)
        }
    return cache

def rescue_candidates(known,report_cache):
    out=[]
    for sample,data in report_cache.items():
        rows=data['rows']
        gene_rows=[r for r in rows if r['gene']==known['gene']]
        other_gene='BRCA2' if known['gene']=='BRCA1' else 'BRCA1'
        other_rows=[r for r in rows if r['gene']==other_gene]
        best,best_tier=match_in_rows(known,gene_rows)
        if best is not None:
            out.append({
                'candidate_sample':sample,
                'candidate_status':'same_gene_match',
                'match_detail':best_tier,
                'excluded':'yes' if data['excluded'] else 'no',
                'report_path':data['report_path'],
                'hit':fmt_hit(best)
            })
            continue
        best,best_tier=match_in_rows(known,other_rows)
        if best is not None:
            out.append({
                'candidate_sample':sample,
                'candidate_status':'other_gene_match',
                'match_detail':best_tier,
                'excluded':'yes' if data['excluded'] else 'no',
                'report_path':data['report_path'],
                'hit':fmt_hit(best)
            })
            continue
        lof1=lof1_brca(gene_rows)
        if lof1:
            out.append({
                'candidate_sample':sample,
                'candidate_status':'same_gene_lof1',
                'match_detail':f"{len(lof1)}_lof1",
                'excluded':'yes' if data['excluded'] else 'no',
                'report_path':data['report_path'],
                'hit':' ; '.join(fmt_hit(r) for r in lof1[:3])
            })
    return out

def main(argv=None):
    argv=get_args() if argv is None else argv
    known_rows=load_known(argv['known'])
    exclude=load_exclude(argv['exclude'])
    report_cache=build_report_cache(argv['workdir'],exclude)
    out_fields=[
        'sample','known_gene','known_hgvsc','known_hgvsp','status','match_detail',
        'hit','report_path','n_brca_rows','n_brca_lof1'
    ]
    results=[]
    wt_lof_rows=[]
    rescue_rows=[]
    final_rows=[]

    for kn in known_rows:
        sample=kn['sample']
        gene=kn['gene']
        rep=os.path.join(argv['workdir'],argv['report_path'].format(sample=sample))
        base={
            'sample':sample,
            'known_gene':gene,
            'known_hgvsc':kn['hgvsc'],
            'known_hgvsp':kn['hgvsp'],
            'report_path':rep,
            'n_brca_rows':'',
            'n_brca_lof1':'',
            'hit':'',
            'match_detail':''
        }

        if sample in exclude:
            base['status']='skip_cnv'
            base['match_detail']='in cnvkit_exclude_normals.list'
            results.append(base)
            final_rows.append({
                'sample':sample,
                'group':'skip_cnv',
                'expected_gene':gene,
                'expected_hgvsc':kn['hgvsc'],
                'expected_hgvsp':kn['hgvsp'],
                'found_status':base['status'],
                'found_match_detail':base['match_detail'],
                'found_variant':'',
                'other_brca_lof1':'',
                'report_path':rep
            })
            continue

        if is_structural(kn['hgvsc']):
            base['status']='skip_cnv'
            base['match_detail']='structural_or_uncertain_hgvs'
            results.append(base)
            final_rows.append({
                'sample':sample,
                'group':'skip_cnv',
                'expected_gene':gene,
                'expected_hgvsc':kn['hgvsc'],
                'expected_hgvsp':kn['hgvsp'],
                'found_status':base['status'],
                'found_match_detail':base['match_detail'],
                'found_variant':'',
                'other_brca_lof1':'',
                'report_path':rep
            })
            continue

        if sample in report_cache:
            brca=report_cache[sample]['rows']
        elif not os.path.exists(rep):
            base['status']='no_report'
            base['match_detail']='missing report csv'
            results.append(base)
            final_rows.append({
                'sample':sample,
                'group':'no_report',
                'expected_gene':gene,
                'expected_hgvsc':kn['hgvsc'],
                'expected_hgvsp':kn['hgvsp'],
                'found_status':base['status'],
                'found_match_detail':base['match_detail'],
                'found_variant':'',
                'other_brca_lof1':'',
                'report_path':rep
            })
            continue
        else:
            brca=load_brca_rows(rep)
        l1=lof1_brca(brca)
        all_lof1=fmt_hits(l1)
        base['n_brca_rows']=str(len(brca))
        base['n_brca_lof1']=str(len(l1))

        if gene=='WT':
            if l1:
                base['status']='wt_unexpected_lof1'
                base['match_detail']=f"{len(l1)} BRCA1/2 LoF_level=1"
                base['hit']=' ; '.join(fmt_hit(r) for r in l1)
                results.append(base)
                append_wt_lof(wt_lof_rows,sample,l1)
                final_rows.append({
                    'sample':sample,
                    'group':'wt_with_lof',
                    'expected_gene':gene,
                    'expected_hgvsc':kn['hgvsc'],
                    'expected_hgvsp':kn['hgvsp'],
                    'found_status':base['status'],
                    'found_match_detail':base['match_detail'],
                    'found_variant':'',
                    'other_brca_lof1':all_lof1,
                    'report_path':rep
                })
            else:
                base['status']='wt_ok'
                base['match_detail']='no BRCA1/2 LoF_level=1'
                results.append(base)
                final_rows.append({
                    'sample':sample,
                    'group':'wt_ok',
                    'expected_gene':gene,
                    'expected_hgvsc':kn['hgvsc'],
                    'expected_hgvsp':kn['hgvsp'],
                    'found_status':base['status'],
                    'found_match_detail':base['match_detail'],
                    'found_variant':'',
                    'other_brca_lof1':'',
                    'report_path':rep
                })
            continue

        if gene not in {'BRCA1','BRCA2'}:
            base['status']='bad_gene'
            base['match_detail']=f"unexpected Gene={gene}"
            results.append(base)
            final_rows.append({
                'sample':sample,
                'group':'bad_gene',
                'expected_gene':gene,
                'expected_hgvsc':kn['hgvsc'],
                'expected_hgvsp':kn['hgvsp'],
                'found_status':base['status'],
                'found_match_detail':base['match_detail'],
                'found_variant':'',
                'other_brca_lof1':all_lof1,
                'report_path':rep
            })
            continue

        if not kn['hgvsc'] and not kn['hgvsp']:
            gene_l1=lof1_brca(brca,gene=gene)
            if gene_l1:
                base['status']='carrier_no_hgvs_lof1'
                base['match_detail']=f"{len(gene_l1)} LoF_level=1 in {gene}"
                base['hit']=' ; '.join(fmt_hit(r) for r in gene_l1)
            else:
                base['status']='carrier_no_hgvs_no_lof1'
                base['match_detail']=f"no LoF_level=1 in {gene}"
            results.append(base)
            final_rows.append({
                'sample':sample,
                'group':'carrier_no_known_hgvs',
                'expected_gene':gene,
                'expected_hgvsc':kn['hgvsc'],
                'expected_hgvsp':kn['hgvsp'],
                'found_status':base['status'],
                'found_match_detail':base['match_detail'],
                'found_variant':base['hit'],
                'other_brca_lof1':all_lof1,
                'report_path':rep
            })
            continue

        status,hit,detail=best_match(kn,brca)
        base['status']=status
        base['match_detail']=detail
        base['hit']=fmt_hit(hit)
        results.append(base)

        if status=='found':
            group='exact_match'
        elif status in {'likely_equivalent','found_other_gene','likely_equivalent_other_gene'}:
            group='close_match'
        elif l1:
            group='missing_expected_other_lof'
        else:
            group='missing_expected_no_lof'
        final_rows.append({
            'sample':sample,
            'group':group,
            'expected_gene':gene,
            'expected_hgvsc':kn['hgvsc'],
            'expected_hgvsp':kn['hgvsp'],
            'found_status':status,
            'found_match_detail':detail,
            'found_variant':base['hit'],
            'other_brca_lof1':all_lof1,
            'report_path':rep
        })

        if status in {'not_found','review_lof1','found_other_gene','likely_equivalent_other_gene'}:
            if kn['hgvsc'] or kn['hgvsp']:
                for cand in rescue_candidates(kn,report_cache):
                    rescue_rows.append({
                        'expected_sample':sample,
                        'expected_gene':gene,
                        'expected_hgvsc':kn['hgvsc'],
                        'expected_hgvsp':kn['hgvsp'],
                        'expected_status':status,
                        'candidate_sample':cand['candidate_sample'],
                        'candidate_status':cand['candidate_status'],
                        'match_detail':cand['match_detail'],
                        'excluded':cand['excluded'],
                        'hit':cand['hit'],
                        'report_path':cand['report_path']
                    })

    with open(argv['output'],'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=out_fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(results)

    if wt_lof_rows:
        with open(argv['wt_lof_out'],'w',newline='') as fh:
            w=csv.DictWriter(fh,fieldnames=list(wt_lof_rows[0].keys()),delimiter='\t')
            w.writeheader()
            w.writerows(wt_lof_rows)
    elif os.path.exists(argv['wt_lof_out']):
        os.remove(argv['wt_lof_out'])

    if rescue_rows:
        rescue_fields=[
            'expected_sample','expected_gene','expected_hgvsc','expected_hgvsp','expected_status',
            'candidate_sample','candidate_status','match_detail','excluded','hit','report_path'
        ]
        with open(argv['rescue_out'],'w',newline='') as fh:
            w=csv.DictWriter(fh,fieldnames=rescue_fields,delimiter='\t')
            w.writeheader()
            w.writerows(rescue_rows)
    elif os.path.exists(argv['rescue_out']):
        os.remove(argv['rescue_out'])

    final_fields=[
        'sample','group','expected_gene','expected_hgvsc','expected_hgvsp',
        'found_status','found_match_detail','found_variant','other_brca_lof1','report_path'
    ]
    with open(argv['final_out'],'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=final_fields,delimiter='\t')
        w.writeheader()
        w.writerows(final_rows)

    counts=Counter(r['status'] for r in results)
    print(f"Wrote {argv['output']} ({len(results)} rows)")
    for k,v in sorted(counts.items()):
        print(f"  {k}: {v}")

    wt_rows=[r for r in results if r['known_gene']=='WT']
    print(f"\nWT BRCA1/2 LoF1 checks ({len(wt_rows)} samples):")
    if not wt_rows:
        print("  (no Gene=WT rows in known TSV)")
    else:
        wt_counts=Counter(r['status'] for r in wt_rows)
        for k,v in sorted(wt_counts.items()):
            print(f"  {k}: {v}")
        if wt_lof_rows:
            print(f"  detail: {argv['wt_lof_out']} ({len(wt_lof_rows)} rows)")

    glance_status={
        'not_found','likely_equivalent','likely_equivalent_other_gene','found_other_gene',
        'review_lof1','wt_unexpected_lof1','carrier_no_hgvs_lof1'
    }
    review=[r for r in results if r['status'] in glance_status]
    if review:
        print("\nNeeds glance:")
        for r in review:
            print(f"  {r['sample']}\t{r['known_gene']}\t{r['status']}\t{r['match_detail']}\t{r['hit'][:120]}")
    if rescue_rows:
        print(f"\nCross-sample rescue: {argv['rescue_out']} ({len(rescue_rows)} rows)")
    print(f"Final summary: {argv['final_out']} ({len(final_rows)} rows)")

if __name__=='__main__':
    main()
