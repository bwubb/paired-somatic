#Author: Brad Wubbenhorst
#
#Check known BRCA1/2 germline exon dels against CNVkit germline .cnr / .call.cns.
#Gene-level only (not exon-resolved). One summary row per known event.
#
#  python check_known_germline_cnvs.py -k cnvkit_exclude_normals.tsv -d data/work -o known_cnv_check.tsv

import argparse
import csv
import math
import os
import sys

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('-k','--known',default='cnvkit_exclude_normals.tsv',help='TSV: GL ID, Gene, Event')
    p.add_argument('-d','--workdir',default='data/work',help='Work directory root')
    p.add_argument('-o','--output',default='known_cnv_check.tsv',help='Summary TSV')
    p.add_argument('--low',type=float,default=-0.4,help='log2 threshold for low bins')
    argv=p.parse_args()
    return vars(argv)

def mean(xs):
    return sum(xs)/len(xs) if xs else float('nan')

def median(xs):
    if not xs:
        return float('nan')
    ys=sorted(xs)
    n=len(ys)
    mid=n//2
    return ys[mid] if n%2 else 0.5*(ys[mid-1]+ys[mid])

def fmt(x,nd=4):
    if x is None or (isinstance(x,float) and (math.isnan(x) or math.isinf(x))):
        return 'NA'
    if isinstance(x,float):
        return f'{x:.{nd}f}'
    return str(x)

def load_known(fp):
    rows=[]
    with open(fp,'r',newline='') as fh:
        reader=csv.reader(fh,delimiter='\t')
        header=next(reader,None)
        if header is None:
            return rows
        #Header like "GL ID, Gene, Event" or bare rows starting with *-GL.
        h0=header[0].strip().lower().replace(' ','_')
        first_is_header=h0 in {'gl_id','sample','sample_id','id'} or h0.startswith('gl')
        if not first_is_header:
            if len(header)>=3:
                rows.append({'sample':header[0].strip(),'gene':header[1].strip().upper(),'event':header[2].strip()})
        for parts in reader:
            if not parts or not parts[0].strip() or parts[0].strip().startswith('#'):
                continue
            if len(parts)<3:
                continue
            rows.append({'sample':parts[0].strip(),'gene':parts[1].strip().upper(),'event':parts[2].strip()})
    return rows

def cnr_gene_stats(cnr_fp,gene,low_cut):
    vals=[]
    depths=[]
    if not os.path.exists(cnr_fp):
        return None
    with open(cnr_fp,'r',newline='') as fh:
        r=csv.DictReader(fh,delimiter='\t')
        for row in r:
            if row.get('gene','')!=gene:
                continue
            try:
                vals.append(float(row['log2']))
                depths.append(float(row['depth']))
            except (KeyError,ValueError):
                continue
    if not vals:
        return {'cnr_status':'no_gene_bins','n_bins':0,'log2_mean':'NA','log2_median':'NA','log2_min':'NA','n_low':'NA','frac_low':'NA','depth_mean':'NA'}
    n_low=sum(1 for v in vals if v<=low_cut)
    return {
        'cnr_status':'ok',
        'n_bins':len(vals),
        'log2_mean':fmt(mean(vals)),
        'log2_median':fmt(median(vals)),
        'log2_min':fmt(min(vals)),
        'n_low':n_low,
        'frac_low':fmt(n_low/len(vals)),
        'depth_mean':fmt(mean(depths),1),
    }

def call_gene_segments(cns_fp,gene):
    hits=[]
    if not os.path.exists(cns_fp):
        return hits
    with open(cns_fp,'r',newline='') as fh:
        r=csv.DictReader(fh,delimiter='\t')
        for row in r:
            if gene not in row.get('gene','').split(','):
                continue
            try:
                start=int(row['start'])
                end=int(row['end'])
                log2=float(row['log2'])
            except (KeyError,ValueError):
                continue
            hits.append({
                'chrom':row['chromosome'],
                'start':start,
                'end':end,
                'length':end-start+1,
                'log2':log2,
                'cn':row.get('cn',''),
                'cn1':row.get('cn1',''),
                'cn2':row.get('cn2',''),
                'probes':row.get('probes',''),
            })
    return hits

def classify(cnr,seg_hits,low_cut):
    if cnr is None:
        return 'missing_cnr'
    if cnr['cnr_status']=='no_gene_bins':
        return 'no_gene_bins'
    loss_segs=[h for h in seg_hits if str(h['cn']).isdigit() and int(h['cn'])<2]
    if loss_segs:
        return 'segment_cn_loss'
    try:
        log2_mean=float(cnr['log2_mean'])
        frac_low=float(cnr['frac_low'])
        log2_min=float(cnr['log2_min'])
    except ValueError:
        return 'cnr_parse_error'
    if log2_mean<=low_cut or frac_low>=0.3:
        return 'cnr_gene_low'
    if log2_min<=low_cut:
        return 'cnr_focal_dip'
    if any(str(h['cn'])=='2' for h in seg_hits):
        return 'segment_diploid_only'
    if not seg_hits:
        return 'no_segment_hit'
    return 'no_clear_loss'

def empty_cnr():
    return {'n_bins':'NA','log2_mean':'NA','log2_median':'NA','log2_min':'NA','n_low':'NA','frac_low':'NA','depth_mean':'NA'}

def main():
    argv=get_args()
    known=load_known(argv['known'])
    if not known:
        raise SystemExit(f'no known rows in {argv["known"]}')
    fields=['sample','gene','event','status','cnr_path','n_bins','log2_mean','log2_median','log2_min','n_low','frac_low','depth_mean','n_segments','seg_cn_values','seg_min_log2','seg_max_length','note']
    out_rows=[]
    counts={}
    for row in known:
        sample=row['sample']
        gene=row['gene']
        cnr_fp=os.path.join(argv['workdir'],sample,'cnvkit',f'{sample}.germline.cnr')
        cns_fp=os.path.join(argv['workdir'],sample,'cnvkit',f'{sample}.germline.call.cns')
        note=''
        if gene not in {'BRCA1','BRCA2'}:
            status='bad_gene'
            cnr=None
            segs=[]
            note=f'gene={gene}'
        else:
            cnr=cnr_gene_stats(cnr_fp,gene,argv['low'])
            segs=call_gene_segments(cns_fp,gene)
            status=classify(cnr,segs,argv['low'])
            if cnr is None:
                note='missing germline.cnr'
            elif not os.path.exists(cns_fp):
                note='missing germline.call.cns'
        counts[status]=counts.get(status,0)+1
        if cnr is None:
            cnr=empty_cnr()
        cn_vals=','.join(sorted({str(s['cn']) for s in segs})) if segs else 'NA'
        seg_min_log2=fmt(min(s['log2'] for s in segs)) if segs else 'NA'
        seg_max_len=max(s['length'] for s in segs) if segs else 'NA'
        out_rows.append({
            'sample':sample,
            'gene':gene,
            'event':row['event'],
            'status':status,
            'cnr_path':cnr_fp if os.path.exists(cnr_fp) else 'MISSING',
            'n_bins':cnr['n_bins'],
            'log2_mean':cnr['log2_mean'],
            'log2_median':cnr['log2_median'],
            'log2_min':cnr['log2_min'],
            'n_low':cnr['n_low'],
            'frac_low':cnr['frac_low'],
            'depth_mean':cnr['depth_mean'],
            'n_segments':len(segs),
            'seg_cn_values':cn_vals,
            'seg_min_log2':seg_min_log2,
            'seg_max_length':seg_max_len,
            'note':note,
        })
    with open(argv['output'],'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t')
        w.writeheader()
        w.writerows(out_rows)
    print(f'Wrote {argv["output"]} ({len(out_rows)} rows)')
    print('Status counts:')
    for k,v in sorted(counts.items(),key=lambda kv:(-kv[1],kv[0])):
        print(f'  {k}\t{v}')

if __name__=='__main__':
    main()
