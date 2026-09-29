#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
  Comprehensive sedaDNA Species Authenticity Evaluation
===============================================================================

This script scores species-level sedaDNA candidates from ngsLCA + bamdam
outputs and summarizes their authenticity confidence.

Scoring dimensions (100 points total):
  A. Terminal damage patterns (35 points)
  B. Coverage uniformity and read support (30 points)
  C. Taxonomic specificity (35 points)

Hard filters:
  - 5' C->T < 5%: damage score becomes 0 and total score loses 10 points
  - Mean fragment length > 100 bp: total score loses 10 points

Confidence tiers:
  >= 65  HIGH
  >= 45  MEDIUM
  >= 30  LOW
  >= 12  VERY_LOW
  <  12  BACKGROUND

Required inputs:
  --bamdam_tsv      species-level bamdam TSV
  --bamdam_subs     species-level bamdam .subs file
  --genus_tsv       genus-level bamdam TSV
  --family_tsv      family-level bamdam TSV

Optional inputs:
  --coverage_dir    directory containing samtools coverage outputs
  --abundance_tsv   optional species abundance table

Outputs:
  --output          simplified TSV report
  --output_detailed detailed TSV report with intermediate metrics
===============================================================================
"""
import os, sys, re, glob, math, csv, argparse
from collections import defaultdict

def safe_float(val, default=0.0):
    try: return float(val)
    except: return default

def get_rank_taxid(taxpath, rank_name):
    for rs in taxpath.split(';'):
        p = rs.split(':')
        if len(p) >= 3 and p[2].strip() == rank_name:
            return p[0].strip()
    return None

def linreg(xs, ys):
    n = len(xs)
    if n < 2: return 0.0, 0.0, 0.0
    mx, my = sum(xs)/n, sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    den = sum((x-mx)**2 for x in xs)
    if den == 0: return 0.0, my, 0.0
    slope = num/den
    intercept = my - slope*mx
    ss_res = sum((y-(slope*x+intercept))**2 for x,y in zip(xs,ys))
    ss_tot = sum((y-my)**2 for y in ys)
    r2 = 1 - (ss_res/ss_tot) if ss_tot > 0 else 0.0
    return slope, intercept, r2

def gini(vals):
    if len(vals) < 2: return 0.0
    sv = sorted(vals)
    n, total = len(sv), sum(sv)
    if total == 0: return 0.0
    ws = sum((i+1)*v for i,v in enumerate(sv))
    return (2*ws - (n+1)*total) / (n*total)

def parse_subs_line(line):
    parts = line.strip().split('\t')
    if len(parts) < 2: return None
    taxid = parts[0]
    kv = re.compile(r'([ACGT][ACGT])(-?\d+):([\d.]+)')
    fd = {}
    for field in parts[1:]:
        for m in kv.finditer(field):
            fd[(m.group(1), int(m.group(2)))] = float(m.group(3))
    p5, p3 = {}, {}
    for (sub, pos), freq in fd.items():
        if pos < 0 and sub in ('CC','CT'):
            denom = fd.get(('CC',pos),0) + fd.get(('CT',pos),0)
            if denom > 0: p5[pos] = fd.get(('CT',pos),0)/denom
        elif pos > 0 and sub in ('GG','GA'):
            denom = fd.get(('GG',pos),0) + fd.get(('GA',pos),0)
            if denom > 0: p3[pos] = fd.get(('GA',pos),0)/denom
    return {'taxid':taxid,'5p':p5,'3p':p3}

def damage_metrics(d5_raw, d3_raw, subs):
    """d5_raw, d3_raw: 直接来自tsv的Damage+1/Damage-1列(fraction 0-1)
       subs: 来自.subs文件的per-position数据，仅用于衰减模式"""
    p5 = subs.get('5p',{}); p3 = subs.get('3p',{})
    if not p5:
        return {'d5':d5_raw,'d3':d3_raw,'asym':0,'slope':0,'r2':0,'signal':0,'plat':0}
    d5v = [p5.get(p,0) for p in range(-15,0)]
    asym = d5_raw / max(d3_raw, 0.001)
    slope, _, r2 = linreg(list(range(1,16)), d5v)
    signal = sum(d5v[:3])/3 - sum(d5v[-3:])/3
    plat = sum(d5v[-5:])/5
    return {'d5':d5_raw,'d3':d3_raw,'asym':asym,'slope':slope,'r2':r2,'signal':signal,'plat':plat}

def cov_metrics(taxid, cov_dir):
    pat = os.path.join(cov_dir, f'{taxid}_*_coverage.txt')
    files = glob.glob(pat)
    default = {'breadth':0,'mdepth':0,'gsize':0,'nctg':0,'ncov':0,'cv':0,'gini':0,'ctgrat':0,'maxsc':0,'has':False}
    if not files: return default
    dps, lens, ncov, gs, tcb, msc = [], [], 0, 0, 0, 0
    try:
        with open(files[0]) as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.strip().split('\t')
                if len(p) < 7: continue
                try: sl, nr, cb, cp, md = int(p[2]),int(p[3]),int(p[4]),float(p[5]),float(p[6])
                except: continue
                lens.append(sl); dps.append(md); gs+=sl; tcb+=cb
                if nr>0: ncov+=1
                if cp>msc: msc=cp
    except: return default
    nctg = len(lens)
    if nctg==0 or gs==0: return default
    breadth = (tcb/gs)*100
    mdepth = sum(d*l for d,l in zip(dps,lens))/gs
    cv = math.sqrt(sum((d-mdepth)**2 for d in dps)/nctg)/mdepth if mdepth>0 and nctg>1 else 0
    g = gini(dps)
    ctr = (ncov/nctg)*100
    return {'breadth':breadth,'mdepth':mdepth,'gsize':gs,'nctg':nctg,'ncov':ncov,'cv':cv,'gini':g,'ctgrat':ctr,'maxsc':msc,'has':True}

def spec_metrics(row, gd, fd, family_groups):
    tr = safe_float(row.get('TotalReads',0))
    ur = safe_float(row.get('UnaggregatedReads',0))
    ani = safe_float(row.get('ANI',1.0))
    ukpr = safe_float(row.get('UniqKmersPerRead',0))
    rgc = safe_float(row.get('AvgReadGC',0))
    fgc = safe_float(row.get('AvgRefGC',0))
    ta  = safe_float(row.get('TotalAlignments',tr))
    dup = safe_float(row.get('Duplicity',0))
    dust = safe_float(row.get('MeanDust',0))
    lca = (ur/tr)*100 if tr>0 else 100.0
    align_ratio = ta/tr if tr>0 else 1.0
    tp = row.get('taxpath','')
    gid = get_rank_taxid(tp,'genus')
    fid = get_rank_taxid(tp,'family')
    gl, fl = 0, 0
    if gid and gid in gd:
        gt=safe_float(gd[gid].get('TotalReads',1)); gu=safe_float(gd[gid].get('UnaggregatedReads',0))
        if gt>0: gl=(gu/gt)*100
    if fid and fid in fd:
        ft=safe_float(fd[fid].get('TotalReads',1)); fu=safe_float(fd[fid].get('UnaggregatedReads',0))
        if ft>0: fl=(fu/ft)*100

    # === 分类独特性: Read Sink 检测 ===
    # 同一科内，该物种reads占比 × LCA降解率 = sink_score (0~1)
    # sink_score高 = 该物种吸收了科内大量reads但无法种水平确认 → 参考偏差"汇集点"
    family_total = 0; family_species_count = 0
    if fid and fid in family_groups:
        family_total = family_groups[fid]['total_reads']
        family_species_count = family_groups[fid]['species_count']
    read_share = tr / family_total if family_total > 0 else 0
    lca_pct = lca / 100.0
    sink_score = read_share * lca_pct

    return {'lca':lca,'ani':ani,'ukpr':ukpr,'gc_bias':abs(rgc-fgc),
            'glca':gl,'flca':fl,'dup':dup,'dust':dust,'align_ratio':align_ratio,
            'read_share':read_share,'sink_score':sink_score,'family_total':family_total,
            'family_species_count':family_species_count}

def ab_metrics(taxid, abd, gsize):
    ab = abd.get(taxid,{})
    reads = safe_float(ab.get('reads',0))
    rel = safe_float(ab.get('relative_abundance',0))
    rpm = reads/(gsize/1e6) if gsize>0 else 0
    return {'rel':rel*100,'reads':reads,'rpm':rpm}

def score(m):
    dmg = m['damage']; cov = m['coverage']; spec = m['specificity']
    mlen = m['length']['mean_length']; tr = m['total_reads']
    ss = {}
    # ============================================================
    # A. 末端损伤模式 Damage (35分) -- 古DNA最核心特征
    #    5' C→T < 5% 直接pass，损伤分为0
    # ============================================================
    ds = 0
    d5 = dmg['d5']*100
    d5_pass = (d5 >= 5)   # 硬门槛：低于5%不算古DNA损伤
    if d5_pass:
        # 5' C→T 损伤率 (0-14分) -- 略微放宽，沉积物中15%+已是很强的古DNA信号
        if d5>=30: ds+=14
        elif d5>=20: ds+=12
        elif d5>=15: ds+=11
        elif d5>=10: ds+=8
        elif d5>=5: ds+=5
        # 损伤衰减模式: 沉积物中衰减信号可能弱甚至为负(噪音干扰)，放宽门槛 (0-10分)
        sig = dmg['signal']*100; r2 = dmg['r2']
        if sig>=4 and r2>0.2: ds+=10
        elif sig>=2 and r2>0.05: ds+=8
        elif sig>=0.5: ds+=6
        elif r2>0.1: ds+=4          # 衰减方向不明确但R²尚可，给基础分
        # 5'/3' 不对称性 (0-6分)
        asym = dmg['asym']
        if asym>4: ds+=6
        elif asym>2.5: ds+=5
        elif asym>1.8: ds+=4
        elif asym>1.3: ds+=3
        elif asym>1.0: ds+=1
        # 3' G→A 损伤 (0-5分)
        d3 = dmg['d3']*100
        if d3>=8: ds+=5
        elif d3>=5: ds+=4
        elif d3>=2: ds+=2
        elif d3>0: ds+=1
    # 3'/5' 对称性惩罚: 3'损伤接近5'→非古DNA末端模式(参考偏差)
    if d5_pass:
        d3_ratio = dmg['d3'] / max(dmg['d5'], 0.001)
        if d3_ratio > 0.6:
            if tr >= 50000: ds = max(0, ds - 7)
            elif tr >= 10000: ds = max(0, ds - 4)
            elif tr >= 2000: ds = max(0, ds - 2)
    ss['damage'] = min(35,ds)
    ss['d5_pass'] = d5_pass

    # ============================================================
    # B. 片段长度 Length (定性，不计入总分)
    #    过长片段(>100bp)基本不可能是古DNA
    # ============================================================
    ls = 0 if mlen > 100 else 1    # 1=正常古DNA长度范围, 0=过长存疑
    ss['length'] = ls

    # ============================================================
    # C. 覆盖均匀度 + Reads支持度 Coverage (30分)
    #    覆盖广度放宽(沉积物aDNA天然低覆盖) + reads数作为可靠性加权
    # ============================================================
    cs = 0
    br = cov['breadth']
    # 放宽广度阈值: 沉积物中0.001%已有意义
    if br>=1: cs+=10
    elif br>=0.1: cs+=8
    elif br>=0.01: cs+=5
    elif br>=0.001: cs+=3
    elif br>0: cs+=1
    gi = cov['gini']; cr = cov.get('ctgrat',0)
    if gi<0.3 and cr>50: cs+=9
    elif gi<0.5 and cr>30: cs+=7
    elif gi<0.7 and cr>10: cs+=5
    elif gi<0.85 and cr>5: cs+=3
    elif cr>1: cs+=1
    md = cov['mdepth']
    if md>=0.01: cs+=5
    elif md>=0.001: cs+=3
    elif md>0: cs+=1
    # Reads数量作为覆盖统计的可靠性因子
    if tr>=20000: cs+=6
    elif tr>=5000: cs+=5
    elif tr>=2000: cs+=4
    elif tr>=500: cs+=3
    elif tr>=100: cs+=2
    # 无coverage文件时给基线分(避免因缺少参考基因组而过度惩罚)
    # 沉积物aDNA很多物种参考基因组不完整是常态，提高基线
    if not cov['has']:
        cs += 6
    ss['coverage'] = min(30,cs)

    # ============================================================
    # D. 分类特异性 Specificity (35分)
    #    A. 分类独特性 (0-15分): Read Sink检测 + LCA降解
    #    B. 序列匹配质量 (0-12分): ANI + GC偏差 + 多重比对
    #    C. Reads质量 (0-8分): Kmer特异性 + PCR重复率 + 序列复杂度
    # ============================================================
    sps = 0

    # --- A. 分类独特性 Taxonomic Uniqueness (0-15分) ---
    # 核心理念: 如果同一科内大量reads都集中到一个物种上，
    # 且该物种LCA降解率极高(无法种水平确认)，这是参考偏差(read sink)
    lca = spec['lca']
    sink = spec['sink_score']       # read_share × lca_pct (0~1)
    family_sp_count = spec.get('family_species_count',0)
    read_share = spec.get('read_share',0)

    uniq = 15
    # Read Sink 惩罚: 吸收了科内大部分reads却无法种水平确认
    if sink > 0.95: uniq -= 12
    elif sink > 0.7: uniq -= 7
    elif sink > 0.4: uniq -= 4
    elif sink > 0.15: uniq -= 2
    # LCA降级率直接惩罚
    if lca > 99: uniq -= 3
    elif lca > 95: uniq -= 2
    elif lca > 80: uniq -= 1
    sps += max(0, uniq)

    # --- B. 序列匹配质量 Sequence Match (0-12分) ---
    # ANI (0-7分) -- 沉积物metagenomics中ANI>=0.97已是很好的匹配
    ani = spec['ani']
    if ani>=0.99: sps+=7
    elif ani>=0.985: sps+=6
    elif ani>=0.98: sps+=5
    elif ani>=0.97: sps+=3
    elif ani>=0.95: sps+=2
    elif ani>=0.90: sps+=1
    # GC偏差 (0-3分)
    gb = spec['gc_bias']
    if gb<0.01: sps+=3
    elif gb<0.02: sps+=2
    elif gb<0.05: sps+=1
    # 多重比对比例 (0-2分)
    ar = spec['align_ratio']
    if ar<=1.05: sps+=2
    elif ar<=1.2: sps+=1

    # --- C. Reads质量 Read Quality (0-8分) ---
    # Kmer特异性 (0-3分)
    ukpr = spec['ukpr']
    if ukpr>=0.8: sps+=3
    elif ukpr>=0.6: sps+=2
    elif ukpr>=0.4: sps+=1
    # PCR重复率 (0-3分)
    dup = spec['dup']
    if dup<5: sps+=3
    elif dup<20: sps+=2
    elif dup<50: sps+=1
    # 序列复杂度 (0-2分)
    dust = spec['dust']
    if dust>=30: sps+=2
    elif dust>=20: sps+=1

    ss['specificity'] = min(35,sps)

    # ============================================================
    # 综合评分
    # ============================================================
    total = max(0, min(100, ss['damage'] + ss['coverage'] + ss['specificity']))
    # 5' C→T < 5%：不满足古DNA损伤基本门槛，扣10分
    if not ss['d5_pass']:
        total = max(0, total - 10)
    # 片段过长直接降级
    long_flag = (mlen > 100)
    if long_flag:
        total = max(0, total - 10)

    if total>=65: conf, lbl = 'HIGH', 'High-confidence ancient DNA'
    elif total>=45: conf, lbl = 'MEDIUM', 'Medium confidence'
    elif total>=30: conf, lbl = 'LOW', 'Low confidence, needs further validation'
    elif total>=12: conf, lbl = 'VERY_LOW', 'Very low confidence, possible false positive'
    else: conf, lbl = 'BACKGROUND', 'Background noise or contamination'
    return total, conf, lbl, ss

def load_abd(path):
    d={}
    if not os.path.exists(path): return d
    with open(path) as f:
        f.readline()
        for l in f:
            p=l.strip().split('\t')
            if len(p)>=5: d[p[0]]={'name':p[1],'reads':int(p[3]) if p[3].isdigit() else 0,'relative_abundance':safe_float(p[4])}
    return d

def load_bd(path):
    d={}
    if not os.path.exists(path): return d
    with open(path) as f:
        f.readline()
        for l in f:
            p=l.strip().split('\t')
            if len(p)>=5: d[p[0]]={'TotalReads':safe_float(p[2]),'UnaggregatedReads':safe_float(p[-2])}
    return d

def load_subs(path):
    d={}
    if not os.path.exists(path): return d
    with open(path) as f:
        for l in f:
            if l.strip():
                ps = parse_subs_line(l)
                if ps: d[ps['taxid']]=ps
    return d

def main():
    ap = argparse.ArgumentParser(description='sedaDNA species authenticity evaluation')
    ap.add_argument('--bamdam_tsv', default='species/compute/bamdam.species.upto_species.tsv')
    ap.add_argument('--bamdam_subs', default='species/compute/bamdam.species.upto_species.subs')
    ap.add_argument('--coverage_dir', default='coverage_analysis')
    ap.add_argument('--abundance_tsv', default='')
    ap.add_argument('--genus_tsv', default='genus/compute/bamdam.genus.upto_genus.tsv')
    ap.add_argument('--family_tsv', default='family/compute/bamdam.family.upto_family.tsv')
    ap.add_argument('--output', default='comprehensive_adna_report.tsv')
    ap.add_argument('--output_detailed', default='comprehensive_adna_detailed.tsv')
    ap.add_argument('--min_reads', type=int, default=50)
    ap.add_argument('--top_n', type=int, default=50)
    args = ap.parse_args()

    print('[Step 1/4] Loading data...', flush=True)
    abd = load_abd(args.abundance_tsv)
    print(f'  Abundance: {len(abd)} species', flush=True)
    gd = load_bd(args.genus_tsv)
    print(f'  Genus data: {len(gd)} genera', flush=True)
    fd = load_bd(args.family_tsv)
    print(f'  Family data: {len(fd)} families', flush=True)
    sd = load_subs(args.bamdam_subs)
    print(f'  Damage profiles: {len(sd)} species', flush=True)

    print('[Step 2/4] Computing metrics...', flush=True)
    if not os.path.exists(args.bamdam_tsv):
        print(f'[ERROR] Missing: {args.bamdam_tsv}', flush=True)
        sys.exit(1)

    # --- 第一遍: 按科(family)分组，计算科内总reads和物种数 ---
    family_groups = defaultdict(lambda: {'total_reads':0,'species_count':0})
    species_rows = []
    with open(args.bamdam_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            taxid = row.get('TaxNodeID','').strip()
            if not taxid: continue
            tr = safe_float(row.get('TotalReads',0))
            if tr < args.min_reads: continue
            species_rows.append(row)
            tp = row.get('taxpath','')
            fid = get_rank_taxid(tp,'family')
            if fid:
                family_groups[fid]['total_reads'] += int(tr)
                family_groups[fid]['species_count'] += 1
    family_groups = dict(family_groups)
    print(f'  Families with species: {len(family_groups)}', flush=True)

    # --- 第二遍: 逐物种计算指标并评分 ---
    results = []
    cnt = 0
    for row in species_rows:
        taxid = row.get('TaxNodeID','').strip()
        if not taxid: continue
        tr = safe_float(row.get('TotalReads',0))
        if tr < args.min_reads: continue

        d5_raw = safe_float(row.get('Damage+1',0))
        d3_raw = safe_float(row.get('Damage-1',0))
        subs = sd.get(taxid, {'5p':{},'3p':{}})
        dmg = damage_metrics(d5_raw, d3_raw, subs)
        mlen = safe_float(row.get('MeanLength',0))
        cov = cov_metrics(taxid, args.coverage_dir)
        spec = spec_metrics(row, gd, fd, family_groups)

        allm = {'damage':dmg,'length':{'mean_length':mlen},'coverage':cov,'specificity':spec,'total_reads':int(tr)}
        total_score, conf, lbl, ss = score(allm)

        results.append({
            'taxid':taxid,'species':row.get('TaxName','').strip('"'),
            'total_reads':int(tr),
            'damage_5p_pct':round(dmg['d5']*100,2),
            'damage_3p_pct':round(dmg['d3']*100,2),
            'damage_asymmetry':round(dmg['asym'],2),
            'damage_decay_slope':round(dmg['slope'],5),
            'damage_decay_r2':round(dmg['r2'],3),
            'damage_signal_5p':round(dmg['signal']*100,2),
            'damage_plateau':round(dmg['plat']*100,2),
            'd5_pass':ss['d5_pass'],
            'mean_length':round(mlen,2),
            'is_long_fragment':(mlen > 100),
            'breadth_1x_pct':round(cov['breadth'],4),
            'mean_depth':round(cov['mdepth'],6),
            'gini_coverage':round(cov['gini'],3),
            'n_contigs':cov['nctg'],'n_contigs_covered':cov['ncov'],
            'contig_covered_pct':round(cov['ctgrat'],2),
            'genome_size_mb':round(cov['gsize']/1e6,2),
            'ani':round(spec['ani'],4),
            'lca_degradation_pct':round(spec['lca'],2),
            'genus_lca_pct':round(spec['glca'],2),
            'family_lca_pct':round(spec['flca'],2),
            'gc_bias':round(spec['gc_bias'],4),
            'uniq_kmer_per_read':round(spec['ukpr'],3),
            'duplicity':round(spec['dup'],2),
            'mean_dust':round(spec['dust'],2),
            'align_ratio':round(spec['align_ratio'],3),
            'read_share':round(spec['read_share'],4),
            'sink_score':round(spec['sink_score'],4),
            'family_reads':spec.get('family_total',0),
            'family_species':spec.get('family_species_count',0),
            'score_damage':ss['damage'],'score_length_flag':ss['length'],
            'score_coverage':ss['coverage'],'score_specificity':ss['specificity'],
            'total_score':total_score,'confidence':conf,'diagnosis':lbl,
            'has_cov':cov['has'],
        })
        cnt+=1
        if cnt%10==0: print(f'  ... processed {cnt} species', flush=True)

    print(f'  Total: {len(results)} species (min_reads >= {args.min_reads})', flush=True)

    print('[Step 3/4] Writing reports...', flush=True)
    results.sort(key=lambda x: x['total_score'], reverse=True)

    sf = ['species','total_reads','damage_5p_pct','damage_3p_pct','damage_asymmetry',
          'damage_decay_r2','damage_signal_5p','d5_pass','mean_length','is_long_fragment',
          'breadth_1x_pct','gini_coverage','contig_covered_pct',
          'ani','lca_degradation_pct','gc_bias','uniq_kmer_per_read','duplicity','mean_dust','align_ratio',
          'read_share','sink_score','family_reads','family_species',
          'score_damage','score_coverage','score_specificity',
          'total_score','confidence','diagnosis']
    sh = ['Species','TotalReads','Damage5p%','Damage3p%','DmgAsymmetry',
          'DmgDecayR2','DmgSignal5p%','D5_Pass(>=5%)','MeanLen','LongFrag(>100bp)',
          'Breadth1x%','GiniCov','ContigCov%',
          'ANI','LCADeg%','GCbias','KmerPerRead','Duplicity','Dust','AlignRatio',
          'ReadShare','SinkScore','FamReads','FamSpp',
          'S_Damage(35)','S_Cov+Reads(30)','S_Specificity(35)',
          'Total(100)','Confidence','Diagnosis']

    with open(args.output,'w') as out:
        out.write('\t'.join(sh)+'\n')
        for r in results[:args.top_n]:
            out.write('\t'.join(str(r.get(f,'')) for f in sf)+'\n')
    print(f'  Simple report: {args.output} ({min(args.top_n,len(results))} rows)', flush=True)

    df = ['taxid','species','total_reads','damage_5p_pct','damage_3p_pct','damage_asymmetry',
          'damage_decay_slope','damage_decay_r2','damage_signal_5p','damage_plateau',
          'd5_pass','mean_length','is_long_fragment','breadth_1x_pct','mean_depth','gini_coverage',
          'n_contigs','n_contigs_covered','contig_covered_pct','genome_size_mb',
          'ani','lca_degradation_pct','genus_lca_pct','family_lca_pct',
          'gc_bias','uniq_kmer_per_read','duplicity','mean_dust','align_ratio',
          'read_share','sink_score','family_reads','family_species',
          'score_damage','score_length_flag','score_coverage','score_specificity',
          'total_score','confidence','diagnosis']
    with open(args.output_detailed,'w') as out:
        out.write('\t'.join(df)+'\n')
        for r in results:
            out.write('\t'.join(str(r.get(f,'')) for f in df)+'\n')
    print(f'  Detailed report: {args.output_detailed} ({len(results)} rows)', flush=True)

    print('[Step 4/4] Summary', flush=True)
    print('='*70, flush=True)
    high=sum(1 for r in results if r['confidence']=='HIGH')
    med=sum(1 for r in results if r['confidence']=='MEDIUM')
    low=sum(1 for r in results if r['confidence']=='LOW')
    vlow=sum(1 for r in results if r['confidence']=='VERY_LOW')
    bg=sum(1 for r in results if r['confidence']=='BACKGROUND')
    print(f'  Total evaluated:     {len(results)}', flush=True)
    print(f'  HIGH confidence:     {high}', flush=True)
    print(f'  MEDIUM confidence:   {med}', flush=True)
    print(f'  LOW confidence:      {low}', flush=True)
    print(f'  VERY_LOW confidence: {vlow}', flush=True)
    print(f'  BACKGROUND noise:    {bg}', flush=True)
    print('='*70, flush=True)
    print('\n  >>> Top 15 Species <<<', flush=True)
    for r in results[:15]:
        print(f"  {r['total_score']:4d}  {r['species']:<35s}  D5={r['damage_5p_pct']:5.1f}%  R2={r['damage_decay_r2']:.2f}  Br={r['breadth_1x_pct']:7.3f}%  Gini={r['gini_coverage']:.2f}  ANI={r['ani']:.3f}  LCA={r['lca_degradation_pct']:5.1f}%  {r['diagnosis']}", flush=True)
    print(f'\n[Done] Reports: {args.output}, {args.output_detailed}', flush=True)

if __name__ == '__main__':
    main()
