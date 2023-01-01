import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os, re, sys

# configurations
configs = {
    'num_processors': 1,
    'genome_dir': './genomes',
    'models': './resources/gwas_filtered_fun.csv.gz',
    'background': '',
    'prefix': 'out', # output prefix
    'data_type': '23andme' # or 'vcf'
}


def match_23andme(genome, snps, chunksize=10**5):
    """
    match 23andme data to snps, return matches
    process 100k records at a time to save memory
    matches is indexed as snp:effect_allele, has a single column [0,1,2] indicting number of matches
    """
    m = []
    for chunk in pd.read_csv(genome, sep='\t', usecols=[0,3], names = ['SNPS','geno'], index_col=0, comment='#', low_memory=False, chunksize=chunksize):
        m.append(snps.join(chunk, how='inner'))
    m = pd.concat(m, axis=0)    
    m['geno'] = (m.geno.str[0] == m[1]).astype(int) + (m.geno.str[1] == m[1]).astype(int)
    m[1] = m.index + ':' + m[1]
    m = m.set_index(1)
    return m

def map_alleles(x):
    """
    map gt to alleles
    """
    alleles = [x['ref']] + x['alt'].split(',')
    gt = re.split(r'[|/:]', x['gt'])
    if len(gt)>1:
        a1, a2 = gt[0:2]
        return alleles[int(a1)], alleles[int(a2)]
    else: # assuming the missing allele is reference
        return alleles[int(gt[0])], alleles[0]    

def match_vcf(genome, snps, chunksize=10**5):
    """
    match vcf data to snps, return matches
    process 100k records at a time to save memory
    matches is indexed as snp:effect_allele, has a single column [0,1,2] indicting number of matches
    """
    m = []
    cols = [0,1,2,3,4,9]
    names = ['chr', 'position', 'rsid', 'ref', 'alt', 'gt']
    for chunk in pd.read_csv(genome, sep='\t', header=None, usecols=cols, names=names, index_col=2, comment='#', low_memory=False, chunksize=chunksize):
        # remove null values in gt
        # chunk = chunk['.' not in chunk['gt']]
        chunk['a1'], chunk['a2'] = zip(*chunk[['ref', 'alt', 'gt']].apply(map_alleles, axis=1))
        m.append(snps.join(chunk[['a1', 'a2']], how='inner'))
    m = pd.concat(m, axis=0)
    m['geno'] = (m[1].str[0] == m['a1']).astype(int) + (m[1].str[0] == m['a2']).astype(int)
    m[1] = m.index + ':' + m[1]
    m = m.set_index(1)
    return pd.DataFrame(m['geno'])

def parse_ref():
    """          
    parse EBI filtered GWAS catalog, return traits and SNPs for matching
    """
    ref = pd.read_csv(configs['models'], sep='\t', compression='gzip', dtype={'PUBMEDID':str})
    traits = []
    for trait, variants in ref.groupby('DISEASE/TRAIT'):
        traits.append( {
            "name": trait,
            "variants": variants['variant'],
            "OR": variants['OR or BETA'].values,
            "literature": set(variants['PUBMEDID']),
            "dates": set(variants['DATE ADDED TO CATALOG'])
        })
    # traits = pd.DataFrame.from_dict(traits) # dataframe version
    snps = ref['variant'].str.split(':', expand=True)[[0,1]]
    snps = snps.drop_duplicates()
    snps = snps.set_index(0)
    
    return traits, snps

def matchOne(genome1, name1):
    """
    match one 23andme or vcf genomes, return matches
    """
    traits, snps = parse_ref()
    if configs['data_type'] == '23andme':
        m1 = match_23andme(genome1,snps)
    elif configs['data_type'] == 'vcf':
        m1 = match_vcf(genome1,snps)
    else:
        print("Unsupported Data type")
        sys.exit(2)
    matched = {}

    for trait in traits:
        v = trait["variants"] # extra row
        v = v.reindex(index=v.values)
        # linear combination
        aligned = v.align(m1, level='variant', join='left')[1]['geno'].values
        missing_sites = sum(np.isnan(aligned))
        aligned[np.isnan(aligned)] = 0
        v1 = np.dot(aligned, trait['OR'])
        matched[trait['name']] = [v1, ','.join(list(trait['dates'])), ','.join(list(trait['literature'])), len(aligned), missing_sites]   
    matched = pd.DataFrame.from_dict(matched, orient='index') 
    matched.fillna(0, inplace=True)
    matched = matched.rename(columns={0:name1, 1:'Publication Dates', 2:'Pubmed IDs', 3:'#sites_required', 4:'#sites_absent'})
    if configs['background']:
        percentiles = pd.read_csv(configs['background'], sep='\t', compression='gzip', index_col=0)
        matched = matched.join(percentiles, how='left')
        matched = matched[matched['max']>0]    
    return matched

def boxplot(df, ax=None, box_width=0.2, whisker_size=20, mean_size=10, median_size = 10 , line_width=1.5, xoffset=0, color=[0,0,0]):
    """Plots a boxplot from existing percentiles.

    Parameters
    ----------
    df: pandas DataFrame
    ax: pandas AxesSubplot
        if to plot on en existing axes
    box_width: float
    whisker_size: float
        size of the bar at the end of each whisker
    mean_size: float
        size of the mean symbol
    color: int or rgb(list)
        If int particular color of property cycler is taken. Example of rgb: [1,0,0] (red)

    Returns
    -------
    f, a, boxes, vlines, whisker_tips, mean, median
    """

    if type(color) == int:
        color = plt.rcParams['axes.prop_cycle'].by_key()['color'][color]

    if ax:
        a = ax
        f = a.get_figure()
    else:
        f, a = plt.subplots()

    boxes = []
    vlines = []
    xn = []
    # legends = df.index
    df = df.reset_index()
    df = df.drop('index', axis=1)
    
    for row in df.iterrows():
        x = row[0] + xoffset
        xn.append(x)

        # box
        y = row[1][25]
        height = row[1][75] - row[1][25]
        box = plt.Rectangle((x - box_width / 2, y), box_width, height)
        a.add_patch(box)
        boxes.append(box)

        # whiskers
        y = (row[1][95] + row[1][5]) / 2
        vl = a.vlines(x, row[1][5], row[1][95])
        vlines.append(vl)

    for b in boxes:
        b.set_linewidth(line_width)
        b.set_facecolor([1, 1, 1, 1])
        b.set_edgecolor(color)
        b.set_zorder(2)

    for vl in vlines:
        vl.set_color(color)
        vl.set_linewidth(line_width)
        vl.set_zorder(1)

    whisker_tips = []
    if whisker_size:
        g, = a.plot(xn, df[5], ls='')
        whisker_tips.append(g)

        g, = a.plot(xn, df[95], ls='')
        whisker_tips.append(g)

    for wt in whisker_tips:
        wt.set_markeredgewidth(line_width)
        wt.set_color(color)
        wt.set_markersize(whisker_size)
        wt.set_marker('_')

    mean = None
    if mean_size:
        g, = a.plot(xn, df['mean'], ls='')
        g.set_marker('.')
        g.set_markersize(mean_size)
        g.set_zorder(20)
        g.set_markerfacecolor('None')
        g.set_markeredgewidth(line_width)
        g.set_markeredgecolor(color)
        mean = g

    median = None
    if median_size:
        g, = a.plot(xn, df['median'], ls='')
        g.set_marker('_')
        g.set_markersize(median_size)
        g.set_zorder(20)
        g.set_markeredgewidth(line_width)
        g.set_markeredgecolor(color)
        median = g
        
    # userdata
    g, = a.plot(xn, df['user'], ls='')
    g.set_marker('o')
    g.set_markersize(median_size)
    g.set_zorder(20)
    g.set_markeredgewidth(line_width)
    g.set_markeredgecolor('red')
        
    a.set_ylim(np.nanmin(df)-1, np.nanmax(df)+1)
    plt.xticks(ticks=np.arange(len(df)))
    # a.set_xticklabels(list(legends), rotation = 90, ha="right")
    return f

def format_match_table(matched):
    """
    format match table, adding abstract
    matched has three columns, user_score, publication dates, PubMedIDs
    """
    tables = []
    i = 0
    for row in matched.iterrows():
        html_code = '<table border=1 width=700>'
        html_code += '<tr><th width=200 align=left>' + str(i) + '. <td align=right><font color="green">' + row[0] + '</font></td></tr>'
        html_code += '<tr><th width=200 align=left>Your Score: <td align=right>' + '{:.2f}'.format(row[1][0]) + '</td></tr>'
        html_code += '<tr><th width=200 align=left># sites required: '+ str(row[1][3]) 
        html_code += '<td align=right>Your genome data misses ' + str(row[1][4]) + ' required sites.</td></tr>'
        html_code += '<tr><th width=200 align=left>Publication Date:<td align=right>' + row[1][1] + '</td></tr>'
        pubmed = ''
        for l in row[1][2].split(','):
            link = 'https://www.ncbi.nlm.nih.gov/pubmed/?term=' + l
            pubmed += '<a href="' + link + '" target="_blank">' + l + '</a>, '
        html_code += '<tr><th width=200 align=left>Links to publication<td align=right>' + pubmed + '</td></tr>'
        html_code += '</table>'
        tables.append(html_code)
        i+=1

    return tables


