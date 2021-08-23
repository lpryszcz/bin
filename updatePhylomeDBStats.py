#!/usr/bin/env python2
desc="""Generate phylomeDB stats file. Check no. of phylomes, trees, algs, proteins,
species and genomes. In addition, checks for web-serwer stats from last week i.e.
the most common query etc.
NOTE: Need to be run from cgenomics.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona/Mizerow, 25/06/2014
"""

import os, sys, time
import locale
locale.setlocale(locale.LC_ALL, 'en_US.utf8')
import MySQLdb
from datetime import datetime
#phylomedb wsgi
sys.path.append('/home/services/web/phylomedb.org/wsgi')
from phylomedb import PUBLIC_PHYLOMES

def phylome_stats( ):
    """Generate database stats."""
    conn = MySQLdb.connect(host="mysqlsrv-tgabaldon.linux.crg.es", user="phyReader", \
                           passwd="phyd10.-Reader", db='phylomedb')
    c = conn.cursor()
    
    #public phylomes are stored by python, not by mysql
    public = len(PUBLIC_PHYLOMES) #19 
    
    #trees
    c.execute("SELECT COUNT(*) FROM tree" ) 
    trees, = c.fetchone()
    
    #ml trees
    c.execute("SELECT COUNT(*) FROM tree WHERE method != 'NJ'")
    mltree,= c.fetchone()

    #algs
    c.execute("SELECT COUNT(*) FROM alignment")
    algs,  = c.fetchone()
    
    #proteins
    c.execute("SELECT COUNT(*) FROM unique_protein")
    prots, = c.fetchone()
    
    #species
    c.execute("SELECT COUNT(*) FROM species")
    specs, = c.fetchone()

    #genomes
    c.execute("SELECT COUNT(*) FROM genome")
    genms, = c.fetchone()

    
    lines  = """<p><strong>Public Phylomes: </strong><a href="?q=phylomes">%s</a></p>
<p><strong>Trees: </strong> %s (%s Max.Likelihood)</p>
<p><strong>Alignments: </strong> %s</p>
<p><strong>Proteins: </strong> %s</p>
<p><strong>Species: </strong> %s</p>
<p><strong>Genomes: </strong> %s</p>
""" % (public, locale.format("%d", trees, grouping=True), locale.format("%d", mltree, grouping=True), \
       locale.format("%d", algs, grouping=True), locale.format("%d", prots, grouping=True), \
       locale.format("%d", specs, grouping=True), locale.format("%d", genms, grouping=True) )
    
    #x-links
    lines += "<p><strong>External IDs: </strong>\n<ul>"
    c.execute("SELECT external_db, count(*) FROM external_id GROUP BY external_db")
    line = "\n <li><strong>%s</strong>: %s [%.2f%s]</li>"
    for db, c in c.fetchall():
        lines += line%(db, locale.format("%d", c, grouping=True), c*100.0/prots, '%')

    lines += "</ul></p>"
    return lines

def webserver_stats(serverlog, days=31):
    """Generate webserver stats i.e. most popular query."""
    query=protein=tree=phylome = ""
    phylomes = {}
    proteins = {}
    trees    = {}
    algs     = {}
    d2       = datetime.now()
    #parse log
    for l in open(serverlog):
        if ' "GET /?q=' not in l:
            continue

        #ip,v1,v2,date,request,status,reply_size,referer,user_agent
        ip,h1,h2,date,timezone,get,request,request_type,status,rSize,referer = l.split()[:11]

        if status!='200':
            continue
        
        #look for entries at mosth one month old
        date = date.strip('[')
        t1   = time.strptime( date,'%d/%b/%Y:%H:%M:%S' )
        d1   = datetime(*t1[:6])
        delta= d2-d1
        if delta.days > days:
            continue

        #/?q=phylome_browser&phyid=26
        #/?q=search_tree&seqid=Phy0002IQG_CANAL&phyid=23&method=best
        #/?q=search_alg&seqid=Phy00086SJ_HUMAN&phyid=1&alg_type=clean "http://www.phylomedb.org/?q=search_tree&seqid=TP53"
        phyid=seqid=method=algtype=""
        if "phyid=" in request:
            phyid = request.split("phyid=")[1]
            phyid = phyid.split()[0]
            phyid = phyid.split("&")[0]
            phyid = phyid.split("?")[0]
            try: 
                phyid = int(phyid)
            except: 
                print "Error: Wrong phyid: %s" % phyid
                continue
            #only public
            if phyid not in PUBLIC_PHYLOMES:
                continue
            if phyid not in phylomes:
                phylomes[phyid] = 0
            phylomes[phyid] += 1

        if "seqid=" in request:
            seqid = request.split("seqid=")[1]
            if not seqid or seqid=="RANDOM": continue
            seqid = seqid.split()[0]
            seqid = seqid.split("&")[0]
            if seqid not in proteins:
                proteins[seqid] = 0
            proteins[seqid] += 1

        if "method=" in request:
            method= request.split("method=")[1]
            method= method.split()[0]
            method= method.split("&")[0]

        if "alg_type=" in request:
            algtype=request.split("alg_type=")[1]
            algtype=algtype.split()[0]
            algtype=algtype.split("&")[0]

        if "search_tree" in request and method and phyid:
            tree = "%s|%s|%s" % ( seqid,phyid,method )
            if tree not in trees:
                trees[tree] = 0
            trees[tree] += 1

        elif "search_alg" in request:
            alg  = "%s|%s|%s" %( seqid,phyid,algtype )
            if alg not in algs:
                algs[alg] = 0
            algs[alg] += 1

    #get most popular
    for ph in phylomes.keys():
        if not ph in PUBLIC_PHYLOMES:
            del phylomes[ph]
    phylome = sorted(phylomes.keys(), key=lambda x: phylomes[x], reverse=True )[0]
    tree    = sorted(trees.keys(),    key=lambda x: trees[x],    reverse=True )[0]
    protein = sorted(proteins.keys(), key=lambda x: proteins[x], reverse=True )[0]
    alg     = sorted(algs.keys(),     key=lambda x: algs[x],     reverse=True )[0]


## <li><strong>tree:           </strong> <a href="/?q=search_tree&seqid=%s&phyid=%s&method=%s">%s (%s)</a></li>
## tree.split('|')[0], tree.split('|')[1], tree.split('|')[2]
    lines = """
<p><strong>The most popular (last %s days): </strong>
<ul>
<li><strong>public phylome: </strong> <a href="/phylome_%s">%s (%s)</a></li>
<li><strong>tree:           </strong> <a href="/?q=search_tree&seqid=%s&phyid=%s">%s (%s)</a></li>
<li><strong>alignment:      </strong> <a href="/?q=search_alg&seqid=%s&phyid=%s&alg_type=%s">%s (%s)</a></li>
<li><strong>protein:        </strong> <a href="/?q=seqinfo&seqid=%s">%s (%s)</a></li>
</ul>
</p>

""" % (days, phylome, phylome, phylomes[phylome], \
       tree.split('|')[0], tree.split('|')[1], \
       tree.split('|')[0], trees[tree], \
       alg.split('|')[0], alg.split('|')[1], alg.split('|')[2], \
       alg.split('|')[0], algs[alg], \
       protein, protein, proteins[protein])
    return lines

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")
    parser.add_argument('--version', action='version', version='1.1')   
    parser.add_argument("-o", "--outfile", 
                       default="/home/services/web/phylomedb.org/tmp/phylome.stats.txt", 
                       help="output file [%(default)s]")
    parser.add_argument("-l", "--serverlog",
                       default="/home/services/web/logs/phylomedb.org_access.log",
                       help="output file [%(default)s]")
    parser.add_argument("-d", "--days", default=30, type=int,
                       help="how many days backwards to look [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    
    lines = ""

    #generate database stats
    lines += phylome_stats()
    
    #generate last week access stats
    lines += webserver_stats(o.serverlog, o.days)
    
    #add current date
    date   = datetime.ctime(datetime.now())
    lines += '<p style="text-align: right;"><span style="font-size: smaller;"><em>(%s)</em></span></p>\n'%date

    #save to file
    #print lines
    out = open(o.outfile, 'w')
    out.write(lines)
    out.close()

if __name__ == "__main__":
    t0=datetime.now()
    main()
    dt=datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
