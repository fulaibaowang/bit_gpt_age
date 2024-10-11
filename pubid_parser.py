# ---
# jupyter:
#   jupytext:
#     cell_markers: '{{{,}}}'
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.6
#   kernelspec:
#     display_name: Python 3.10.8
#     language: python
#     name: py3.10.8
# ---

import pubstats.__init__ as ps
import requests
import pandas as pd
import numpy as np
import json
import difflib
import time
import sys
import os
from bs4 import BeautifulSoup
# from py2cytoscape import cyrest
import subprocess as sb
from subprocess import Popen, PIPE, STDOUT
import matplotlib
import matplotlib.pyplot as plt
import paramiko
import html

# {{{
pmc="/nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/pmc"
pbmc_folders=os.listdir(pmc)
pbmc_folders=[ s for s in pbmc_folders if os.path.isdir(f"{pmc}/{s}") ] 
pmc_files=pd.DataFrame()
for f in  pbmc_folders:
    print(f)
    sys.stdout.flush()
    files=os.listdir(f"{pmc}/{f}")
    pmc_files_=pd.DataFrame()
    for f_ in files:
        content=[]
        with open(f"{pmc}/{f}/{f_}", "r" ) as txt :
            for l in txt:
                l=l.split("\t")[0]
                content.append(l)
        target=f_.replace("filelist.txt","tar.gz") 
        if f == "manuscript":
            target=f"https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/{target}"
        else: 
            target=f"https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/{f}/txt/{target}"
        tmp=pd.DataFrame({"folder":f,"file":f_, "files":content, "tar.gz": target  })
        pmc_files_=pd.concat([pmc_files_,tmp])
    pmc_files_.to_csv(f"{pmc}/{f}.tsv", sep="\t",index=None)
    pmc_files=pd.concat([pmc_files,pmc_files_])
    
pmc_files=pmc_files[pmc_files["files"]!="Article File"]
pmc_files["subfolder"]=pmc_files["files"].apply(lambda x: x.split("/")[0])
def get_ftp_pmc(x):
    try:
        r=x.split("/")[1].split(".txt")[0] 
    except:
        print(x)
        r=x.split("/")[1].split(".txt")[0] 
    return r

pmc_files["pmc"]=pmc_files["files"].apply(lambda x: get_ftp_pmc(x) )
pmc_files.head()
# }}}

tmp=pmc_files[pmc_files["pmc"]=="PMC6225988"]
print(tmp["tar.gz"].tolist())
tmp

# {{{
BASEURL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ESEARCH=BASEURL+"esearch.fcgi?db=pubmed&rettype=json&retmode=json&term="
ESUMMARY=BASEURL+"esummary.fcgi?db=pubmed&rettype=json&retmode=json&id="

title="Sequencing the Immunoglobulin Heavy-Chain Locus (IgH) in Turquoise Killifish ( Nothobranchius furzeri)"

def get_pubmedid(title,author,doi,ESEARCH=ESEARCH ):
    
    if doi:
        url=ESEARCH+doi#+"[Title]"
        r = requests.post(url = url, verify=False)
        r = json.loads(r.content)
        data=r['esearchresult']
        idlist=data['idlist']
        if len(idlist) == 1:
            return idlist[0]
        
    url=ESEARCH+'"'+title+'"'#+"[Title]"
    r = requests.post(url = url, verify=False)
    r = json.loads(r.content)
    data=r['esearchresult']
    idlist=data['idlist']
    if len(idlist) != 1:
        url=url+" "+author
        r = requests.post(url = url, verify=False)
        r = json.loads(r.content)
        data=r['esearchresult']
        idlist=data['idlist']
        if len(idlist) != 1:
            print(len(idlist), doi, title)
            return None
        else:
            return idlist[0]

    else:
        return idlist[0]
        
    
# print(idlist)
# for article_id in idlist:
#     url = esummary+article_id
#     r = requests.post(url = url, verify=False)
#     r = json.loads(r.content)
# }}}

get_pubmedid("ncreasing number of long-lived ancestors marks a decade of healthspan extension and healthier metabolomics profiles", "van den Berg, Niels", "10.1038/s41467-023-40245-6")

[ "Deletion of endogenous Tau proteins is not detrimental in Drosophila", "26976084", "10.1038/srep23102"],\
[ "Dietary Restriction: Theory Fails to Satiate - Response",

# {{{
# pubmedids=[]
# pubmedid=None
# title=None
# author=None
# titleid=None
# doi=None
# with open("export_age_all_publ_20230817.txt", "r") as f:
#     for l in f.readlines():
#         if ( l[:2] == "%A" ) and ( not author ):
#             author=l[3:].split("\n")[0] #-4]
#         if l[:9] == '%F OTHER:':
#             pubmedid=l.split("OTHER: ")[1].split("\n")[0]
#         if l[:2] == "%T" :
#             title=l[3:].split("\n")[0]   #-4]
#         if l[:2] == "%R" :
#             doi=l[3:].split("\n")[0]   
#         if ( l == "\n" ) and title:
#             # print(pubmedid,titleid,doi, title )
#             if not pubmedid:
#                 pubmedid=get_pubmedid(title, author,doi)
#             pubmedids.append([pubmedid,doi, title ])
#             pubmedid=None
#             title=None
#             titleid=None
#             author=None    
#             doi=None
# }}}

# {{{
# pubmedids=pd.DataFrame(pubmedids,columns=["PMID","DOI","Title"])
# pubmedids.columns=["PMID","DOI","Title"]
# pubmedids.head()
# pubmedids.to_csv("pubmedids.tsv", sep="\t",index=None)
# }}}

pubmedids=pd.read_csv("pubmedids.tsv",sep="\t")
pubmedids.head()

ids=pubmedids.dropna(subset=["PMID"])
noids=pubmedids[ ~pubmedids.index.isin(ids.index.tolist() ) ]
print(len(noids))
ids_values=ids["PMID"].tolist()
noids

pmc_files.head()

# {{{
# article_id=ids_values[0]
# print(article_id)
article_ids=[]
for article_id in ids_values:
    # works=False
    # while not works:
    try:
        url=ESUMMARY+article_id
        r = requests.post(url = url, verify=False)
        r = json.loads(r.content)
        data=r['result'][article_id]
        values=data["articleids"]
        pmc=[ s["value"] for s in values if s["idtype"]=="pmc" ]
        pmcid=[ s["value"] for s in values if s["idtype"]=="pmcid" ]
        if len(pmc) == 1 :
            r=pmc[0]
        elif len(pmcid) == 1:
            r=pmc[0]
        else:
            r=None
        article_ids.append([ article_id, r] )
    except:
        print("failed for ", url)
        # time.sleep(2)
        sys.stdout.flush()
        
article_ids=pd.DataFrame(article_ids, columns=["PMID","PMCID"])
# }}}

pubmedids=pd.merge(pubmedids,article_ids,on=["PMID"],how="left")
pubmedids.to_csv("pubmedids.tsv", sep="\t",index=None)

pubmedids=pd.read_csv("pubmedids.tsv",sep="\t")
pubmedids.head()

# {{{
bash='''
#!/bin/bash

cd /nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/
mkdir -p pmc_papers
cd pmc_papers

#################


'''

pubmedids_valid=pubmedids[["PMCID"]].dropna()["PMCID"].tolist()

papers=pmc_files.dropna(subset=["tar.gz"])
papers=papers[papers["pmc"].isin(pubmedids_valid)]
print(len(papers))

for target in list( set(papers["tar.gz"].tolist() )) :
    if target:
        wget=f"wget {target}\n"

        bash=bash+wget

        files=papers[papers["tar.gz"]==target]["files"].tolist()
        files=list(set(files))
        files=" ".join(files)
        targz=target.split("/")[-1]
        untar=f"tar -zxvf {targz} {files}\n"

        bash=bash+untar

        rm=f"rm -rf {targz}\n"

        ash="\n#################\n\n"

        bash=bash+rm+ash

print(bash)

with open("download.extract.pmc.sh", "w" ) as sh:
    sh.write(bash)
        
# }}}

pmc_files.head()


