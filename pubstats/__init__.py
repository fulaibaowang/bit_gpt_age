import sys
import requests
import pandas as pd
import numpy as np
import json
import difflib
import time
from bs4 import BeautifulSoup
# from py2cytoscape import cyrest
import subprocess as sb
from subprocess import Popen, PIPE, STDOUT
import matplotlib
import matplotlib.pyplot as plt
import paramiko

import warnings
warnings.filterwarnings('ignore')




baseurl="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
esearch=baseurl+"esearch.fcgi?db=pubmed&rettype=json&retmode=json&term="
esummary=baseurl+"esummary.fcgi?db=pubmed&rettype=json&retmode=json&id="

#bad_names={"Lpez-Otn C":'Lopez-Otin C',"Cochem HM":"Cocheme HM",\
#          "Grnke S":"Gronke S","Sgr CM":"Sgro CM","Khl I":'Kuhl I',\
#          "Prez-Prez R":'Perez-Perez R',"Sphr H":'Spahr H',"Hllberg BM":'Hallberg BM',\
#          "Hagstrm E":'Hagstrom E',"Cmara Y":'Camara Y',"Srensen L":'Sorensen L'}


# list_of_last_authors=["Larsson NG","Partridge L","Antebi A","Langer T","Demetriades C",\
#                      "Graef M","Tessarz P","Valenzano DR","Wickstrom SA","Denzel MS",\
#                      "Stewart JB","Matic I"]
# files_tag="MPI-AGE"
# institute_sub_strings=[ "planck institute for biology of ag", "planck institute for the biology of ag"]
# # we use the variable kolle to keep a list of substrings which we will use 
# # for identifying safe authors. ie. if any affiliation of any of the authors
# # contains one of these substrings we will define all authors on the paper
# # as safe authors. based on these safe authors we then expande a social network
# # which we define as safe.
# kolle=[ "planck institute for biology of ag", "planck institute for the biology of ag", "koeln", "koln", "kln", "koeln", "cologne", "cecad", "ageing", "aging"]


#list_of_last_authors=["Hoehn M","Backes H", "Graf R", "Steculorum SM", "Korotkova T", "Bruning JC",\
#                     "Kornfeld JW", "Wunderlich TF", "Fenselau H", "Tittgemeyer M"]
#files_tag="MPI-MET"
#institute_sub_strings=[ "cecad", "institute for neurological research", "max planck institute for metabolism research", "cellular stress #responses", "cologne" ]
#kolle=[ "planck institute for biology of ag", "institute for neurological research", "planck institute for the biology of ag", "koeln", #"koln", "kln", "koeln", "cologne", "cecad", "ageing", "aging", "cecad", "max planck institute for metabolism research", "cellular stress #responses" ]


def get_aff(x):
    """
    Get affiliations from text like output
    """
    pos=[]
    for s, i in zip( x, range( len(x) ) ):
        if s == "(":
            p=x[i:].split(")")[0]+")"
            pos.append(p)
    return pos

#def fix_non_ascii(x):
#    """
#    Fixes non ascii strings
#    """
#    
#    try:
#        res=str(x)
#    except:
#        res=x.encode('ascii','ignore')
#        res=str(res)
#    return res

def fix_non_ascii(x):
    """
    Function for fixing non ascii text
    """ 

    try:
        res=str(x)
    except:
        res=x.encode('ascii','ignore')
        res=str(res)
    res=res.replace("\xfc","u")
    return res

def clean_encoding(x):
    """
    Cleans utf8 encoding
    """
    
    #x=str(x)
    try:
        x=x.decode('utf8')
    except:
        x=x
    try:
        x=x.encode('ascii','ignore')
    except:
        x=str(x).encode('ascii','ignore')
    x=str(x)
    return x

def get_affiliations_json(article_id):
    """
    Gets the abstract from a paper in bad json format 
    and collects the affiliations of the authors.
    
    :param article_id: an article pubmed id
    :returns: a dictionary of the form { author : affiliation }
    
    """
    
    efetch=baseurl+"efetch.fcgi?db=pubmed&id="
    url=efetch+article_id

    r = requests.post(url = url, verify=False)
    
    affiliations={}
    
    text=str(r.content)
    text_=text.split(" ::= ")[-1]
    if "names std {\n" in text_:
        text=text_.split("names std {\n")[-1].split("from journal {\n")[0]
        text=text.replace("\n","")
        text=text.split("name ml")
        text=text[1:]
        for l in text:
            author=l.split(",")[0][2:-1]
            if '"' in author:
                author=author.split('"')[0]
            if 'affil str' in l:
                affiliation=l.split('affil str "')[1].split('"')[0]
            else:
                affiliation=str(None)
            affiliations[author]=affiliation

    if "names ml {" in text_:
        text=text_.split("names ml {\n")[1].split("from journal {\n")[0].split("{")[0]
        text=text.replace("\n","")
        text=text.replace('"','').split("}")[0].rstrip(" ")
        text=text.split(",")
        
        for l in text:
            author=l.lstrip(" ")
            affiliations[author]=str(None)
            
    return affiliations

def clean_address(x):
    """
    Given an affiliation field it cleans up and extracts the address / affiliation
    """
    add=None
    inst=None
    x_=[]
    bea=BeautifulSoup(x)

    # the affiliation might be under <institution>
    # or under <addr-line>        
    if "<institution>" in x:
        inst=bea.findAll('institution')
        inst=[ clean_encoding(s.contents[0]) for s in inst ]
        inst="; ".join(inst)
        x_.append(inst)
    if "<addr-line>" in x:
        add=bea.findAll('addr-line')
        add=[ clean_encoding(s.contents[0]) for s in add ]
        add="; ".join(add)
        x_.append(add)

    # if we fail to clean/extract the affiliation 
    # we simply make sure we have ascii endoding of
    # the given input
    if len(x_)>0:
        res=", ".join(x_)
    else:
        res=clean_encoding(x)
    return res

#def pmc_affiliations(pmc_id="PMC1682176",known_affiliations=affiliations, tmpid=None):
def pmc_affiliations(pmc_id,known_affiliations, tmpid=None):
    """
    Retrieves afiliations from a pmc available paper.
    
    :param pmc_id: a pmc id
    :param known_affiliations: a dictionary of the form { author:[afilliations] }
    :param tmpid: used for verbose printing of id and target URL
    
    :returns: a dictionary of the form { author : [affiliations] }
    """
    
    FAILURE=False
    
    # get already known list of authors from known_affiliations argument
    list_of_authors=known_affiliations.keys()
    
    # request URL
    # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=1682176
    pmc_fetch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id="
    pmc_id=pmc_id.upper().lstrip("PMC")
    url=pmc_fetch+pmc_id
    
    #if tmpid:
    #    print tmpid
    #    print url
    
    # request
    r = requests.post(url = url)
    cont = clean_encoding(r.content)
    bs = BeautifulSoup(cont)
    
    # get contributors
    contributors=bs.article.front.findAll("contrib")
    
    authors={} # a dictionary of the form { author : { rid : [ list, of, indexes ] } }
    for c in contributors:
        # define each author by last name and first name
        try:
            first_name=c.findAll("given-names")[0].contents[0]
            last_name=c.surname.contents[0]
        except:
            #print c.findAll("given-names")
            #print c.findAll("given-names")[0].contents
            #print c.surname.contents[0]
            print(pmc_id, "Failiing to properly filter given names")
            first_name=c.findAll("given-names")[0].contents[0]
            last_name=c.surname.contents[0]        
        
        author=last_name+" "+first_name[0]
         
        # match the defined name by matching it to the known list of authors
        if author not in list_of_authors:
            possibilitites=difflib.get_close_matches(author,list_of_authors)
            if len(possibilitites) > 0:
                author=possibilitites[0]
            #else:
            #    FAILURE=True
            #    print "Coult not find ", author, "in\n", possibilitites
                
        # get the affiliations indexes for the current contributor
        affs=c.findAll("xref")
        
        affs_={} # a dictionary of the form { rid : [ list, of, indexes ]}
        for a in affs:
            rid=str(a["rid"]).upper()
            sups=a.findAll("sup")
            if len(sups) > 0:
                sups=[  s.contents[0] for s in sups ]
            else:
                try:
                    sups=[ a.contents[0] ]
                except:
                    sups=[]
            affs_[rid]=sups
        
        #affs=[ s["rid"].upper() for s in affs ]
        #affs=c.findAll("sup")
        #print affs
        #affs=[ s.contents[0] for s in affs ]
        #print affs
        #print author, affs_
        authors[author]=affs_           
    
    # get affiliations
    affiliations=bs.article.front.findAll("aff")
    affs={}
    for c in affiliations:

        # we either have an affiliation / index for each id (try) or
        # a single affiliation / index for all authors (except)
        try:
            rid=c["id"].upper()#.findAll("label")[0].contents[0]
        except:
            rid="all".upper()
        
        # we either have the affiliations under "sup"
        # or under "label"
        # or "insitution"
        sups=c.findAll("sup")
        labels=c.findAll("label")

        
        if len(sups) > 0:
            # if affiliatons under sup we extract the affiliations here
            sups=str(c).split("</aff>")[0].split("<sup>")[1:] # list of affiliations
            sups_={} # dictionary of the form { id : affiliation }
            for s in sups:
                s=s.split("</sup>") # list of the type [ id , affiliation ]
                sup_aff=str(s[1]).rstrip(" ").rstrip(",") # clean affiliation
                sups_[str(s[0])]=clean_address(sup_aff) # clean affiliation
            #insts=[ s.contents[0] for s in sups ]
            affs[str(rid)]=sups_
            
        elif len(labels) > 0:
            labels=str(c).split("</aff>")[0].split("<label>")[1:]
            sups_={}
            for s in labels:
                s=s.split("</label>")
                sup_aff=str(s[1]).rstrip(" ").rstrip(",")
                sups_[str(s[0])]=clean_address(sup_aff)
            #insts=[ s.contents[0] for s in sups ]
            affs[str(rid)]=sups_
            
        else:
            institution=c.findAll("institution")
            institution=[ s.contents[0] for s in institution ]
            institution="; ".join(institution)
            address=c.findAll("addr-line")
            address=[ s.contents[0] for s in address ]
            address="; ".join(address)
            affs[str(rid)]=str(", ".join([institution, address]).encode('ascii','ignore'))
    
    
    affiliations_dic={} # a dictionary of the form { author : affiliation }
    bad_authors=[] # authors with no affiliation
    
    # for all the authors we extracted from the PMC paper
    # we gone match the respect affiliations id/index
    # to a real affiliation address 
    for author in authors.keys():
        # if we have any index / id for this author we extract the affiliation
        # else, we see if we have one unique affiliation for all
        # otherwise we define it as none for this author
        if len(authors[author]) > 0:
            if type(authors[author]) == type({}):
                affiliations=[]
                keys=[ s for s in authors[author].keys() if s in affs.keys() ]
                for key in keys:
                    affs_vals=affs[key]
                    #print type(affs_vals), affs_vals
                    if type(affs_vals) != type({}):
                        affiliations.append(affs_vals)
                    #print affs_vals
                    #print authors[author][key]
                    else:
                        values=[ v for v in authors[author][key] if v in affs_vals.keys() ]
                        #print "vals", values
                        if len(values) > 0:
                            for value in values :
                                #print value
                                #print affs_vals[value]
                                affiliations.append(affs_vals[value])
                        else:
                            for v in affs_vals.keys():
                                affiliations.append(affs_vals[v])
                                
                affiliations="; ".join(affiliations)
            else:    
                ids=[ s.upper() for s in authors[author] if s in affs.keys() ]
                affiliations=[ affs[s] for s in ids ]
                affiliations="; ".join(affiliations)
        else:
            try:
                affiliations=affs["ALL"]
            except:
                #print "Failure on ", c
                #print url
                FAILURE=True
                affiliations=str(None)
                bad_authors.append(c)
        
        # we do a final clean to the final affiliation 
        # and add the { author : affiliation } to the affiliations_dic
        endaff=affiliations.rstrip(" ").rstrip(",")
        affiliations_dic[author]=endaff
        
        
        #if (len(endaff) == 0) | (endaff == ', ; '):
        #    print authors, affs
    #print affiliations_dic
    
    # check how many affiliations have len < 5 or are None
    audit = [ s for s in affiliations_dic.values() if len(s) <= 5 ]
    audit_b = [ s for s in affiliations_dic.values() if str(s) == "None" ]
    
    if (FAILURE) | (len(audit) > 0) :
        print("Issues with "+url)
        if (len(audit) > 0) :
            print("Some affiliations have less than 5 characteres")
        if FAILURE :
            try:
                print("Issues with the following authors:\n", "; ".join(bad_authors) )
            except:
                print("!!!\t", pmc_id, "Issues with the following authors:\n", bad_authors)
                
        #print "contrib\n", "\n".join(bs.article.front.findAll("contrib"))
        #print "aff\n", "\n".join(bs.article.front.findAll("aff"))
        breakit="\n\n##################################################################\n"
        print("\n\n", affiliations_dic, breakit)
        sys.stdout.flush()
        
    if len(affiliations_dic.keys()) ==0 :
        affiliations_dic=known_affiliations
        
    return affiliations_dic

def get_affiliations_text(article_id,author,article_values):
    """
    Gets the abstract from a paper in text format 
    and collects the affiliations of the authors.
    
    :param article_id: an article pubmed id
    :returns: a dictionary of the form { author : affiliation }
    
    """
    
    # URL
    efetch=baseurl+"efetch.fcgi?db=pubmed&rettype=json&retmode=json&id="
    url=efetch+article_id
    
    # request
    try:
        r = requests.post(url = url, verify=False)
    except:
        print(url)
        r = requests.post(url = url, verify=False)
    r = r.content.split("\n")
    
    r = [ s for s in r if len(s) != 0 ]

    # get first and last authors strings respectively 
    author_first=[ s for s in r if article_values['sortfirstauthor'] in s ][0]
    author_last=[ s for s in r if author.split(" ")[0] in s ][0]

    # get affiliations strings
    affiliations=[ s for s in r if s[0] == "(" ]
    affil_first=get_aff(author_first)
    affil_last=get_aff(author_last)

    affil_first=[ s for s in affiliations if s.split(")")[0]+")" in affil_first ]
    affil_last=[ s for s in affiliations if s.split(")")[0]+")" in affil_last ]

    return {article_values["affiliation_first"]:"; ".join(affil_first),\
           article_values["affiliation_last"]:"; ".join(affil_last) }

def FIXLAST_AUTHOR(df,funny_authors): #s=funny_authors
    """
    This function makes sure that if our author is last
    we do not have issues with non ascii characteres ...
    """
    author=df["author"]
    last=df["lastauthor"]
    index_author=[ a for a in funny_authors.keys() if author in a ]
    if len(index_author) == 0:
        if author in funny_authors.keys():
            if last in funny_authors[author]:
                last=author
    else:
        if index_author[0] in funny_authors.keys():
            if last in funny_authors[index_author[0]]:
                last=author
    return last


def build_publications_df(list_of_last_authors, funny_authors):
    baseurl="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    esearch=baseurl+"esearch.fcgi?db=pubmed&rettype=json&retmode=json&term="
    esummary=baseurl+"esummary.fcgi?db=pubmed&rettype=json&retmode=json&id="
    
    queries_df=pd.DataFrame() # dataframe with the list of ids for publications of each author
    publications_df=pd.DataFrame()

    for author in list_of_last_authors:
        print("\n", author)
        sys.stdout.flush()

        # fecth publications for author

        # URL
        author_=author.replace(" ","+")
        url=esearch+author_+"[author]&RetMax=1000"

        # request
        try:
            r = requests.post(url = url, verify=False)
        except:
            print(url)
            r = requests.post(url = url, verify=False)
        r = json.loads(r.content)

        # extract data from json
        data=r['esearchresult']
        querytranslation=data['querytranslation']
        idlist=data['idlist']
        count=int(data["count"])

        qdf=pd.DataFrame({"author":author,\
                          "querytranslation":querytranslation,\
                          "idlist":idlist,\
                          "count":count})
        queries_df=pd.concat([queries_df,qdf])



        authordf=pd.DataFrame() # dataframe with detailed description of publications
        for article_id in idlist:

            # fetch information on each paper

            # URL
            url=esummary+article_id

            # query often fails by server overload
            # so we do it slowly 
            works=False
            while not works:
                try:
                    # request
                    r = requests.post(url = url, verify=False)
                    r = json.loads(r.content)
                    works=True
                except:
                    print("sleeping for ", url)
                    time.sleep(2)
                    sys.stdout.flush()


            article_values={} # dictionary containing all the relevant data for each article

            pid=r['result'].keys()[0]
            data=r['result'][article_id]

            authors_=[ s["name"] for s in data["authors"] ]

            #print authors_

            # we make sure we are realy finding our author
            authorcheck=False

            funny_authors_={}
            for a in funny_authors.keys():
                funny_authors_[ a.upper() ]=[ s.upper() for s in funny_authors[a] ]
            authors_cap=[ s.upper() for s in authors_ ]



            # some author have non ascii names ..
            if author.upper() in funny_authors_.keys():
                for a in funny_authors_[author.upper()]:
                    if a in authors_cap:
                        authorcheck=True

            # we check if our author is in the list of authors in the paper metadata
            elif author.upper() in authors_cap:
                authorcheck=True

            # if we didn't find our author we move on and forget this paper
            if not authorcheck:
                print("Failed author check:\t", author, ",", article_id, ",", authors_)
                sys.stdout.flush()
                continue

            article_values["pubmed id"]=article_id

            # we check if this article is publicaly available
            try:
                pmc=data["articleids"]
                pmc=[ s["value"] for s in pmc if s["idtype"] == "pmc" ][0]
                article_values["pmc"]=pmc
            except:
                article_values["pmc"]=str(None)

            # this are the values we want to extract from the paper metadata
            values=['sortfirstauthor','lastauthor','source',\
                    'fulljournalname', 'pubdate','epubdate',\
                    'sortpubdate','pubtype','title','issn']                
            for value in values:
                v=data[value]
                if type(v) == list:
                    v=", ".join([ fix_non_ascii(s) for s in v ])
                article_values[value]=fix_non_ascii(v)

            # fetch paper metadata for affiliations
            affiliations=get_affiliations_json(article_id)

            # if we didn't find the exact autor in the affiliations we match it to the next possible name        
            if article_values['sortfirstauthor'] not in affiliations.keys():
                new_name=difflib.get_close_matches(article_values['sortfirstauthor'], affiliations.keys())
                if len(new_name) > 0:
                    article_values['sortfirstauthor']=new_name[0] 

            afirst=affiliations[article_values['sortfirstauthor']]
            if author in affiliations.keys():
                aa=affiliations[author]
            elif author.upper() in affiliations.keys():
                aa=affiliations[author.upper()]
            elif author.upper() in funny_authors_.keys() :
                i=funny_authors_.keys().index(author.upper())
                i=funny_authors.keys()[i]
                if funny_authors[i][0] in affiliations.keys():
                    aa=affiliations[funny_authors[i][0]]
                elif funny_authors[i][1] in affiliations.keys():
                    aa=affiliations[funny_authors[i][1]]


            # if we don't have affiliations for the first or last author and the
            # paper is publicaly available we read the paper for extracting the affiliations
            if ( ( afirst == str(None) ) | ( aa == str(None) ) ) & ( article_values["pmc"] != str(None) ):
                try:
                    affiliations=pmc_affiliations(article_values["pmc"],affiliations,tmpid=article_values["pmc"])
                except:
                    print(article_id,  "coult not run pmc_affiliations")
                    sys.stdout.flush()

                try:
                    afirst=affiliations[article_values['sortfirstauthor']]
                    if author in affiliations.keys():
                        aa=affiliations[author]
                    elif author.upper() in affiliations.keys():
                        aa=affiliations[author.upper()]
                    else:
                        i=funny_authors_.keys().index(author.upper())
                        i=funny_authors.keys()[i]
                        if funny_authors[i][0] in affiliations.keys():
                            aa=affiliations[funny_authors[i][0]]
                        else:
                            aa=affiliations[funny_authors[i][1]]
                except:
                    print("ISSUES:\n")
                    print(article_values["pmc"],"\n", article_values['sortfirstauthor'], "\n", affiliations.keys(),"\n")
                    sys.stdout.flush()

            article_values["affiliation_first"]=afirst          
            article_values["affiliation_author"]=aa

            # we create a string with all the authors and respective affiliations
            # author A aff: affiliations || author B aff: affiliations || author i aff: affiliations ..
            affiliations_str=[]
            for a in affiliations.keys():
                r=" aff: ".join([a,affiliations[a]])
                affiliations_str.append(r)
            affiliations_str=" || ".join(affiliations_str)
            article_values["affiliations"]=affiliations_str

            if author.upper() in funny_authors_.keys():
                i=funny_authors_.keys().index(author.upper())
                i=funny_authors.keys()[i]

                article_values["author"]=funny_authors[i][0]
            else:
                article_values["author"]=author

            artdf=pd.DataFrame(article_values, index=[0]) # dataframe for this article
            authordf=pd.concat([authordf,artdf]) # dataframe for this author

        publications_df=pd.concat([publications_df,authordf]) # dataframe for all authors
    return queries_df, publications_df



def check_institute(x, list_of_strings):
    """
    Given list of strings it looks for this strings in a given target string.
    """
    res="No"
    for l in list_of_strings:
        if l in str(x).lower():
            res="Yes"
    return res

def check_author_last(df, funny_authors):#=funny_authors
    """
    Checks if query author is last author
    """
    author=df["author"]
    last_author=df["lastauthor"]
    if author in funny_authors.keys():
        if last_author in funny_authors[author]:
            res="Yes"
        else:
            res="No"
    
    elif author in last_author:
        res="Yes"
    else:
        res="No"
    return res

def first_and_author(df):
    """
    Checks if first author and query author are from reference institution.
    """   
    first=df["1st author age"]
    author=df['author age']
    last=df["author is last"]
    if (first == "Yes") & (author == "Yes") & (last=="Yes"):
        res="Yes"
    else:
        res="No"
    return res

def first_author(df,funny_authors=None):
    """
    Checks if query author is first author.
    """
    first=df['sortfirstauthor']
    author=df['author']
    if author in funny_authors.keys():
        if first in funny_authors[author]:
            first=author
    if first==author:
        res="Yes"
    else:
        res="No"
    return res

def fix_pubdate(x):
    """
    Fixes publication date string to a timestamp
    """
    d=x.split(" ")
    year=d[0]
    if len(d) == 3:
        month=d[1]
        day=d[2]
        if "-" in day:
            day=day.split("-")[0]
        time_stamp=pd.Timestamp(year+"-"+month+"-"+day)
    elif len(d) == 2:
        month=d[1]
        if "-" in month:
            month=month.split("-")[0]
        try:
            time_stamp=pd.Timestamp(year+"-"+month)
        except:
            time_stamp=pd.Timestamp(year)
    elif len(d) == 1:
        time_stamp=pd.Timestamp(year+"-01")

    return time_stamp

def AUTHOR_IS_TRUE(df,safe_dic, kolle): # c=safe_dic, =kolle
    """
    Given a list of coauthors (safe_dic) and a list of substrings for the affiliations 
    check if the coauthors are present of the substring in the affiliations.
    """
    
    author=df["author"]
    authors=df["affiliations"].split(" || ")
    author_affiliation=str(df['affiliation_author'])
    authors=[ s.split(" aff: ")[0] for s in authors ]
    valid_authors=safe_dic[author]
    authors=[ s for s in authors if s in valid_authors ]
    res="No"
    for k in kolle:
        if k in author_affiliation.lower():
            res="Yes"
            break
    if len(authors) > 0:
        res="Yes"
    return res

    
def AUTHOR_IS_RE_TRUE(df,aff_check_dic,aff_wrong_dic): #=aff_check_dic, =aff_wrong_dic
    """
    Given a list of correct affiliations and a list of wrong affiliations
    we define if we are realy looking at the author of interest or not.
    """
    author=df["author"]
    author_affiliation=str(df['affiliation_author'])
    stat=df["true author"]
    if stat != "Yes":
        if author in aff_check_dic.keys():
            if author_affiliation in aff_check_dic[author]:
                stat="Yes"
    if stat != "Yes":
        if author in aff_wrong_dic.keys():
            if author_affiliation in aff_wrong_dic[author]:
                stat="Wrong aff."
    return stat

def STATS_A(publications_df, last_author=["Yes","No"], minY=2008, maxY=2019):
    """
    Given a publications dataframe it generates a table for the respective time frame
    :returns stats_a: a pandas dataframe with 
                "author", "npublications","ncorrect", "nwrongAff", "unknown"
    :returns stats_a_check: a pandas dataframe with the publications 
                where true author == No
    """
    
    
    stats_a=[]
    stats_a_check=pd.DataFrame()
    
    if (minY) | (maxY) :  
        minY=pd.Timestamp(str(minY)) 
        maxY=pd.Timestamp(str(maxY))
        tf=publications_df[ (publications_df["pubdate"]<maxY) & (publications_df["pubdate"]>=minY) ]
    
    for author in list(set(publications_df["author"].tolist())):

        tmp=tf[(tf["author"]==author) & \
               (tf["author is last"].isin(last_author) ) ]# & \
                            #(publications_df["pubdate"]<maxY) & \
                            #(publications_df["pubdate"]>=minY) ] 
            
        npublications=len(tmp)
        ncorrect=len(tmp[tmp["true author"]=="Yes"])
        nwrongAff=len(tmp[~tmp["true author"].isin(["Yes","No"])])
        # if affiliation is in defined as wrong the value for true author will be eg. Wrong aff.
        unknown=len(tmp[tmp["true author"]=="No"])

        stats_a_check=pd.concat([stats_a_check,tmp[tmp["true author"]=="No"]])

        stats_a.append([author, npublications, ncorrect, nwrongAff, unknown])
        #print author, npublications, ncorrect, nwrongAff, unknown
    stats_a=pd.DataFrame(stats_a,columns=["author", "npublications",\
                                          "ncorrect", "nwrongAff",\
                                          "unknown"])
    stats_a=stats_a.sort_values(by=["ncorrect"],ascending=False)
    stats_a

    for c in stats_a_check.columns.tolist():
        stats_a_check[c]= stats_a_check[c].apply( lambda x: clean_encoding(x) )
        
    return stats_a, stats_a_check

def useStatCheck(fin,publications_df):
    """
    Use the 'pubmed id','true author' from an excell file to 
    define if a publication belongs to an author or not.
    """
    
    stats_a_check_=pd.read_excel(fin)
    stats_a_check_=stats_a_check_[['pubmed id','true author']]
    stats_a_check_.index=stats_a_check_['pubmed id'].tolist()
    stats_a_check_=stats_a_check_.to_dict()["true author"]

    def Correct_stats_a(df, dic=stats_a_check_):
        pid=int(df['pubmed id'])
        v=df['true author']
        if pid in stats_a_check_.keys():
            v=stats_a_check_[pid]
        return v

    publications_df["true author"]=publications_df.apply(Correct_stats_a, axis=1)
    return publications_df

def create_links(x):
    r="https://www.ncbi.nlm.nih.gov/pubmed/?term="+str(x)
    return r
     
def FindBadAuthors(df):
    """
    Uses the information on wrong authorship attribution to 
    generate a list of authors which belong to a sosia. 
    
    :returns: a dictionary of the form { author : [ bad authors ] }
    """
    bddic={}
    for author in list(set(df["author"].tolist())):
        tmp=df[df["author"]==author]
        bad_authors=tmp[~tmp["true author"].isin(["Yes","No", "NoIF"])]["affiliations"].tolist()
        bad_authors=" || ".join(bad_authors)
        bad_authors=bad_authors.split(" || ")
        bad_authors=[ s.split(" aff: ")[0] for s in bad_authors ]
        bad_authors=[ s for s in bad_authors if s != author ]
        bddic[author]=bad_authors
        
    return bddic

def NegNetworkCorrection(df, dic):# =bddic
    """
    Uses information on secondary authors that belong to a sosia 
    for defining papers that belong to the author sosia.
    
    :param dic: a dictionary of the form { author : [ bad authors ] }
    """
    trueauthor=df["true author"]
    if trueauthor == "No":
        author=df["author"]
        affs=df["affiliations"]
        authors=affs.split(" || ")
        authors=[ s.split(" aff: ")[0] for s in authors ]
        authors=[ s for s in authors if s != author ]
        badauthors=[ s for s in authors if s in dic[author] ]
        if len(badauthors) > 0:
            trueauthor="Wrong col."
    return trueauthor
    
# we further manually curate unknown publications and expande negative social networks
def create_links(x):
    r="https://www.ncbi.nlm.nih.gov/pubmed/?term="+str(x)
    return r
    
def stats_b(publications_df, t, list_of_last_authors):
    """
    Creates stats on a publication dataframe
    
    :returns pub_counts: number of publications 
    :returns pub_if: sum of publications impact factor
    :returns pub_if_m: mean of publications impact factor
    """
    
    pub_counts=pd.DataFrame({ "author":list_of_last_authors } )
    pub_if=pd.DataFrame({ "author":list_of_last_authors } )
    pub_if_m=pd.DataFrame({ "author":list_of_last_authors } )

    #pub_counts["2008"]=0
    #pub_counts["2009"]=0

    for y in range(2008,2019):
        year_=pd.Timestamp(str(y+1)) 
        year=pd.Timestamp(str(y)) 
        pdf=publications_df[(publications_df["pubdate"]<year_) & \
                            (publications_df["pubdate"]>=year) &\
                            (publications_df["Journal Impact Factor"].astype(str)!="nan") & \
                            (publications_df["Journal Impact Factor"].astype(float)>0) & \
                           (publications_df["true author"]=="Yes")]
        
        if t == "all":
            pdf=pdf[pdf["author is last"].isin(["Yes","No"])]
        elif t == "first_or_last":
            pdf=pdf[ (pdf["author is last"].isin(["Yes"])) | \
                    (pdf["check_first_author"].isin(["Yes"]))]
        elif t == "first":
            pdf=pdf[ pdf["check_first_author"].isin(["Yes"]) ]
        elif t == "last":
            pdf=pdf[ (pdf["author is last"].isin(["Yes"]))]
                     
        year_df=[]
        sum_df=[]
        mean_df=[]
        for author in list(set(pdf["author"].tolist())):
            pdf_=pdf[pdf["author"]==author]
            sum_of_if=pdf_[["Journal Impact Factor"]].astype(float)
            sum_of_if=sum_of_if.dropna()
            sum_of_if=sum_of_if["Journal Impact Factor"].tolist()
            mean_of_if=np.mean(sum_of_if)
            sum_of_if=sum(sum_of_if)
            year_df.append([author, len(pdf_)])
            sum_df.append([author,sum_of_if])
            mean_df.append([author,mean_of_if])

        year_df=pd.DataFrame(year_df,columns=[ "author", str(y) ])
        pub_counts=pd.merge(pub_counts,year_df,on=["author"],how="outer")

        sum_df=pd.DataFrame(sum_df,columns=[ "author", str(y) ])
        pub_if=pd.merge(pub_if,sum_df,on=["author"],how="outer")

        mean_df=pd.DataFrame(mean_df,columns=[ "author", str(y) ])
        pub_if_m=pd.merge(pub_if_m,mean_df,on=["author"],how="outer")
        
    all_cols=pub_counts.columns.tolist()
    year_cols=[ s for s in all_cols if s != "author" ]
        
    pub_counts["sum"]=pub_counts[year_cols].sum(axis=1)
    pub_if["sum"]=pub_if[year_cols].sum(axis=1)
    pub_if_m["sum"]=pub_if_m[year_cols].sum(axis=1)

    pub_counts=pub_counts.sort_values(by=["sum"], ascending=False)
    pub_if=pub_if.sort_values(by=["sum"], ascending=False)
    pub_if_m=pub_if_m.sort_values(by=["sum"], ascending=False)

    all_cols.append("sum")

    pub_counts=pub_counts[all_cols]
    pub_if=pub_if[all_cols]
    pub_if_m=pub_if_m[all_cols]
    
    pub_counts.reset_index(inplace=True,drop=True)
    pub_if.reset_index(inplace=True,drop=True)
    pub_if_m.reset_index(inplace=True,drop=True)
        
    return pub_counts, pub_if, pub_if_m

def GetNetworkDFs(publications_df_file, files_tag, output_folder,\
                  min_number_of_publications=2, min_number_of_interactions=2,\
                 minY=2008, maxY=2019):

    publications_df=pd.read_table(publications_df_file)
    publications_df["pubdate"]=publications_df["pubdate"].apply(lambda x: fix_pubdate(x) )
    if (minY) | (maxY) :  
        print("Doing from", str(minY), "to", str(int(maxY)-1))
        sys.stdout.flush()
        minY=pd.Timestamp(str(minY)) 
        maxY=pd.Timestamp(str(maxY))
        publications_df=publications_df[ (publications_df["pubdate"]<maxY) & (publications_df["pubdate"]>=minY) ]

    df=publications_df[publications_df['true author']=="Yes"]
    people=df['affiliations'].tolist()
    people=" || ".join(people)
    people=people.split(" || ")
    people=[ s.split(" aff: ")[0] for s in people ]
    people=[[x,people.count(x)] for x in set(people)]
    people=pd.DataFrame(people)
    people.columns=["author","#"]
    people=people.sort_values(by=["#"],ascending=False)
    pdf=people[people["#"]>=min_number_of_publications]
    names=pdf["author"].tolist()

    df=publications_df[publications_df['true author']=="Yes"]
    people=df['affiliations'].tolist()
    def get_paper_names(x):
        x=x.split(" || ")
        x=[ s.split(" aff: ")[0] for s in x ]
        x=", ".join(x)
        return x
    people=[ get_paper_names(x) for x in people ]

    def get_coauthors(x, people=people, names=names):
        r=[ s for s in people if x in s ]
        r=", ".join(r)
        r=r.split(", ")
        r=[ s for s in r if s in names ]
        r=[": ".join([str(n),str(r.count(n))]) for n in set(r)]
        r=", ".join(r)
        return r

    pdf["coauthors"]=pdf["author"].apply(lambda x: get_coauthors(x) )


    cdf=pd.DataFrame()
    for a in pdf["author"].tolist():
        coauthors=pdf[pdf["author"]==a]["coauthors"].tolist()[0]
        coauthors=coauthors.split(", ")
        coauthors=[ s.split(": ") for s in coauthors ]
        coauthors=pd.DataFrame(coauthors, columns=["coauthor","n"])
        coauthors["author"]=a
        cdf=pd.concat([cdf,coauthors]) 

    vclist=[]
    vc=cdf[["coauthor","author"]]
    vc=list(vc.as_matrix())
    vc=[ list(s) for s in vc ]
    for v in vc:
        if ( v not in vclist ) & ( [v[1],v[0]] not in vclist ):
            vclist.append(v)

    # Partridge L, Isaacs AM

    check=[ s for s in vclist if "Partridge L" in s ]
    check=[ s for s in check if "Isaacs AM" in s ]  

    cdf_=cdf.copy()
    cdf=pd.DataFrame()
    for v in vclist:
        tmp=cdf_[ ( cdf_["coauthor"]==v[0] ) & ( cdf_["author"]==v[1] )]
        cdf=pd.concat([cdf,tmp])

    def KEEPERS(df):
        co=df["coauthor"]
        au=df["author"]
        if co==au:
            return False
        else:
            return True

    cdf["keep"]=cdf.apply(KEEPERS, axis=1)

    cdf=cdf[cdf["keep"]==True] 
    cdf=cdf.drop(["keep"],axis=1)
    cdf=cdf.sort_values(by=["n"],ascending=False)
    cdf.reset_index(inplace=True, drop=True)

    cdf=cdf[cdf["n"].astype(int)>=min_number_of_interactions]

    network_table=cdf.copy()
    network_table["n"]="pub"
    network_table.to_csv(output_folder+"network.authors."+files_tag+".sif",sep="\t", index=None, header=False)

    edges_table=cdf.copy()
    edges_table["shared name"]=edges_table["coauthor"]+" (pub) "+edges_table["author"]
    edges_table=edges_table[["shared name","n"]]
    edges_table.to_csv(output_folder+"edges.authors."+files_tag+".csv",index=None)

    nodes_table=pdf[["author","#"]]
    nodes_table.columns=["shared name","n"]
    nodes_table.to_csv(output_folder+"nodes.authors."+files_tag+".csv",index=None)
    
    def CheckAuthorNode(authors,author):
        if author in authors:
            return True
        else:
            return False

    def AuthorIFNode(author,publications_df):
        df=publications_df.copy()
        df["FOUND"]=df["affiliations"].apply(lambda x: CheckAuthorNode(x,author) )
        df=df[df["FOUND"]==True]
        IFs=df['Journal Impact Factor'].tolist()
        IFs=[ float(s) for s in IFs]
        IFs=[ s for s in IFs if str(s) not in ["nan"]]
        sumIFs=sum(IFs)
        if sumIFs == 0:
            medIFs=0.00
        else:
            medIFs=np.median(IFs)
        return str(sumIFs)+", "+str(medIFs)

    nodes_table["sum IF"]=nodes_table["shared name"].apply(lambda x: AuthorIFNode(x,publications_df))
    nodes_table["median IF"]=nodes_table["sum IF"].apply(lambda x: float(x.split(", ")[1] ))
    nodes_table["sum IF"]=nodes_table["sum IF"].apply(lambda x: float(x.split(", ")[0] ))
    
    def CheckAuthorEdge(authors,author):
        author=author.split(" (pub) ")
        if (author[0] in authors) & (author[1] in authors) :
            return True
        else:
            return False

    def AuthorIFEdge(author,publications_df):
        df=publications_df.copy()
        df["FOUND"]=df["affiliations"].apply(lambda x: CheckAuthorEdge(x,author) )
        df=df[df["FOUND"]==True]
        IFs=df['Journal Impact Factor'].tolist()
        IFs=[ float(s) for s in IFs]
        IFs=[ s for s in IFs if str(s) not in ["nan"]]
        sumIFs=sum(IFs)
        sumIFs=sum(IFs)
        if sumIFs == 0:
            medIFs=0
        else:
            medIFs=np.median(IFs)    
        return str(sumIFs)+", "+str(medIFs)

    edges_table["sum IF"]=edges_table["shared name"].apply(lambda x: AuthorIFEdge(x,publications_df))
    edges_table["median IF"]=edges_table["sum IF"].apply(lambda x: float(x.split(", ")[1] ))
    edges_table["sum IF"]=edges_table["sum IF"].apply(lambda x: float(x.split(", ")[0] ))
    
    return network_table, nodes_table, edges_table


# def PlotNetwork(network_file,remote_folder,output_folder,files_tag,\
#                network_table, nodes_table, edges_table,\
#                hostname="192.168.50.110",minY=2008, maxY=2019):

#     mypath=network_file
#     remotepath=remote_folder+network_file.split("/")[-1]
#     #remotepath="/Users/JBoucas/network.MPI-AGE.sif"

#     #ssh = paramiko.SSHClient()
#     #ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     #ssh.connect(hostname)
#     #ftp_client=ssh.open_sftp()
#     #ftp_client.put(mypath,remotepath)
    
#     # scp network file

#     call="scp "+mypath+" "+hostname+":"+remotepath
#     call=call.split(" ")
#     out=sb.call(call)

#     # connect to cytoscape and load network file
    
#     cytoscape=cyrest.cyclient(host=hostname)
#     nets=cytoscape.network.list()["networks"]
#     for n in nets:
#         cytoscape.network.destroy(network="SUID:"+str(n) )

#     cytoscape.vizmap.apply(styles="default")
#     cytoscape.network.load_file(remotepath)
    
    
#     # upload tables
    
#     nodes_table["n"]=nodes_table["n"].astype(float)
#     time.sleep(2)
#     res=cytoscape.table.loadTableData(nodes_table[["shared name","n", "sum IF", "median IF"]],table="node", df_key="shared name",\
#                                   table_key_column="shared name")
#     edges_table["n"]=edges_table["n"].astype(float)
#     res=cytoscape.table.loadTableData(edges_table[["shared name","n", "sum IF", "median IF"]],table="edge", df_key="shared name",\
#                                   table_key_column="shared name")
#     defaults_dic={"NODE_SHAPE":"ellipse",\
#                    "NODE_SIZE":"60",\
#                    "NODE_FILL_COLOR":"#AAAAAA",\
#                    "EDGE_TRANSPARENCY":"120"}
#     defaults_list=cytoscape.vizmap.simple_defaults(defaults_dic)

#     NODE_LABEL=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL",mappingType="passthrough",mappingColumn="name")

#     cytoscape.vizmap.create_style(title="dataStyle",defaults=defaults_list,mappings=[NODE_LABEL])
#     #sleep(2)
#     cytoscape.vizmap.apply(styles="dataStyle")

#     # node sized
    
#     nodes_table["n"]=nodes_table["n"].astype(float)
#     max_node=100#nodes_table[["n"]].quantile(0.99)[0]
#     min_node=5#nodes_table[["n"]].quantile(0.10)[0]
#     mean_node=(float(max_node)+float(min_node)) / 2.00

#     NODE_SIZE=cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_SIZE',\
#                                                  mappingType="continuous",\
#                                                  mappingColumn='n',\
#                                                  lower=[min_node,30.00],\
#                                                  center=[mean_node,(225.00+30.00)/2.00],\
#                                                  upper=[max_node,225.00],table="node")
#     cytoscape.vizmap.update_style(title="dataStyle",mappings=[NODE_SIZE])
#     cytoscape.vizmap.apply(styles="dataStyle")


#     # edge size
    
#     edges_table["n"]=edges_table["n"].astype(float)
#     max_edge=edges_table[["n"]].quantile(0.99)[0]
#     min_edge=edges_table[["n"]].quantile(0.10)[0]
#     mean_edge=(float(max_edge)+float(min_edge)) / 2.00

#     EDGE_SIZE=cytoscape.vizmap.mapVisualProperty(visualProperty='EDGE_WIDTH',\
#                                                  mappingType="continuous",\
#                                                  mappingColumn='n',\
#                                                  lower=[min_edge,0.75],\
#                                                  center=[mean_edge,(175.00+.75)/2.00],\
#                                                  upper=[max_edge,175.00],table="edge")
    
#     cytoscape.vizmap.update_style(title="dataStyle",mappings=[EDGE_SIZE])
#     cytoscape.vizmap.apply(styles="dataStyle")
#     time.sleep(2)
    
#     cytoscape.layout.force_directed(defaultSpringCoefficient=".000004", defaultSpringLength="5")
#     time.sleep(2)

#     # saving files
    
#     cytoscape.network.deselect(edgeList="all",nodeList="all")
#     time.sleep(2)

#     cytoscape.view.export(options="PNG", \
#                           OutputFile=remote_folder+"network.authors.bw."+files_tag+".png")
#     cytoscape.view.export(options="PDF", \
#                           OutputFile=remote_folder+"network.authors.bw."+files_tag+".pdf")
#     cytoscape.session.save_as(session_file=remote_folder+"network.authors.bw."+files_tag+".cys")
    
    
#     # adding colors to the nodes and edges
    
#     min_node_IF = min(nodes_table['median IF'].tolist())
#     med_node_IF = np.median(nodes_table['median IF'].tolist())
#     max_node_IF = nodes_table[['median IF']].quantile(0.75)[0]

#     cmap = matplotlib.cm.get_cmap("bwr")
#     norm = matplotlib.colors.Normalize(vmin=min_node_IF, vmax=max_node_IF)
#     min_color=matplotlib.colors.rgb2hex(cmap(norm(min_node_IF)))
#     center_color=matplotlib.colors.rgb2hex(cmap(norm(med_node_IF)))
#     max_color=matplotlib.colors.rgb2hex(cmap(norm(max_node_IF)))  

#     NODE_FILL_COLOR=cytoscape.vizmap.mapVisualProperty(visualProperty='NODE_FILL_COLOR',\
#                                                        mappingType="continuous",\
#                                                        mappingColumn='median IF',\
#                                                        lower=[min_node_IF,min_color],\
#                                                        center=[med_node_IF,center_color],\
#                                                        upper=[max_node_IF,max_color])
#     cytoscape.vizmap.update_style(title="dataStyle",mappings=[NODE_FILL_COLOR])
#     cytoscape.vizmap.apply(styles="dataStyle")
#     time.sleep(2)

#     min_node_IF = min(edges_table['median IF'].tolist())
#     med_node_IF = np.median(edges_table['median IF'].tolist())
#     max_node_IF = edges_table[['median IF']].quantile(0.75)[0]

#     cmap = matplotlib.cm.get_cmap("bwr")
#     norm = matplotlib.colors.Normalize(vmin=min_node_IF, vmax=max_node_IF)
#     min_color=matplotlib.colors.rgb2hex(cmap(norm(min_node_IF)))
#     center_color=matplotlib.colors.rgb2hex(cmap(norm(med_node_IF)))
#     max_color=matplotlib.colors.rgb2hex(cmap(norm(max_node_IF)))  

#     NODE_FILL_COLOR=cytoscape.vizmap.mapVisualProperty(visualProperty='EDGE_STROKE_UNSELECTED_PAINT',\
#                                                        mappingType="continuous",\
#                                                        mappingColumn='median IF',\
#                                                        lower=[min_node_IF,min_color],\
#                                                        center=[med_node_IF,center_color],\
#                                                        upper=[max_node_IF,max_color],\
#                                                       table="edge")
#     cytoscape.vizmap.update_style(title="dataStyle",mappings=[NODE_FILL_COLOR])
#     cytoscape.vizmap.apply(styles="dataStyle")
#     time.sleep(2)
    
#     cytoscape.network.deselect(edgeList="all",nodeList="all")
#     time.sleep(2)

#     cytoscape.view.export(options="PNG", \
#                           OutputFile=remote_folder+"network.authors.color."+files_tag+".png")
#     cytoscape.view.export(options="PDF", \
#                           OutputFile=remote_folder+"network.authors.color."+files_tag+".pdf")
#     cytoscape.session.save_as(session_file=remote_folder+"network.authors.color."+files_tag+".cys")
    
#     ssh = paramiko.SSHClient()
#     ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     ssh.connect(hostname)
#     ftp_client=ssh.open_sftp()
    
#     for f in ["network.authors.bw."+files_tag+".png",\
#              "network.authors.bw."+files_tag+".pdf",\
#              "network.authors.bw."+files_tag+".cys",\
#              "network.authors.color."+files_tag+".png",\
#              "network.authors.color."+files_tag+".pdf",\
#              "network.authors.color."+files_tag+".cys"]:
#         ftp_client.get(remote_folder+f,output_folder+f)
#         ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("rm "+remote_folder+f )
#     ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("rm "+remotepath )
    
#     return cytoscape