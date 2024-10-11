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


output_folder="/nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/parser_out"

baseurl="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
esearch=baseurl+"esearch.fcgi?db=pubmed&rettype=json&retmode=json&term="
esummary=baseurl+"esummary.fcgi?db=pubmed&rettype=json&retmode=json&id="

# {{{
list_of_last_authors=["Larsson NG","Partridge L","Antebi A","Langer T","Schaefer A",\
                      "Demetriades C","Graef M","Tessarz P","Valenzano DR","Wickstrom SA","Denzel MS",\
                     "Stewart JB","Matic I", "Pernas L", "Deelen J", "Frentz Z", "Huppertz I", "Jachimowicz RD",\
                     "Panier S"]
list_of_last_authors=["Langer T","Schaefer A",\
                      "Demetriades C","Graef M","Tessarz P","Valenzano DR","Wickstrom SA","Denzel MS",\
                     "Stewart JB","Matic I", "Pernas L", "Deelen J", "Frentz Z", "Huppertz I", "Jachimowicz RD",\
                     "Panier S"]

files_tag="MPI-AGE"
institute_sub_strings=[ "planck institute for biology of ag", "planck institute for the biology of ag"]
# we use the variable kolle to keep a list of substrings which we will use 
# for identifying safe authors. ie. if any affiliation of any of the authors
# contains one of these substrings we will define all authors on the paper
# as safe authors. based on these safe authors we then expande a social network
# which we define as safe.
kolle=[ "planck institute for biology of ag", "planck institute for the biology of ag", "koeln", "koln", "kln", "koeln", "cologne", "cecad", "ageing", "aging"]

# list_of_last_authors=["Hoehn M","Backes H", "Graf R", "Steculorum SM", "Korotkova T", "Bruning JC",\
#                      "Kornfeld JW", "Wunderlich TF", "Fenselau H", "Tittgemeyer M"]
# files_tag="MPI-MET"
# institute_sub_strings=[ "cecad", "institute for neurological research", "max planck institute for metabolism research", "cellular stress responses", "cologne" ]
# kolle=[ "planck institute for biology of ag", "institute for neurological research", "planck institute for the biology of ag", "koeln", "koln", "kln", "koeln", "cologne", "cecad", "ageing", "aging", "cecad", "max planck institute for metabolism research", "cellular stress responses" ]


funny_authors={ "Wickstrom SA" : [ u'Wickstr\xf6m SA' , "Wickstrom SA" , 'Wickstrm SA' ], "Pernas L": ["Pernas L","Pernas LF"] }
# }}}

# {{{
pmc="/nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/pmc"
pbmc_folders=os.listdir(pmc)
pmc_files=pd.DataFrame()
for f in  pbmc_folders:
    files=os.listdir(f"{pmc}/{f}")
    content=[]
    for f_ in files:
        with open(f"{pmc}/{f}/{f_}", "r" ) as txt :
            for l in txt:
                l=l.split("\t")[0]
                content.append(l)
    target=f_.replace("filelist.txt","tar.gz") 
    if f == "manuscript":
        target=f"https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/txt/{target}"
    else: 
        target=f"https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/{f}/txt/{target}"
    tmp=pd.DataFrame({"folder":f, "files":content, "tar.gz": target })
    pmc_files=pd.concat([pmc_files,tmp])
    
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

# {{{
# text='<?xml version="1.0" ?>\n<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2023//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_230101.dtd">\n<PubmedArticleSet>\n<PubmedArticle><MedlineCitation Status="MEDLINE" Owner="NLM" IndexingMethod="Automated"><PMID Version="1">37160881</PMID><DateCompleted><Year>2023</Year><Month>05</Month><Day>11</Day></DateCompleted><DateRevised><Year>2023</Year><Month>06</Month><Day>07</Day></DateRevised><Article PubModel="Electronic"><Journal><ISSN IssnType="Electronic">2041-1723</ISSN><JournalIssue CitedMedium="Internet"><Volume>14</Volume><Issue>1</Issue><PubDate><Year>2023</Year><Month>May</Month><Day>09</Day></PubDate></JournalIssue><Title>Nature communications</Title><ISOAbbreviation>Nat Commun</ISOAbbreviation></Journal><ArticleTitle>Artificial Hsp104-mediated systems for re-localizing protein aggregates.</ArticleTitle><Pagination><StartPage>2663</StartPage><MedlinePgn>2663</MedlinePgn></Pagination><ELocationID EIdType="pii" ValidYN="Y">2663</ELocationID><ELocationID EIdType="doi" ValidYN="Y">10.1038/s41467-023-37706-3</ELocationID><Abstract><AbstractText>Spatial Protein Quality Control (sPQC) sequesters misfolded proteins into specific, organelle-associated inclusions within the cell to control their toxicity. To approach the role of sPQC in cellular fitness, neurodegenerative diseases and aging, we report on the construction of Hsp100-based systems in budding yeast cells, which can artificially target protein aggregates to non-canonical locations. We demonstrate that aggregates of mutant huntingtin (mHtt), the disease-causing agent of Huntington\'s disease can be artificially targeted to daughter cells as well as to eisosomes and endosomes with this approach. We find that the artificial removal of mHtt inclusions from mother cells protects them from cell death suggesting that even large mHtt inclusions may be cytotoxic, a trait that has been widely debated. In contrast, removing inclusions of endogenous age-associated misfolded proteins does not significantly affect the lifespan of mother cells. We demonstrate also that this approach is able to manipulate mHtt inclusion formation in human cells and has the potential to be useful as an alternative, complementary approach to study the role of sPQC, for example in aging and neurodegenerative disease.</AbstractText><CopyrightInformation>&#xa9; 2023. The Author(s).</CopyrightInformation></Abstract><AuthorList CompleteYN="Y"><Author ValidYN="Y"><LastName>Fischbach</LastName><ForeName>Arthur</ForeName><Initials>A</Initials><Identifier Source="ORCID">0000-0001-6804-6564</Identifier><AffiliationInfo><Affiliation>Institute for Biomedicine, Sahlgrenska Academy, Centre for Ageing and Health-AgeCap, University of Gothenburg, Gothenburg, Sweden. arthur.fischbach@age.mpg.de.</Affiliation></AffiliationInfo><AffiliationInfo><Affiliation>Max-Planck Research Group Chromatin and Ageing, Max Planck Institute for Biology of Ageing, Cologne, Germany. arthur.fischbach@age.mpg.de.</Affiliation></AffiliationInfo></Author><Author ValidYN="Y"><LastName>Johns</LastName><ForeName>Angela</ForeName><Initials>A</Initials><Identifier Source="ORCID">0000-0003-3074-1208</Identifier><AffiliationInfo><Affiliation>Institute for Biomedicine, Sahlgrenska Academy, Centre for Ageing and Health-AgeCap, University of Gothenburg, Gothenburg, Sweden.</Affiliation></AffiliationInfo></Author><Author ValidYN="Y"><LastName>Schneider</LastName><ForeName>Kara L</ForeName><Initials>KL</Initials><Identifier Source="ORCID">0000-0001-6870-2939</Identifier><AffiliationInfo><Affiliation>Institute for Biomedicine, Sahlgrenska Academy, Centre for Ageing and Health-AgeCap, University of Gothenburg, Gothenburg, Sweden.</Affiliation></AffiliationInfo></Author><Author ValidYN="Y"><LastName>Hao</LastName><ForeName>Xinxin</ForeName><Initials>X</Initials><Identifier Source="ORCID">0000-0001-5758-6290</Identifier><AffiliationInfo><Affiliation>Institute for Biomedicine, Sahlgrenska Academy, Centre for Ageing and Health-AgeCap, University of Gothenburg, Gothenburg, Sweden.</Affiliation></AffiliationInfo></Author><Author ValidYN="Y"><LastName>Tessarz</LastName><ForeName>Peter</ForeName><Initials>P</Initials><Identifier Source="ORCID">0000-0002-6953-9835</Identifier><AffiliationInfo><Affiliation>Max-Planck Research Group Chromatin and Ageing, Max Planck Institute for Biology of Ageing, Cologne, Germany.</Affiliation></AffiliationInfo><AffiliationInfo><Affiliation>Cologne Excellence Cluster on Stress Responses in Ageing-Associated Diseases (CECAD), Cologne, Germany.</Affiliation></AffiliationInfo></Author><Author ValidYN="Y"><LastName>Nystr&#xf6;m</LastName><ForeName>Thomas</ForeName><Initials>T</Initials><Identifier Source="ORCID">0000-0001-5489-2903</Identifier><AffiliationInfo><Affiliation>Institute for Biomedicine, Sahlgrenska Academy, Centre for Ageing and Health-AgeCap, University of Gothenburg, Gothenburg, Sweden. thomas.nystrom@cmb.gu.se.</Affiliation></AffiliationInfo></Author></AuthorList><Language>eng</Language><PublicationTypeList><PublicationType UI="D016428">Journal Article</PublicationType><PublicationType UI="D013485">Research Support, Non-U.S. Gov\'t</PublicationType></PublicationTypeList><ArticleDate DateType="Electronic"><Year>2023</Year><Month>05</Month><Day>09</Day></ArticleDate></Article><MedlineJournalInfo><Country>England</Country><MedlineTA>Nat Commun</MedlineTA><NlmUniqueID>101528555</NlmUniqueID><ISSNLinking>2041-1723</ISSNLinking></MedlineJournalInfo><ChemicalList><Chemical><RegistryNumber>0</RegistryNumber><NameOfSubstance UI="D066329">Protein Aggregates</NameOfSubstance></Chemical></ChemicalList><CitationSubset>IM</CitationSubset><MeshHeadingList><MeshHeading><DescriptorName UI="D006801" MajorTopicYN="N">Humans</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D066329" MajorTopicYN="Y">Protein Aggregates</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D019636" MajorTopicYN="Y">Neurodegenerative Diseases</DescriptorName><QualifierName UI="Q000235" MajorTopicYN="N">genetics</QualifierName></MeshHeading><MeshHeading><DescriptorName UI="D000375" MajorTopicYN="N">Aging</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D008136" MajorTopicYN="N">Longevity</DescriptorName></MeshHeading><MeshHeading><DescriptorName UI="D016923" MajorTopicYN="N">Cell Death</DescriptorName></MeshHeading></MeshHeadingList><CoiStatement>The authors declare no competing interests.</CoiStatement></MedlineCitation><PubmedData><History><PubMedPubDate PubStatus="received"><Year>2022</Year><Month>4</Month><Day>12</Day></PubMedPubDate><PubMedPubDate PubStatus="accepted"><Year>2023</Year><Month>3</Month><Day>28</Day></PubMedPubDate><PubMedPubDate PubStatus="medline"><Year>2023</Year><Month>5</Month><Day>11</Day><Hour>6</Hour><Minute>42</Minute></PubMedPubDate><PubMedPubDate PubStatus="pubmed"><Year>2023</Year><Month>5</Month><Day>10</Day><Hour>6</Hour><Minute>42</Minute></PubMedPubDate><PubMedPubDate PubStatus="entrez"><Year>2023</Year><Month>5</Month><Day>10</Day><Hour>1</Hour><Minute>24</Minute></PubMedPubDate></History><PublicationStatus>epublish</PublicationStatus><ArticleIdList><ArticleId IdType="pubmed">37160881</ArticleId><ArticleId IdType="pmc">PMC10169802</ArticleId><ArticleId IdType="doi">10.1038/s41467-023-37706-3</ArticleId><ArticleId IdType="pii">10.1038/s41467-023-37706-3</ArticleId></ArticleIdList><ReferenceList><Reference><Citation>Gems D, Partridge L. Genetics of longevity in model organisms: debates and paradigm shifts. Annu. Rev. Physiol. 2013;75:621&#x2013;644. doi: 10.1146/annurev-physiol-030212-183712.</Citation><ArticleIdList><ArticleId IdType="doi">10.1146/annurev-physiol-030212-183712</ArticleId><ArticleId IdType="pubmed">23190075</ArticleId></ArticleIdList></Reference><Reference><Citation>Kirkwood TB. Understanding the odd science of aging. Cell. 2005;120:437&#x2013;447. doi: 10.1016/j.cell.2005.01.027.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.cell.2005.01.027</ArticleId><ArticleId IdType="pubmed">15734677</ArticleId></ArticleIdList></Reference><Reference><Citation>Vijg J, Campisi J. Puzzles, promises and a cure for ageing. Nature. 2008;454:1065&#x2013;1071. doi: 10.1038/nature07216.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature07216</ArticleId><ArticleId IdType="pmc">PMC2774752</ArticleId><ArticleId IdType="pubmed">18756247</ArticleId></ArticleIdList></Reference><Reference><Citation>Steinkraus KA, Kaeberlein M, Kennedy BK. Replicative aging in yeast: the means to the end. Annu. Rev. Cell Dev. Biol. 2008;24:29&#x2013;54. doi: 10.1146/annurev.cellbio.23.090506.123509.</Citation><ArticleIdList><ArticleId IdType="doi">10.1146/annurev.cellbio.23.090506.123509</ArticleId><ArticleId IdType="pmc">PMC2730916</ArticleId><ArticleId IdType="pubmed">18616424</ArticleId></ArticleIdList></Reference><Reference><Citation>Kaeberlein M. Lessons on longevity from budding yeast. Nature. 2010;464:513&#x2013;519. doi: 10.1038/nature08981.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature08981</ArticleId><ArticleId IdType="pmc">PMC3696189</ArticleId><ArticleId IdType="pubmed">20336133</ArticleId></ArticleIdList></Reference><Reference><Citation>Aguilaniu H, Gustafsson L, Rigoulet M, Nystrom T. Asymmetric inheritance of oxidatively damaged proteins during cytokinesis. Science. 2003;299:1751&#x2013;1753. doi: 10.1126/science.1080418.</Citation><ArticleIdList><ArticleId IdType="doi">10.1126/science.1080418</ArticleId><ArticleId IdType="pubmed">12610228</ArticleId></ArticleIdList></Reference><Reference><Citation>Erjavec N, Larsson L, Grantham J, Nystrom T. Accelerated aging and failure to segregate damaged proteins in Sir2 mutants can be suppressed by overproducing the protein aggregation-remodeling factor Hsp104p. Genes Dev. 2007;21:2410&#x2013;2421. doi: 10.1101/gad.439307.</Citation><ArticleIdList><ArticleId IdType="doi">10.1101/gad.439307</ArticleId><ArticleId IdType="pmc">PMC1993872</ArticleId><ArticleId IdType="pubmed">17908928</ArticleId></ArticleIdList></Reference><Reference><Citation>McFaline-Figueroa JR, et al. Mitochondrial quality control during inheritance is associated with lifespan and mother&#x2013;daughter age asymmetry in budding yeast. Aging Cell. 2011;10:885&#x2013;895. doi: 10.1111/j.1474-9726.2011.00731.x.</Citation><ArticleIdList><ArticleId IdType="doi">10.1111/j.1474-9726.2011.00731.x</ArticleId><ArticleId IdType="pmc">PMC3173513</ArticleId><ArticleId IdType="pubmed">21726403</ArticleId></ArticleIdList></Reference><Reference><Citation>Lai CY, Jaruga E, Borghouts C, Jazwinski SM. A mutation in the ATP2 gene abrogates the age asymmetry between mother and daughter cells of the yeast Saccharomyces cerevisiae. Genetics. 2002;162:73&#x2013;87. doi: 10.1093/genetics/162.1.73.</Citation><ArticleIdList><ArticleId IdType="doi">10.1093/genetics/162.1.73</ArticleId><ArticleId IdType="pmc">PMC1462265</ArticleId><ArticleId IdType="pubmed">12242224</ArticleId></ArticleIdList></Reference><Reference><Citation>Hughes AL, Gottschling DE. An early age increase in vacuolar pH limits mitochondrial function and lifespan in yeast. Nature. 2012;492:261&#x2013;265. doi: 10.1038/nature11654.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature11654</ArticleId><ArticleId IdType="pmc">PMC3521838</ArticleId><ArticleId IdType="pubmed">23172144</ArticleId></ArticleIdList></Reference><Reference><Citation>Sinclair DA, Guarente L. Extrachromosomal rDNA circles&#x2014;a cause of aging in yeast. Cell. 1997;91:1033&#x2013;1042. doi: 10.1016/S0092-8674(00)80493-6.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/S0092-8674(00)80493-6</ArticleId><ArticleId IdType="pubmed">9428525</ArticleId></ArticleIdList></Reference><Reference><Citation>Kruegel U, et al. Elevated proteasome capacity extends replicative lifespan in Saccharomyces cerevisiae. PLoS Genet. 2011;7:e1002253. doi: 10.1371/journal.pgen.1002253.</Citation><ArticleIdList><ArticleId IdType="doi">10.1371/journal.pgen.1002253</ArticleId><ArticleId IdType="pmc">PMC3169524</ArticleId><ArticleId IdType="pubmed">21931558</ArticleId></ArticleIdList></Reference><Reference><Citation>Janssens GE, et al. Protein biogenesis machinery is a driver of replicative aging in yeast. eLife. 2015;4:e08527. doi: 10.7554/eLife.08527.</Citation><ArticleIdList><ArticleId IdType="doi">10.7554/eLife.08527</ArticleId><ArticleId IdType="pmc">PMC4718733</ArticleId><ArticleId IdType="pubmed">26422514</ArticleId></ArticleIdList></Reference><Reference><Citation>Moreno, D. F. et al. Proteostasis collapse, a hallmark of aging, hinders the chaperone-Start network and arrests cells in G1. Elife8, 10.7554/eLife.48240 (2019).</Citation><ArticleIdList><ArticleId IdType="pmc">PMC6744273</ArticleId><ArticleId IdType="pubmed">31518229</ArticleId></ArticleIdList></Reference><Reference><Citation>Hill SM, Hao X, Liu B, Nystrom T. Life-span extension by a metacaspase in the yeast Saccharomyces cerevisiae. Science. 2014;344:1389&#x2013;1392. doi: 10.1126/science.1252634.</Citation><ArticleIdList><ArticleId IdType="doi">10.1126/science.1252634</ArticleId><ArticleId IdType="pubmed">24855027</ArticleId></ArticleIdList></Reference><Reference><Citation>Hanzen S, et al. Lifespan control by redox-dependent recruitment of chaperones to misfolded proteins. Cell. 2016;166:140&#x2013;151. doi: 10.1016/j.cell.2016.05.006.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.cell.2016.05.006</ArticleId><ArticleId IdType="pubmed">27264606</ArticleId></ArticleIdList></Reference><Reference><Citation>Hill S, et al. Asymmetric inheritance of aggregated proteins and age reset in yeast are regulated by Vac17-dependent vacuolar functions. Cell Rep. 2016;16:826&#x2013;838. doi: 10.1016/j.celrep.2016.06.016.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.celrep.2016.06.016</ArticleId><ArticleId IdType="pmc">PMC4963537</ArticleId><ArticleId IdType="pubmed">27373154</ArticleId></ArticleIdList></Reference><Reference><Citation>Dobson CM. Protein folding and misfolding. Nature. 2003;426:884&#x2013;890. doi: 10.1038/nature02261.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature02261</ArticleId><ArticleId IdType="pubmed">14685248</ArticleId></ArticleIdList></Reference><Reference><Citation>DeSantis ME, et al. Operational plasticity enables hsp104 to disaggregate diverse amyloid and nonamyloid clients. Cell. 2012;151:778&#x2013;793. doi: 10.1016/j.cell.2012.09.038.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.cell.2012.09.038</ArticleId><ArticleId IdType="pmc">PMC3496281</ArticleId><ArticleId IdType="pubmed">23141537</ArticleId></ArticleIdList></Reference><Reference><Citation>Tyedmers J, Mogk A, Bukau B. Cellular strategies for controlling protein aggregation. Nat. Rev. Mol. Cell Biol. 2010;11:777&#x2013;788. doi: 10.1038/nrm2993.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nrm2993</ArticleId><ArticleId IdType="pubmed">20944667</ArticleId></ArticleIdList></Reference><Reference><Citation>Hill SM, Hanzen S, Nystrom T. Restricted access: spatial sequestration of damaged proteins during stress and aging. EMBO Rep. 2017;18:377&#x2013;391. doi: 10.15252/embr.201643458.</Citation><ArticleIdList><ArticleId IdType="doi">10.15252/embr.201643458</ArticleId><ArticleId IdType="pmc">PMC5331209</ArticleId><ArticleId IdType="pubmed">28193623</ArticleId></ArticleIdList></Reference><Reference><Citation>Babazadeh R, et al. Syntaxin 5 is required for the formation and clearance of protein inclusions during proteostatic stress. Cell Rep. 2019;28:2096&#x2013;2110.e2098. doi: 10.1016/j.celrep.2019.07.053.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.celrep.2019.07.053</ArticleId><ArticleId IdType="pubmed">31433985</ArticleId></ArticleIdList></Reference><Reference><Citation>Escusa-Toret S, Vonk WI, Frydman J. Spatial sequestration of misfolded proteins by a dynamic chaperone pathway enhances cellular fitness during stress. Nat. Cell Biol. 2013;15:1231&#x2013;1243. doi: 10.1038/ncb2838.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/ncb2838</ArticleId><ArticleId IdType="pmc">PMC4121856</ArticleId><ArticleId IdType="pubmed">24036477</ArticleId></ArticleIdList></Reference><Reference><Citation>Kaganovich D, Kopito R, Frydman J. Misfolded proteins partition between two distinct quality control compartments. Nature. 2008;454:1088&#x2013;U1036. doi: 10.1038/nature07195.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature07195</ArticleId><ArticleId IdType="pmc">PMC2746971</ArticleId><ArticleId IdType="pubmed">18756251</ArticleId></ArticleIdList></Reference><Reference><Citation>Miller SBM, et al. Compartment-specific aggregases direct distinct nuclear and cytoplasmic aggregate deposition. EMBO J. 2015;34:778&#x2013;797. doi: 10.15252/embj.201489524.</Citation><ArticleIdList><ArticleId IdType="doi">10.15252/embj.201489524</ArticleId><ArticleId IdType="pmc">PMC4369314</ArticleId><ArticleId IdType="pubmed">25672362</ArticleId></ArticleIdList></Reference><Reference><Citation>Zhou C, et al. Organelle-based aggregation and retention of damaged proteins in asymmetrically dividing. Cells Cell. 2014;159:530&#x2013;542. doi: 10.1016/j.cell.2014.09.026.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.cell.2014.09.026</ArticleId><ArticleId IdType="pmc">PMC6726438</ArticleId><ArticleId IdType="pubmed">25417105</ArticleId></ArticleIdList></Reference><Reference><Citation>Specht S, Miller SBM, Mogk A, Bukau B. Hsp42 is required for sequestration of protein aggregates into deposition sites in Saccharomyces cerevisiae. J. Cell Biol. 2011;195:617&#x2013;629. doi: 10.1083/jcb.201106037.</Citation><ArticleIdList><ArticleId IdType="doi">10.1083/jcb.201106037</ArticleId><ArticleId IdType="pmc">PMC3257523</ArticleId><ArticleId IdType="pubmed">22065637</ArticleId></ArticleIdList></Reference><Reference><Citation>Shiber A, Breuer W, Brandeis M, Ravid T. Ubiquitin conjugation triggers misfolded protein sequestration into quality control foci when Hsp70 chaperone levels are limiting. Mol. Biol. Cell. 2013;24:2076&#x2013;2087. doi: 10.1091/mbc.e13-01-0010.</Citation><ArticleIdList><ArticleId IdType="doi">10.1091/mbc.e13-01-0010</ArticleId><ArticleId IdType="pmc">PMC3694792</ArticleId><ArticleId IdType="pubmed">23637465</ArticleId></ArticleIdList></Reference><Reference><Citation>Samant RS, Livingston CM, Sontag EM, Frydman J. Distinct proteostasis circuits cooperate in nuclear and cytoplasmic protein quality control. Nature. 2018;563:407&#x2013;411. doi: 10.1038/s41586-018-0678-x.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/s41586-018-0678-x</ArticleId><ArticleId IdType="pmc">PMC6707801</ArticleId><ArticleId IdType="pubmed">30429547</ArticleId></ArticleIdList></Reference><Reference><Citation>Spokoini R, et al. Confinement to organelle-associated inclusion structures mediates asymmetric inheritance of aggregated protein in budding yeast. Cell Rep. 2012;2:738&#x2013;747. doi: 10.1016/j.celrep.2012.08.024.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.celrep.2012.08.024</ArticleId><ArticleId IdType="pubmed">23022486</ArticleId></ArticleIdList></Reference><Reference><Citation>Johnston JA, Ward CL, Kopito RR. Aggresomes: a cellular response to misfolded proteins. J. Cell Biol. 1998;143:1883&#x2013;1898. doi: 10.1083/jcb.143.7.1883.</Citation><ArticleIdList><ArticleId IdType="doi">10.1083/jcb.143.7.1883</ArticleId><ArticleId IdType="pmc">PMC2175217</ArticleId><ArticleId IdType="pubmed">9864362</ArticleId></ArticleIdList></Reference><Reference><Citation>Christian Wigley W, et al. Dynamic association of proteasomal machinery with the centrosome. J. Cell Biol. 1999;145:481&#x2013;490. doi: 10.1083/jcb.145.3.481.</Citation><ArticleIdList><ArticleId IdType="doi">10.1083/jcb.145.3.481</ArticleId><ArticleId IdType="pmc">PMC2185077</ArticleId><ArticleId IdType="pubmed">10225950</ArticleId></ArticleIdList></Reference><Reference><Citation>Liu B, et al. The polarisome is required for segregation and retrograde transport of protein aggregates. Cell. 2010;140:257&#x2013;267. doi: 10.1016/j.cell.2009.12.031.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.cell.2009.12.031</ArticleId><ArticleId IdType="pubmed">20141839</ArticleId></ArticleIdList></Reference><Reference><Citation>Saarikangas J, Barral Y. Protein aggregates are associated with replicative aging without compromising protein quality control. Elife. 2015;4:e06197. doi: 10.7554/eLife.06197.</Citation><ArticleIdList><ArticleId IdType="doi">10.7554/eLife.06197</ArticleId><ArticleId IdType="pmc">PMC4635334</ArticleId><ArticleId IdType="pubmed">26544680</ArticleId></ArticleIdList></Reference><Reference><Citation>King GA, et al. Meiotic cellular rejuvenation is coupled to nuclear remodeling in budding yeast. eLife. 2019;8:e47156. doi: 10.7554/eLife.47156.</Citation><ArticleIdList><ArticleId IdType="doi">10.7554/eLife.47156</ArticleId><ArticleId IdType="pmc">PMC6711709</ArticleId><ArticleId IdType="pubmed">31397671</ArticleId></ArticleIdList></Reference><Reference><Citation>Unal E, Kinde B, Amon A. Gametogenesis eliminates age-induced cellular damage and resets life span in yeast. Science. 2011;332:1554&#x2013;1557. doi: 10.1126/science.1204349.</Citation><ArticleIdList><ArticleId IdType="doi">10.1126/science.1204349</ArticleId><ArticleId IdType="pmc">PMC3923466</ArticleId><ArticleId IdType="pubmed">21700873</ArticleId></ArticleIdList></Reference><Reference><Citation>D&#xfc;nkler A, et al. Type V myosin focuses the polarisome and shapes the tip of yeast cells. J. Cell Biol. 2021;220:e202006193. doi: 10.1083/jcb.202006193.</Citation><ArticleIdList><ArticleId IdType="doi">10.1083/jcb.202006193</ArticleId><ArticleId IdType="pmc">PMC7933982</ArticleId><ArticleId IdType="pubmed">33656555</ArticleId></ArticleIdList></Reference><Reference><Citation>Krobitsch S, Lindquist S. Aggregation of huntingtin in yeast varies with the length of the polyglutamine expansion and the expression of chaperone proteins. Proc. Natl Acad. Sci. USA. 2000;97:1589&#x2013;1594. doi: 10.1073/pnas.97.4.1589.</Citation><ArticleIdList><ArticleId IdType="doi">10.1073/pnas.97.4.1589</ArticleId><ArticleId IdType="pmc">PMC26479</ArticleId><ArticleId IdType="pubmed">10677504</ArticleId></ArticleIdList></Reference><Reference><Citation>Winter Georg E, et al. Phthalimide conjugation as a strategy for in vivo target protein degradation. Science. 2015;348:1376&#x2013;1381. doi: 10.1126/science.aab1433.</Citation><ArticleIdList><ArticleId IdType="doi">10.1126/science.aab1433</ArticleId><ArticleId IdType="pmc">PMC4937790</ArticleId><ArticleId IdType="pubmed">25999370</ArticleId></ArticleIdList></Reference><Reference><Citation>Sakamoto KM, et al. Protacs: Chimeric molecules that target proteins to the Skp1&#x2013;Cullin&#x2013;F box complex for ubiquitination and degradation. Proc. Natl Acad. Sci. USA. 2001;98:8554. doi: 10.1073/pnas.141230798.</Citation><ArticleIdList><ArticleId IdType="doi">10.1073/pnas.141230798</ArticleId><ArticleId IdType="pmc">PMC37474</ArticleId><ArticleId IdType="pubmed">11438690</ArticleId></ArticleIdList></Reference><Reference><Citation>Banik SM, et al. Lysosome-targeting chimaeras for degradation of extracellular proteins. Nature. 2020;584:291&#x2013;297. doi: 10.1038/s41586-020-2545-9.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/s41586-020-2545-9</ArticleId><ArticleId IdType="pmc">PMC7727926</ArticleId><ArticleId IdType="pubmed">32728216</ArticleId></ArticleIdList></Reference><Reference><Citation>Comyn SA, Young BP, Loewen CJ, Mayor T. Prefoldin promotes proteasomal degradation of cytosolic proteins with missense mutations by maintaining substrate solubility. PLOS Genet. 2016;12:e1006184. doi: 10.1371/journal.pgen.1006184.</Citation><ArticleIdList><ArticleId IdType="doi">10.1371/journal.pgen.1006184</ArticleId><ArticleId IdType="pmc">PMC4957761</ArticleId><ArticleId IdType="pubmed">27448207</ArticleId></ArticleIdList></Reference><Reference><Citation>Masser AE, Kandasamy G, Kaimal JM, Andr&#xe9;asson C. Luciferase NanoLuc as a reporter for gene expression and protein levels in Saccharomyces cerevisiae. Yeast. 2016;33:191&#x2013;200. doi: 10.1002/yea.3155.</Citation><ArticleIdList><ArticleId IdType="doi">10.1002/yea.3155</ArticleId><ArticleId IdType="pmc">PMC5069653</ArticleId><ArticleId IdType="pubmed">26860732</ArticleId></ArticleIdList></Reference><Reference><Citation>Mackay RG, Helsen CW, Tkach JM, Glover JR. The C-terminal extension of Saccharomyces cerevisiae Hsp104 plays a role in oligomer assembly. Biochemistry. 2008;47:1918&#x2013;1927. doi: 10.1021/bi701714s.</Citation><ArticleIdList><ArticleId IdType="doi">10.1021/bi701714s</ArticleId><ArticleId IdType="pubmed">18197703</ArticleId></ArticleIdList></Reference><Reference><Citation>Schirmer EC, Queitsch C, Kowal AS, Parsell DA, Lindquist S. The ATPase activity of Hsp104, effects of environmental conditions and mutations. J. Biol. Chem. 1998;273:15546&#x2013;15552. doi: 10.1074/jbc.273.25.15546.</Citation><ArticleIdList><ArticleId IdType="doi">10.1074/jbc.273.25.15546</ArticleId><ArticleId IdType="pubmed">9624144</ArticleId></ArticleIdList></Reference><Reference><Citation>Saarikangas J, Barral Y. Protein aggregation as a mechanism of adaptive cellular responses. Curr. Genet. 2016;62:711&#x2013;724. doi: 10.1007/s00294-016-0596-0.</Citation><ArticleIdList><ArticleId IdType="doi">10.1007/s00294-016-0596-0</ArticleId><ArticleId IdType="pubmed">27032776</ArticleId></ArticleIdList></Reference><Reference><Citation>Ho CT, et al. Cellular sequestrases maintain basal Hsp70 capacity ensuring balanced proteostasis. Nat. Commun. 2019;10:4851. doi: 10.1038/s41467-019-12868-1.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/s41467-019-12868-1</ArticleId><ArticleId IdType="pmc">PMC6813348</ArticleId><ArticleId IdType="pubmed">31649258</ArticleId></ArticleIdList></Reference><Reference><Citation>Song J, et al. Essential genetic interactors of SIR2 required for spatial sequestration and asymmetrical inheritance of protein aggregates. PLOS Genet. 2014;10:e1004539. doi: 10.1371/journal.pgen.1004539.</Citation><ArticleIdList><ArticleId IdType="doi">10.1371/journal.pgen.1004539</ArticleId><ArticleId IdType="pmc">PMC4117435</ArticleId><ArticleId IdType="pubmed">25079602</ArticleId></ArticleIdList></Reference><Reference><Citation>Wang J, et al. Progressive aggregation despite chaperone associations of a mutant SOD1-YFP in transgenic mice that develop ALS. Proc. Natl Acad. Sci. USA. 2009;106:1392&#x2013;1397. doi: 10.1073/pnas.0813045106.</Citation><ArticleIdList><ArticleId IdType="doi">10.1073/pnas.0813045106</ArticleId><ArticleId IdType="pmc">PMC2631083</ArticleId><ArticleId IdType="pubmed">19171884</ArticleId></ArticleIdList></Reference><Reference><Citation>M&#xfc;nch C, Bertolotti A. Exposure of hydrophobic surfaces initiates aggregation of diverse ALS-causing superoxide dismutase-1 mutants. J. Mol. Biol. 2010;399:512&#x2013;525. doi: 10.1016/j.jmb.2010.04.019.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.jmb.2010.04.019</ArticleId><ArticleId IdType="pmc">PMC2927901</ArticleId><ArticleId IdType="pubmed">20399791</ArticleId></ArticleIdList></Reference><Reference><Citation>Rothbauer U, et al. A versatile nanotrap for biochemical and functional studies with fluorescent fusion proteins. Mol. Cell Proteom. 2008;7:282&#x2013;289. doi: 10.1074/mcp.M700342-MCP200.</Citation><ArticleIdList><ArticleId IdType="doi">10.1074/mcp.M700342-MCP200</ArticleId><ArticleId IdType="pubmed">17951627</ArticleId></ArticleIdList></Reference><Reference><Citation>Huh WK, et al. Global analysis of protein localization in budding yeast. Nature. 2003;425:686&#x2013;691. doi: 10.1038/nature02026.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature02026</ArticleId><ArticleId IdType="pubmed">14562095</ArticleId></ArticleIdList></Reference><Reference><Citation>Ruan L, et al. Cytosolic proteostasis through importing of misfolded proteins into mitochondria. Nature. 2017;543:443&#x2013;446. doi: 10.1038/nature21695.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature21695</ArticleId><ArticleId IdType="pmc">PMC5793917</ArticleId><ArticleId IdType="pubmed">28241148</ArticleId></ArticleIdList></Reference><Reference><Citation>Danielli L, Li X, Tuller T, Daniel R. Quantifying the distribution of protein oligomerization degree reflects cellular information capacity. Sci. Rep. 2020;10:17689. doi: 10.1038/s41598-020-74811-5.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/s41598-020-74811-5</ArticleId><ArticleId IdType="pmc">PMC7573690</ArticleId><ArticleId IdType="pubmed">33077848</ArticleId></ArticleIdList></Reference><Reference><Citation>Carcamo WC, et al. Induction of cytoplasmic rods and rings structures by inhibition of the CTP and GTP synthetic pathway in mammalian cells. PLoS ONE. 2011;6:e29690. doi: 10.1371/journal.pone.0029690.</Citation><ArticleIdList><ArticleId IdType="doi">10.1371/journal.pone.0029690</ArticleId><ArticleId IdType="pmc">PMC3248424</ArticleId><ArticleId IdType="pubmed">22220215</ArticleId></ArticleIdList></Reference><Reference><Citation>Chen K, et al. Glutamine analogs promote cytoophidium assembly in human and Drosophila cells. J. Genet. Genomics. 2011;38:391&#x2013;402. doi: 10.1016/j.jgg.2011.08.004.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.jgg.2011.08.004</ArticleId><ArticleId IdType="pubmed">21930098</ArticleId></ArticleIdList></Reference><Reference><Citation>Narayanaswamy R, et al. Widespread reorganization of metabolic enzymes into reversible assemblies upon nutrient starvation. Proc. Natl Acad. Sci. USA. 2009;106:10147&#x2013;10152. doi: 10.1073/pnas.0812771106.</Citation><ArticleIdList><ArticleId IdType="doi">10.1073/pnas.0812771106</ArticleId><ArticleId IdType="pmc">PMC2691686</ArticleId><ArticleId IdType="pubmed">19502427</ArticleId></ArticleIdList></Reference><Reference><Citation>Suresh HG, et al. Prolonged starvation drives reversible sequestration of lipid biosynthetic enzymes and organelle reorganization in Saccharomyces cerevisiae. Mol. Biol. Cell. 2015;26:1601&#x2013;1615. doi: 10.1091/mbc.E14-11-1559.</Citation><ArticleIdList><ArticleId IdType="doi">10.1091/mbc.E14-11-1559</ArticleId><ArticleId IdType="pmc">PMC4436773</ArticleId><ArticleId IdType="pubmed">25761633</ArticleId></ArticleIdList></Reference><Reference><Citation>Noree C, Sato BK, Broyer RM, Wilhelm JE. Identification of novel filament-forming proteins in Saccharomyces cerevisiae and Drosophila melanogaster. J. cell Biol. 2010;190:541&#x2013;551. doi: 10.1083/jcb.201003001.</Citation><ArticleIdList><ArticleId IdType="doi">10.1083/jcb.201003001</ArticleId><ArticleId IdType="pmc">PMC2928026</ArticleId><ArticleId IdType="pubmed">20713603</ArticleId></ArticleIdList></Reference><Reference><Citation>Shen Q-J, et al. Filamentation of metabolic enzymes in Saccharomyces cerevisiae. J. Genet. Genomics. 2016;43:393&#x2013;404. doi: 10.1016/j.jgg.2016.03.008.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.jgg.2016.03.008</ArticleId><ArticleId IdType="pmc">PMC4920916</ArticleId><ArticleId IdType="pubmed">27312010</ArticleId></ArticleIdList></Reference><Reference><Citation>Sathyanarayanan U, et al. ATP hydrolysis by yeast Hsp104 determines protein aggregate dissolution and size in vivo. Nat. Commun. 2020;11:5226. doi: 10.1038/s41467-020-19104-1.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/s41467-020-19104-1</ArticleId><ArticleId IdType="pmc">PMC7568574</ArticleId><ArticleId IdType="pubmed">33067463</ArticleId></ArticleIdList></Reference><Reference><Citation>Mogk A, Ruger-Herreros C, Bukau B. Cellular functions and mechanisms of action of small heat shock proteins. Annu Rev. Microbiol. 2019;73:89&#x2013;110. doi: 10.1146/annurev-micro-020518-115515.</Citation><ArticleIdList><ArticleId IdType="doi">10.1146/annurev-micro-020518-115515</ArticleId><ArticleId IdType="pubmed">31091419</ArticleId></ArticleIdList></Reference><Reference><Citation>Grousl, T. et al. A prion-like domain in Hsp42 drives chaperone-facilitated aggregation of misfolded proteins. J. Cell Biol.10.1083/jcb.201708116 (2018).</Citation><ArticleIdList><ArticleId IdType="pmc">PMC5881502</ArticleId><ArticleId IdType="pubmed">29362223</ArticleId></ArticleIdList></Reference><Reference><Citation>Gourlay CW, Carpp LN, Timpson P, Winder SJ, Ayscough KR. A role for the actin cytoskeleton in cell death and aging in yeast. J. Cell Biol. 2004;164:803&#x2013;809. doi: 10.1083/jcb.200310148.</Citation><ArticleIdList><ArticleId IdType="doi">10.1083/jcb.200310148</ArticleId><ArticleId IdType="pmc">PMC2172293</ArticleId><ArticleId IdType="pubmed">15024029</ArticleId></ArticleIdList></Reference><Reference><Citation>Gheysen D, et al. Assembly and release of HIV-1 precursor Pr55gag virus-like particles from recombinant baculovirus-infected insect cells. Cell. 1989;59:103&#x2013;112. doi: 10.1016/0092-8674(89)90873-8.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/0092-8674(89)90873-8</ArticleId><ArticleId IdType="pubmed">2676191</ArticleId></ArticleIdList></Reference><Reference><Citation>Votteler J, et al. Designed proteins induce the formation of nanocage-containing extracellular vesicles. Nature. 2016;540:292&#x2013;295. doi: 10.1038/nature20607.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature20607</ArticleId><ArticleId IdType="pmc">PMC5729044</ArticleId><ArticleId IdType="pubmed">27919066</ArticleId></ArticleIdList></Reference><Reference><Citation>Melentijevic I, et al. C. elegans neurons jettison protein aggregates and mitochondria under neurotoxic stress. Nature. 2017;542:367&#x2013;371. doi: 10.1038/nature21362.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature21362</ArticleId><ArticleId IdType="pmc">PMC5336134</ArticleId><ArticleId IdType="pubmed">28178240</ArticleId></ArticleIdList></Reference><Reference><Citation>Kim YE, et al. Soluble oligomers of PolyQ-expanded huntingtin target a multiplicity of key cellular factors. Mol. Cell. 2016;63:951&#x2013;964. doi: 10.1016/j.molcel.2016.07.022.</Citation><ArticleIdList><ArticleId IdType="doi">10.1016/j.molcel.2016.07.022</ArticleId><ArticleId IdType="pubmed">27570076</ArticleId></ArticleIdList></Reference><Reference><Citation>Arrasate M, Mitra S, Schweitzer ES, Segal MR, Finkbeiner S. Inclusion body formation reduces levels of mutant huntingtin and the risk of neuronal death. Nature. 2004;431:805&#x2013;810. doi: 10.1038/nature02998.</Citation><ArticleIdList><ArticleId IdType="doi">10.1038/nature02998</ArticleId><ArticleId IdType="pubmed">15483602</ArticleId></ArticleIdList></Reference><Reference><Citation>Winzeler EA, et al. Functional characterization of the S. cerevisiae genome by gene deletion and parallel analysis. Science. 1999;285:901&#x2013;906. doi: 10.1126/science.285.5429.901.</Citation><ArticleIdList><ArticleId IdType="doi">10.1126/science.285.5429.901</ArticleId><ArticleId IdType="pubmed">10436161</ArticleId></ArticleIdList></Reference><Reference><Citation>Tong AH, et al. Systematic genetic analysis with ordered arrays of yeast deletion mutants. Science. 2001;294:2364&#x2013;2368. doi: 10.1126/science.1065810.</Citation><ArticleIdList><ArticleId IdType="doi">10.1126/science.1065810</ArticleId><ArticleId IdType="pubmed">11743205</ArticleId></ArticleIdList></Reference><Reference><Citation>Egilmez NK, Chen JB, Jazwinski SM. Preparation and partial characterization of old yeast cells. J. Gerontol. 1990;45:B9&#x2013;B17. doi: 10.1093/geronj/45.1.B9.</Citation><ArticleIdList><ArticleId IdType="doi">10.1093/geronj/45.1.B9</ArticleId><ArticleId IdType="pubmed">2404060</ArticleId></ArticleIdList></Reference></ReferenceList></PubmedData></PubmedArticle></PubmedArticleSet>'
# affiliatons={}
# authors=text.split("<AuthorList")[1].split("</AuthorList>")[0].split(">", 1)[1]
# authors=authors.split("</Author>")
# authors=[ s.split("<Author")[1].split(">", 1)[1] for s in authors if s ]


# def field_parser(x, i, f=None):
#     if not f:
#         f="/"+i
#     return x.split(f"<{i}>")[1].split(f"<{f}>")[0]
    

# def get_author(x):
#     lastname=field_parser(x, "LastName")
#     initials=field_parser(x, "Initials")
#     affiliations=x.split("</Affiliation></AffiliationInfo>")
#     affiliations=[ s.split("<AffiliationInfo><Affiliation>")[1] for s in affiliations if s ]
#     affiliations="; ".join(affiliations)
#     author=f"{lastname} {initials}"
#     return author, affiliations
    
        
# for a in authors:
#     author, affiliations=get_author(a)
#     affiliatons[author]=affiliations
# affiliatons
# }}}

# {{{

html.unescape('K&#xfc;hl I')
# }}}

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
    
    text=html.unescape(str(r.content))
    
    if "Error: External viewer error: Empty Response." in text:
        print("\n\n!!\n\n", article_id, text, "\n\n!!\n\n")
        return {}
    
    affiliations={}
    try:
        if '<AuthorList Type="authors"' in text:
            authors=text.split('<AuthorList Type="authors"')[1].split("</AuthorList>")[0].split(">", 1)[1]
        else:
            authors=text.split("<AuthorList")[1].split("</AuthorList>")[0].split(">", 1)[1]
            
        authors=authors.split("</Author>")
        authors=[ s.split("<Author")[1].split(">", 1)[1] for s in authors if s ]
    except:
        print(text)
        authors=text.split("<AuthorList")[1].split("</AuthorList>")[0].split(">", 1)[1]
        authors=authors.split("</Author>")
        authors=[ s.split("<Author")[1].split(">", 1)[1] for s in authors if s ]
        

    def field_parser(x, i, f=None):
        if not f:
            f="/"+i
        return x.split(f"<{i}>")[1].split(f"<{f}>")[0]

    def get_author(x):
        if ( "LastName" in x ) and ("Initials" in x ):
            try:
                lastname=field_parser(x, "LastName")
                initials=field_parser(x, "Initials")
                if "Affiliation" in x:
                    affiliations=x.split("</Affiliation></AffiliationInfo>")
                    affiliations=[ s.split("<AffiliationInfo><Affiliation>")[1] for s in affiliations if s ]
                    affiliations="; ".join(affiliations)
                else:
                    affiliations=str(None)
                author=f"{lastname} {initials}"
            except:
                print(x)
                lastname=field_parser(x, "LastName")
                initials=field_parser(x, "Initials")
                affiliations=x.split("</Affiliation></AffiliationInfo>")
                affiliations=[ s.split("<AffiliationInfo><Affiliation>")[1] for s in affiliations if s ]
                affiliations="; ".join(affiliations)
                author=f"{lastname} {initials}"
            return author, affiliations
        else:
            print("! no LastName or no Initials", x)
            return None, None
            
        

    for a in authors:
        author, aff=get_author(a)
        affiliations[author]=aff
            
    # if str(article_id) =="21413227" :
    #     print(text,"\n", affiliations)
            
    return affiliations

# {{{
queries_df=pd.DataFrame() # dataframe with the list of ids for publications of each author
publications_df=pd.DataFrame()

for author in list_of_last_authors:
    print(author)
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
    for article_id in idlist: #["21413227"]:
    
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
        
        pid=list(r['result'].keys())[0]
        data=r['result'][article_id]
        
        authors_=[ s["name"] for s in data["authors"] ]
        
        # we make sure we are realy finding our author
        authorcheck=False
        
        # some author have non ascii names ..
        if author in funny_authors.keys():
            for a in funny_authors[author]:
                if a in authors_:
                    authorcheck=True
                    
        # we check if our author is in the list of authors in the paper metadata
        elif author in authors_:
            authorcheck=True
        
        # if we didn't find our author we move on and forget this paper
        if not authorcheck:
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
                v=", ".join([ ps.fix_non_ascii(s) for s in v ])
            article_values[value]=ps.fix_non_ascii(v)
        
        # fetch paper metadata for affiliations
        # affiliations=ps.get_affiliations_json(article_id)
        affiliations=get_affiliations_json(article_id)

        if not affiliations:
            print(f"No affiliations for {url}")
            continue
        
        # if we didn't find the exact author in the affiliations we match it to the next possible name        
        if article_values['sortfirstauthor'] not in list(affiliations.keys()):
            new_name=difflib.get_close_matches(article_values['sortfirstauthor'], affiliations.keys())
            if len(new_name) > 0:
                article_values['sortfirstauthor']=new_name[0]
            else:
                print("!!", url, "\n\t", article_values['sortfirstauthor'], "\n\t" ,list(affiliations.keys()) )
        
        # print(article_values, affiliations)
        afirst=affiliations[article_values['sortfirstauthor']]
        
        if author in affiliations.keys():
            aa=affiliations[author]
        elif author in funny_authors.keys():
            for a in funny_authors[author]:
                if a in affiliations.keys():
                    aa=affiliations[a]
         
        
        # if we don't have affiliations for the first or last author and the
        # paper is publicaly available we read the paper for extracting the affiliations
        if ( ( afirst == str(None) ) | ( aa == str(None) ) ) & ( article_values["pmc"] != str(None) ):
            try:
                affiliations=ps.pmc_affiliations(article_values["pmc"],affiliations,tmpid=article_values["pmc"])
            except:
                print(article_id,  "coult not run pmc_affiliations")
            
            try:
                afirst=affiliations[article_values['sortfirstauthor']]
                aa=affiliations[author]
            except:
                print("ISSUES:\n")
                print(article_values["pmc"],"\n", article_values['sortfirstauthor'], "\n", affiliations.keys(),"\n")
            
        article_values["affiliation_first"]=afirst          
        article_values["affiliation_author"]=aa
        
        # we create a string with all the authors and respective affiliations
        # author A aff: affiliations || author B aff: affiliations || author i aff: affiliations ..
        affiliations_str=[]
        for a in affiliations.keys():
            r=" aff: ".join([ str(a),str(affiliations[a]) ])
            affiliations_str.append(r)
        affiliations_str=" || ".join(affiliations_str)
        article_values["affiliations"]=affiliations_str
        
        article_values["author"]=author
        
        artdf=pd.DataFrame(article_values, index=[0]) # dataframe for this article
        authordf=pd.concat([authordf,artdf]) # dataframe for this author
    
    publications_df=pd.concat([publications_df,authordf]) # dataframe for all authors
publications_df=publications_df.reset_index(inplace=False, drop=True)
# }}}
publications_df

# {{{
# for a in publications_df["affiliation_author"]:
#     try:
#         ps.check_institute(str(a))
#     except:
#         print("!error!", a, tryp)
# }}}

# {{{
publications_df.reset_index(inplace=True, drop=True)
# for c in ['affiliation_author', 'affiliation_first', 'affiliations', 'author', 'epubdate', 'fulljournalname', 'issn', 'lastauthor', 'pmc', 'pubmed id', 'pubtype', 'sortfirstauthor', 'sortpubdate', 'source', 'title']:
#         publications_df[c]= publications_df[c].apply( lambda x: ps.clean_encoding(x) )

# fix publication date
publications_df["pubdate"]=publications_df["pubdate"].apply(lambda x: ps.fix_pubdate(x) )
        
publications_df["lastauthor"]=publications_df.apply(ps.FIXLAST_AUTHOR, args=(funny_authors,), axis=1)

# check if author affiliation includes the institute of interest
publications_df["author age"]=publications_df["affiliation_author"].apply(lambda x: ps.check_institute(x, institute_sub_strings) )

# check if first author affiliation includes institute of interest 
publications_df["1st author age"]=publications_df["affiliation_first"].apply(lambda x: ps.check_institute(x, institute_sub_strings) )

# check if target author is the last author
publications_df["author is last"]=publications_df.apply(ps.check_author_last, args=(funny_authors,), axis=1)

# check if author is first author
publications_df["first_and_author"]=publications_df.apply(ps.first_and_author, axis=1)
# }}}

publications_df.head()


publications_df=pd.merge(publications_df,pmc_files,on=["pmc"], how="left")

# {{{
bash='''
#!/bin/bash

cd /nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/
mkdir -p pmc_papers
cd pmc_papers

#################


'''

papers=publications_df.dropna(subset=["tar.gz"])

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

        bash=bash+ash

print(bash)

with open("download.extract.pmc.sh", "w" ) as sh:
    sh.write(bash)
        
        
    
# }}}
list( set(publications_df["tar.gz"].tolist() ))[:2]



