# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 09:26:57 2022

@author: luol2

Pipeline: first gene NER, then species assignment
input: species NER bioc xml file
output: gene ner and species assignment results bioc xml file
"""
import argparse
import os
import io
import time
import sys
import re
import shutil
from src_python.GeneNER import model_ner,ner_tag
from src_python.SpeAss import model_sa,sa_tag

import tensorflow as tf

import bioc
import stanza
nlp_token = stanza.Pipeline(model_dir='gnorm_trained_models/stanza', lang='en', processors={'tokenize': 'spacy'},package='None', download_method=None) #package='craft' ;./gnorm_trained_models/stanza

def NER_BioC(infolder,infile,outpath,nn_model):
         
    with open(infolder+"/"+infile, 'r',encoding='utf-8') as fin:
        with open(outpath+"/"+infile,'w', encoding='utf8') as fout:
            collection = bioc.load(fin)
            
            Total_n=len(collection.documents)
            print('Total number of sub-documents:', Total_n)
            pmid_n=0
            for document in collection.documents:
                print("Processing:{0}%".format(round(pmid_n * 100 / Total_n)), end="\r")
                pmid_n+=1
                # print(document.id)
                mention_num_new=0
                for passage in document.passages:
                    if passage.text!='' and (not passage.text.isspace()) and passage.infons['type']!='ref': # have text and is not ref
                        passage_offset=passage.offset
                        tag_result=ner_tag.ML_Tag(passage.text,nn_model,nlp_token)
                        mention_num=0
                        for ele in tag_result:
                            bioc_note = bioc.BioCAnnotation()
                            bioc_note.id = str(mention_num)
                            mention_num+=1
                            bioc_note.infons['type'] = ele[2]
                            start = int(ele[0])
                            last = int(ele[1])
                            loc = bioc.BioCLocation(offset=str(passage_offset+start), length= str(last-start))
                            bioc_note.locations.append(loc)
                            bioc_note.text = passage.text[start:last]
                            passage.annotations.append(bioc_note)
                    #update id
                    for temp_annotation in passage.annotations:
                        temp_annotation.id=str(mention_num_new)
                        mention_num_new+=1
            bioc.dump(collection, fout, pretty_print=True)
            
def NER_PubTator(infolder,infile,outpath,nn_model):
    with open(infolder+"/"+infile, 'r',encoding='utf-8') as fin:
        with open(outpath+"/"+infile,'w', encoding='utf-8') as fout:
            title=''
            abstract=''
            all_text=fin.read().strip().split('\n\n')
            Total_n=len(all_text)
            print('Total number of sub-documents:', Total_n)
            pmid_n=0
            for doc in all_text:
                print("Processing:{0}%".format(round(pmid_n * 100 / Total_n)), end="\r")
                pmid_n+=1
                lines = doc.split('\n')
                seg=lines[0].split('|t|')
                pmid=seg[0]
                title=seg[1]
                seg=lines[1].split('|a|')
                abstract=seg[1]
                
                intext=title+' '+abstract
                tag_result=ner_tag.ML_Tag(intext,nn_model,nlp_token)
                fout.write(doc+'\n')
                for ele in tag_result:
                    ent_start = ele[0]
                    ent_last = ele[1]
                    ent_mention = intext[int(ele[0]):int(ele[1])]
                    ent_type=ele[2]
                    fout.write(pmid+"\t"+ent_start+"\t"+ent_last+"\t"+ent_mention+"\t"+ent_type+"\n")
                fout.write('\n')
                title=''
                abstract=''
                
def geneNER(infolder, outpath, modelfile):
    
    print('loading NER models........')    
        
    if modelfile.lower().find('bioformer')>=0:
        vocabfiles={'labelfile':'./vocab/GeneNER_label.vocab',
                        'checkpoint_path':'./gnorm_trained_models/bioformer-cased-v1.0/', #bioformer-cased-v1.0
                        'lowercase':False,
                        } 
    else:
        vocabfiles={'labelfile':'./vocab/GeneNER_label.vocab',
                    'checkpoint_path':'./gnorm_trained_models/BiomedNLP-PubMedBERT-base-uncased-abstract/',
                    'lowercase':True,
                    } 
        
    nn_model=model_ner.HUGFACE_NER(vocabfiles)
    nn_model.build_encoder()
    nn_model.build_softmax_decoder()
    nn_model.load_model(modelfile)
        
    #tagging text
    print("begin GeneNER tagging........")
    start_time=time.time()
    
    for infile in os.listdir(infolder):
        if os.path.isfile(outpath+"/"+infile):
            print(infile+' has exsited.')
        else:
            print('processing:',infile)             
            fin = open(infolder+"/"+infile, 'r',encoding='utf-8')
            input_format=""
            for line in fin:
                pattern_bioc = re.compile('.*<collection>.*')
                pattern_pubtator = re.compile('^([^\|]+)\|[^\|]+\|(.*)')
                if pattern_pubtator.search(line):
                    input_format="PubTator"
                    break
                elif pattern_bioc.search(line):
                    input_format="BioC"
                    break
            fin.close()
            if(input_format == "PubTator"):
                NER_PubTator(infolder,infile,outpath,nn_model)
            elif(input_format == "BioC"):
                NER_BioC(infolder,infile,outpath,nn_model)    
    
    print('tag done:',time.time()-start_time)
    

#SA for bioc format
def SA_BioC(infolder,infile,outpath,nn_model,virus_set,prefix_dict):
    
    #BioC xml to pubtator
    # pmid|t|text1
    #pmid|a|text2
    #pmid sid eid entity_txt entity_type entity_id (gene is blank)
    fin = open(infolder+"/"+infile, 'r',encoding='utf-8')
    # fout_pubtator=open(outpath+'tmp/input_xml.pubtator','w', encoding='utf-8')
    fin_pubtator0=io.StringIO() #none *species
    fin_pubtator1=io.StringIO() #one *species
    fin_pubtator2=io.StringIO() #two or more species
    collection = bioc.load(fin)
    fin.close()
    ori_ann_index={}  #{'pmid':{'ent.id':'ent_s-ent_e'}}
    species_count={} #{pmid:{speid:num}}
    gene_set=['Gene','FamilyName']
    final_sa_results={} #{'pmid':{'entity_id':species_id}}
    for document in collection.documents:
        doc_pmid=document.id
        doc_title=''
        doc_abstract=''
        doc_annotation=[]
        _ann_index={}
        _species_num={} #{*speciesid:num}
        _gene_num=0
        _passage_num=0              
        if len(document.passages)<=2: #abstract xml or PMC only have title
            for passage in document.passages:
                passage_offset=passage.offset
                _passage_num+=1              
                #print(passage_offset,type(passage_offset))
                #if passage.infons['type']=='title' or passage.infons['type']=='front':
                if _passage_num==1:
                    doc_title=passage.text
                    for temp_annotation in passage.annotations:
                        if temp_annotation.infons['type'] in gene_set:
                            _gene_num+=1
                        ent_start=temp_annotation.locations[0].offset-passage_offset
                        ent_end=ent_start+temp_annotation.locations[0].length
                        #print(ent_start,ent_end)
                        _ann_index[temp_annotation.id]=str(ent_start)+'-'+str(ent_end)
                        # print(temp_annotation.infons)
                        if 'Identifier' in temp_annotation.infons.keys():
                            # print(temp_annotation.infons.keys['Identifier'])
                            species_ID=temp_annotation.infons['Identifier']
                            if species_ID.find('*')>=0:
                                if species_ID not in _species_num.keys():
                                    _species_num[species_ID]=1
                                else:
                                    _species_num[species_ID]+=1
                            doc_annotation.append(doc_pmid+'\t'+temp_annotation.id+'\t'+str(ent_start)+'\t'+str(ent_end)+'\t'+temp_annotation.text+'\t'+temp_annotation.infons['type']+'\t'+species_ID)
                        else:
                            doc_annotation.append(doc_pmid+'\t'+temp_annotation.id+'\t'+str(ent_start)+'\t'+str(ent_end)+'\t'+temp_annotation.text+'\t'+temp_annotation.infons['type'])
                        
                #elif passage.infons['type']=='abstract' or passage.infons['type']=='paragraph':
                else:
                    doc_abstract=passage.text
                    for temp_annotation in passage.annotations:
                        if temp_annotation.infons['type'] in gene_set:
                            _gene_num+=1
                        ent_start=len(doc_title)+1+temp_annotation.locations[0].offset-passage_offset
                        ent_end=ent_start+temp_annotation.locations[0].length
                        #print(ent_start,ent_end)
                        _ann_index[temp_annotation.id]=str(ent_start)+'-'+str(ent_end)
                        if 'Identifier' in temp_annotation.infons.keys():
                            # print(temp_annotation.infons.keys['Identifier'])
                            species_ID=temp_annotation.infons['Identifier']
                            if species_ID.find('*')>=0:
                                if species_ID not in _species_num.keys():
                                    _species_num[species_ID]=1
                                else:
                                    _species_num[species_ID]+=1
                            doc_annotation.append(doc_pmid+'\t'+temp_annotation.id+'\t'+str(ent_start)+'\t'+str(ent_end)+'\t'+temp_annotation.text+'\t'+temp_annotation.infons['type']+'\t'+species_ID)
                        else:
                            doc_annotation.append(doc_pmid+'\t'+temp_annotation.id+'\t'+str(ent_start)+'\t'+str(ent_end)+'\t'+temp_annotation.text+'\t'+temp_annotation.infons['type'])
            
            if len(_species_num)>=2 and _gene_num>0:
                fin_pubtator2.write(doc_pmid+'|t|'+doc_title+'\n')
                fin_pubtator2.write(doc_pmid+'|a|'+doc_abstract+'\n')
                for ele in doc_annotation:
                    fin_pubtator2.write(ele+'\n')
                fin_pubtator2.write('\n')
            elif len(_species_num)==1 and _gene_num>0: #可以直接给结果
                fin_pubtator1.write(doc_pmid+'|t|'+doc_title+'\n')
                fin_pubtator1.write(doc_pmid+'|a|'+doc_abstract+'\n')
                major_speicesid,=_species_num
                fin_pubtator1.write(major_speicesid[1:]+'\n')
                for ele in doc_annotation:
                    fin_pubtator1.write(ele+'\n')
                fin_pubtator1.write('\n')
            elif len(_species_num)==0 and _gene_num>0:
                fin_pubtator0.write(doc_pmid+'|t|'+doc_title+'\n')
                fin_pubtator0.write(doc_pmid+'|a|'+doc_abstract+'\n')
                for ele in doc_annotation:
                    fin_pubtator0.write(ele+'\n')
                fin_pubtator0.write('\n')        
        
        else: # full text xml
            for passage in document.passages:
                passage_annotation=[]
                _species_num_passage={}
                _gene_num_passage=0
                passage_offset=passage.offset
                #print(passage_offset,type(passage_offset))
                if passage.text!='' and (not passage.text.isspace()) and passage.infons['type']!='ref':
                    doc_title=passage.text
                    for temp_annotation in passage.annotations:
                        if temp_annotation.infons['type'] in gene_set:
                            _gene_num_passage+=1
                        ent_start=temp_annotation.locations[0].offset-passage_offset
                        ent_end=ent_start+temp_annotation.locations[0].length
                        #print(ent_start,ent_end)
                        _ann_index[temp_annotation.id]=str(ent_start)+'-'+str(ent_end)
                        # print(temp_annotation.infons)
                        if 'Identifier' in temp_annotation.infons.keys():
                            # print(temp_annotation.infons.keys['Identifier'])
                            species_ID=temp_annotation.infons['Identifier']
                            if species_ID.find('*')>=0:
                                if species_ID not in _species_num.keys():
                                    _species_num[species_ID]=1
                                else:
                                    _species_num[species_ID]+=1
                                if species_ID not in _species_num_passage.keys():
                                    _species_num_passage[species_ID]=1
                                else:
                                    _species_num_passage[species_ID]+=1
                            passage_annotation.append(doc_pmid+'\t'+temp_annotation.id+'\t'+str(ent_start)+'\t'+str(ent_end)+'\t'+temp_annotation.text+'\t'+temp_annotation.infons['type']+'\t'+species_ID)
                        else:
                            passage_annotation.append(doc_pmid+'\t'+temp_annotation.id+'\t'+str(ent_start)+'\t'+str(ent_end)+'\t'+temp_annotation.text+'\t'+temp_annotation.infons['type'])
                        
            
                if len(_species_num_passage)>=2 and _gene_num_passage>0:
                    fin_pubtator2.write(doc_pmid+'|t|'+doc_title+'\n')
                    fin_pubtator2.write(doc_pmid+'|a|'+doc_abstract+'\n')
                    for ele in passage_annotation:
                        fin_pubtator2.write(ele+'\n')
                    fin_pubtator2.write('\n')
                elif len(_species_num_passage)==1 and _gene_num_passage>0: #可以直接给结果
                    fin_pubtator1.write(doc_pmid+'|t|'+doc_title+'\n')
                    fin_pubtator1.write(doc_pmid+'|a|'+doc_abstract+'\n')
                    major_speicesid,=_species_num_passage
                    fin_pubtator1.write(major_speicesid[1:]+'\n')
                    for ele in passage_annotation:
                        fin_pubtator1.write(ele+'\n')
                    fin_pubtator1.write('\n')
                elif len(_species_num_passage)==0 and _gene_num_passage>0:
                    fin_pubtator0.write(doc_pmid+'|t|'+doc_title+'\n')
                    fin_pubtator0.write(doc_pmid+'|a|'+doc_abstract+'\n')
                    for ele in passage_annotation:
                        fin_pubtator0.write(ele+'\n')
                    fin_pubtator0.write('\n')        
        # print(ori_ann_index)
        
        ori_ann_index[doc_pmid]=_ann_index
        species_count[doc_pmid]=_species_num
    
            
    cache_geneid={} #{pmid:{gene1:{id1:num,id2:num}}}
    
    if fin_pubtator2.getvalue()!='':
        #pubtator format  ML tagging
        # print(fin_pubtator2.getvalue())
        ml_out= sa_tag.ml_tag_main(fin_pubtator2,nlp_token, nn_model)   
        #print(ml_out.getvalue())
        fin_result=io.StringIO(ml_out.getvalue())
        all_in=fin_result.read().strip().split('\n\n')
        #print('+2 species:',len(all_in))
        fin_result.close()
        
        prefix_speid_allset=set(prefix_dict.keys())
    
        for doc in all_in:
            lines=doc.split('\n')
            pmid=lines[0].split('|t|')[0]
            _prefix_str2id_dict={}
            doc_species=list(species_count[pmid].keys())
            for _spe_ele in doc_species:
                if _spe_ele[1:] in prefix_speid_allset:
                    for ele in prefix_dict[_spe_ele[1:]]:
                        _prefix_str2id_dict[ele]=_spe_ele[1:]
            
            for i in range(2,len(lines)):
                segs=lines[i].split('\t')
                if pmid not in final_sa_results.keys():
                    final_sa_results[pmid]={segs[1]:'Focus:'+segs[-1]}
                else:
                    final_sa_results[pmid][segs[1]]='Focus:'+segs[-1]
                    
                if segs[5] in gene_set:
                    if segs[4][0:2] in _prefix_str2id_dict: #prefix rule
                        #print('prefix rule:', pmid)
                        # print(_prefix_str2id_dict)
                        if pmid not in final_sa_results.keys():
                            final_sa_results[pmid]={segs[1]:'Focus:'+_prefix_str2id_dict[segs[4][0:2]]}
                        else:
                            final_sa_results[pmid][segs[1]]='Focus:'+_prefix_str2id_dict[segs[4][0:2]]
                    if pmid not in cache_geneid.keys():
                        cache_geneid[pmid]={segs[4]:{'Focus:'+segs[-1]:1}}
                    else:
                        if segs[4] not in cache_geneid[pmid].keys():
                            cache_geneid[pmid][segs[4]]={'Focus:'+segs[-1]:1}
                        else:
                            if segs[-1] not in cache_geneid[pmid][segs[4]].keys():
                                cache_geneid[pmid][segs[4]]['Focus:'+segs[-1]]=1
                            else:
                                cache_geneid[pmid][segs[4]]['Focus:'+segs[-1]]+=1
        
    #print(final_sa_results)
    
    #one species  
    if fin_pubtator1.getvalue()!='':
        fin_result=io.StringIO(fin_pubtator1.getvalue())
        all_in=fin_result.read().strip().split('\n\n')
        fin_result.close()
        #print('1 species:',len(all_in))
        for doc in all_in:
            lines=doc.split('\n')
            pmid=lines[0].split('|t|')[0]
            major_speicesid=lines[2]
            for i in range(3,len(lines)):
                segs=lines[i].split('\t')
                if len(segs)>=7:#species
                    if pmid not in final_sa_results.keys():
                        final_sa_results[pmid]={segs[1]:segs[-1]}
                    else:
                        final_sa_results[pmid][segs[1]]=segs[-1]
                else:#gene
                    marjor_species='Focus:'+major_speicesid
                    if pmid not in final_sa_results.keys():
                        final_sa_results[pmid]={segs[1]:marjor_species}
                    else:
                        final_sa_results[pmid][segs[1]]=marjor_species
                    if pmid not in cache_geneid.keys():
                        cache_geneid[pmid]={segs[4]:{marjor_species:1}}
                    else:
                        if segs[4] not in cache_geneid[pmid].keys():
                            cache_geneid[pmid][segs[4]]={marjor_species:1}
                        else:
                            if segs[-1] not in cache_geneid[pmid][segs[4]].keys():
                                cache_geneid[pmid][segs[4]][marjor_species]=1
                            else:
                                cache_geneid[pmid][segs[4]][marjor_species]+=1

    
    #no species
    fin_result=io.StringIO(fin_pubtator0.getvalue())
    all_in=fin_result.read().strip().split('\n\n')
    fin_result.close()
    #print('no species:',len(all_in))
    for doc in all_in:
        lines=doc.split('\n')
        pmid=lines[0].split('|t|')[0]

        for i in range(2,len(lines)):
            segs=lines[i].split('\t')
            if (pmid in cache_geneid.keys()) and (segs[4] in cache_geneid[pmid].keys()):#same gene in doc
                marjor_species = max(zip(cache_geneid[pmid][segs[4]].values(), cache_geneid[pmid][segs[4]].keys()))
                if pmid not in final_sa_results.keys():
                    final_sa_results[pmid]={segs[1]:marjor_species[1]}
                else:
                    final_sa_results[pmid][segs[1]]=marjor_species[1]
            else: #marjor species in doc
                if (pmid in species_count.keys()) and len(species_count[pmid])>0:#marjor species in doc
                    marjor_species = max(zip(species_count[pmid].values(), species_count[pmid].keys()))

                    if pmid not in final_sa_results.keys():
                        final_sa_results[pmid]={segs[1]:'Focus:'+marjor_species[1][1:]}
                    else:
                        final_sa_results[pmid][segs[1]]='Focus:'+marjor_species[1][1:]
                else:#no any species in doc,assign human
                    if pmid not in final_sa_results.keys():
                        final_sa_results[pmid]={segs[1]:'Focus:9606'}
                    else:
                        final_sa_results[pmid][segs[1]]='Focus:9606'

    
    
    # print(final_sa_results)
    fin = open(infolder+"/"+infile, 'r',encoding='utf-8')
    fout_xml=open(outpath+"/"+infile,'w', encoding='utf8')
    collection = bioc.load(fin)
    for document in collection.documents:
        doc_pmid=document.id
        # print(final_sa_results[doc_pmid])
        # print(doc_pmid)
        for passage in document.passages:
            for temp_annotation in passage.annotations:
                if 'Identifier' not in temp_annotation.infons.keys():
                    if temp_annotation.id in final_sa_results[doc_pmid].keys():
                        if final_sa_results[doc_pmid][temp_annotation.id][6:] in virus_set:
                            temp_annotation.infons['Identifier']=final_sa_results[doc_pmid][temp_annotation.id]+',9606'
                            # print('!!! virus:', doc_pmid)
                        else:
                            temp_annotation.infons['Identifier']=final_sa_results[doc_pmid][temp_annotation.id]
                    else: #same text bug
                        if (doc_pmid in cache_geneid.keys()) and (temp_annotation.text in cache_geneid[doc_pmid].keys()):#same gene in doc
                            marjor_species = max(zip(cache_geneid[doc_pmid][temp_annotation.text].values(), cache_geneid[doc_pmid][temp_annotation.text].keys()))
                            temp_annotation.infons['Identifier']=marjor_species[1]
                        else: 
                            
                            temp_annotation.infons['Identifier']='Focus:9606'
    bioc.dump(collection, fout_xml, pretty_print=True)
    fin.close()
    fout_xml.close()
    

#SA for PubTator format
def SA_PubTator(infolder,infile,outpath,nn_model,virus_set,prefix_dict):
    

    # pmid|t|text1
    #pmid|a|text2
    #pmid entity_id sid eid entity_txt entity_type  (gene is blank)
    fin = open(infolder+"/"+infile, 'r',encoding='utf-8')
    # fout_pubtator=open(outpath+'tmp/input_xml.pubtator','w', encoding='utf-8')
    fin_pubtator2=io.StringIO() #two or more species
    all_in_ori=fin.read().strip().split('\n\n')
    fin.close()
    species_gene_count={} #{pmid:{'spec':_species_num;'gene':_gene_num}}
    gene_set=['Gene','FamilyName']
    ML_results={} #{'pmid':{'sid-eid':species_id}}
    
    prefix_speid_allset=set(prefix_dict.keys())
    
    for document in all_in_ori:
        lines=document.split('\n')
        doc_pmid=lines[0].split('|t|')[0]
        doc_title=lines[0].split('|t|')[1]
        doc_abstract=lines[1].split('|a|')[1]
        doc_annotation=[]
        _species_num=set() #(*speciesid)
        _gene_num=0
        _ML_gene_num=0
        _entity_num=0
        _prefix_str2id_dict={} #{prestr:id}
        for i in range(2,len(lines)):
            segs=lines[i].split('\t')
            if segs[4] in gene_set:
                _gene_num+=1
            if len(segs)>=6: #species
                doc_annotation.append(segs[0]+'\t'+str(_entity_num)+'\t'+'\t'.join(segs[1:]))
                species_ID=segs[-1]
                if species_ID.find('*')>=0:
                    _species_num.add(species_ID)
                    if species_ID[1:] in prefix_speid_allset:
                        for ele in prefix_dict[species_ID[1:]]:
                            _prefix_str2id_dict[ele]=species_ID[1:]
            else: #gene
                if segs[3][0:2] in _prefix_str2id_dict:#prefix rule
                    if _prefix_str2id_dict[segs[3][0:2]] in virus_set:
                        doc_annotation.append(segs[0]+'\t'+str(_entity_num)+'\t'+'\t'.join(segs[1:])+'\tFocus:'+_prefix_str2id_dict[segs[3][0:2]]+',9606')
                        if doc_pmid not in ML_results.keys():
                            ML_results[doc_pmid]={segs[1]+'-'+segs[2]:_prefix_str2id_dict[segs[3][0:2]]+',9606'}
                        else:
                            ML_results[doc_pmid][segs[1]+'-'+segs[2]]=_prefix_str2id_dict[segs[3][0:2]]+',9606'

                        # print('!!! prefixr and virus:', doc_pmid)
                    else:
                        doc_annotation.append(segs[0]+'\t'+str(_entity_num)+'\t'+'\t'.join(segs[1:])+'\tFocus:'+_prefix_str2id_dict[segs[3][0:2]])
                        if doc_pmid not in ML_results.keys():
                            ML_results[doc_pmid]={segs[1]+'-'+segs[2]:_prefix_str2id_dict[segs[3][0:2]]}
                        else:
                            ML_results[doc_pmid][segs[1]+'-'+segs[2]]=_prefix_str2id_dict[segs[3][0:2]]
                    # print('prefix rule!!',_prefix_str2id_dict)
                    # print(doc_pmid)
                else:
                    doc_annotation.append(segs[0]+'\t'+str(_entity_num)+'\t'+'\t'.join(segs[1:]))
                    if segs[4] in gene_set:
                        _ML_gene_num+=1
            _entity_num+=1   
                
        if len(_species_num)>=2 and _ML_gene_num>0:
            fin_pubtator2.write(doc_pmid+'|t|'+doc_title+'\n')
            fin_pubtator2.write(doc_pmid+'|a|'+doc_abstract+'\n')
            for ele in doc_annotation:
                fin_pubtator2.write(ele+'\n')
            fin_pubtator2.write('\n')
            
        species_gene_count[doc_pmid]={'spec':_species_num,'gene':_gene_num}
            
    if fin_pubtator2.getvalue()!='':
        #pubtator format  ML tagging
        #print(fin_pubtator2.getvalue())
        ml_out= sa_tag.ml_tag_main(fin_pubtator2,nlp_token, nn_model)   
        #print(ml_out.getvalue())
        fin_result=io.StringIO(ml_out.getvalue())
        all_in=fin_result.read().strip().split('\n\n')
        #print('+2 species:',len(all_in))
        fin_result.close()
        for doc in all_in:
            lines=doc.split('\n')
            pmid=lines[0].split('|t|')[0]

            for i in range(2,len(lines)):
                segs=lines[i].split('\t')
                if pmid not in ML_results.keys():
                    ML_results[pmid]={segs[2]+'-'+segs[3]:segs[-1]}
                else:
                    ML_results[pmid][segs[2]+'-'+segs[3]]=segs[-1]
    
    #output
    fout_pubtator=open(outpath+"/"+infile,'w', encoding='utf8')
    for doc in all_in_ori:
        lines=doc.split('\n')
        pmid=lines[0].split('|t|')[0]
        fout_pubtator.write(lines[0]+'\n'+lines[1]+'\n')
        if len(species_gene_count[pmid]['spec'])>1 and species_gene_count[pmid]['gene']>0: # ML
            for i in range(2,len(lines)):
                segs=lines[i].split('\t')
                if len(segs)>=6: #species
                    fout_pubtator.write(lines[i]+'\n')
                else:#gene
                    if ML_results[pmid][segs[1]+'-'+segs[2]] in virus_set:
                        fout_pubtator.write(lines[i]+'\tFocus:'+ML_results[pmid][segs[1]+'-'+segs[2]]+',9606'+'\n')
                        # print('!!! virus:', pmid)
                    else:
                        fout_pubtator.write(lines[i]+'\tFocus:'+ML_results[pmid][segs[1]+'-'+segs[2]]+'\n')
            fout_pubtator.write('\n')
                    
        elif len(species_gene_count[pmid]['spec'])==1 and species_gene_count[pmid]['gene']>0: #only one species
            for i in range(2,len(lines)):
                segs=lines[i].split('\t')
                if len(segs)>=6: #species
                    fout_pubtator.write(lines[i]+'\n')
                else:#gene
                    major_species,=species_gene_count[pmid]['spec']
                    if major_species[1:] in virus_set:
                        fout_pubtator.write(lines[i]+'\tFocus:'+major_species[1:]+',9606'+'\n')
                        # print('!!! virus:', pmid)
                    fout_pubtator.write(lines[i]+'\tFocus:'+major_species[1:]+'\n')
            fout_pubtator.write('\n')
            
        elif len(species_gene_count[pmid]['spec'])==0 and species_gene_count[pmid]['gene']>0:#no species
            for i in range(2,len(lines)):
                segs=lines[i].split('\t')
                if len(segs)>=6: #species
                    fout_pubtator.write(lines[i]+'\n')
                else:#gene
                    fout_pubtator.write(lines[i]+'\tFocus:9606'+'\n')
            fout_pubtator.write('\n')
            
        else:
            for i in range(2,len(lines)):
                fout_pubtator.write(lines[i]+'\n')
            fout_pubtator.write('\n')   
    fout_pubtator.close()

    
#SA main    
def speciesAss(infolder,outpath, modelfile):

    if modelfile.lower().find('bioformer')>=0:
        model_type='bioformer'
    else:
        model_type='pubmedbert'
    
    print('loading SA models........')     
    if model_type=='bioformer':
        
        vocabfiles={'labelfile':'./vocab/SpeAss_IO_label.vocab',
                    'checkpoint_path':'./gnorm_trained_models/bioformer-cased-v1.0/',
                    'lowercase':False,
                    }
    else:
        vocabfiles={'labelfile':'./vocab/SpeAss_IO_label.vocab',
                    'checkpoint_path':'./gnorm_trained_models/BiomedNLP-PubMedBERT-base-uncased-abstract/',
                    'lowercase':True,
                    }
        
    nn_model=model_sa.HUGFACE_NER(vocabfiles)
    nn_model.build_encoder()
    nn_model.build_softmax_decoder()
    nn_model.load_model(modelfile)
    
    dict_filename={'prefix':'./Dictionary/SPPrefix.txt',
                   'virus':'./Dictionary/SP_Virus2HumanList.txt'}
    fin=open(dict_filename['virus'],'r',encoding='utf-8')
    virus_set=set(fin.read().strip().split('\n'))
    fin.close()
    
    prefix_dict={}#{id:[prefix1,prefix2]}
    fin=open(dict_filename['prefix'],'r',encoding='utf-8')
    for line in fin:
        seg= line.strip().split('\t')
        if seg[0] not in prefix_dict.keys():
            prefix_dict[seg[0]]=seg[1].split('|')
        else:
            prefix_dict[seg[0]].extend(seg[1].split('|'))
    fin.close()
    

    
    print("begin species assignment........")
    start_time=time.time() 
    
    for infile in os.listdir(infolder): 
        if os.path.isfile(outpath+"/"+infile):
            print(infile+' has exsited.')
        else:
            print('Processing:',infile)
            fin=open(infolder+"/"+infile, 'r',encoding='utf-8')
            file_format=""
            for line in fin:
                pattern_bioc = re.compile('.*<collection>.*')
                pattern_pubtator = re.compile('^([^\|]+)\|[^\|]+\|(.*)')
                if pattern_pubtator.search(line):
                    file_format="PubTator"
                    break
                elif pattern_bioc.search(line):
                    file_format="BioC"
                    break
            fin.close()
            if(file_format == "PubTator"):
                SA_PubTator(infolder,infile,outpath,nn_model,virus_set,prefix_dict)
            elif(file_format == "BioC"):
                SA_BioC(infolder,infile,outpath,nn_model,virus_set,prefix_dict)
        

    print('species assignment done:',time.time()-start_time)

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='run GeneNER and species assignment, python GeneNER_SpeAss_run.py -i input -n NERmodel -s SAmodel -r neroutput -a saoutput')
    parser.add_argument('--infolder', '-i', help="input folder",default='./example/input/')
    parser.add_argument('--NERmodel', '-n', help="trained deep learning NER model file",default='')
    parser.add_argument('--SAmodel', '-s', help="trained deep learning species assignment model file",default='')
    parser.add_argument('--NERoutpath', '-r', help="output folder to save the NER tagged results",default='./example/ner_output/')
    parser.add_argument('--SAoutpath', '-a', help="output folder to save the SA tagged results",default='./example/sa_output/')
    parser.add_argument('--NUM_THREADS', '-t', help="Number of threads",default='3')
    args = parser.parse_args()
       
    
    if args.NUM_THREADS.isdigit() == False:
        args.NUM_THREADS='3'
        
    tf.config.threading.set_inter_op_parallelism_threads(int(args.NUM_THREADS))
    tf.config.threading.set_intra_op_parallelism_threads(int(args.NUM_THREADS))

    if args.NERmodel!='' and args.SAmodel!='':
        
        #pipleline
        print('==============\n| GeneNER and SpeAss |\n==============')
        
        #creat output folder
        
        if args.infolder[-1]!='/':
            args.infolder+='/'
        if not os.path.exists(args.infolder):
            os.makedirs(args.infolder)
            
        if args.NERoutpath[-1]!='/':
            args.NERoutpath+='/'
        if not os.path.exists(args.NERoutpath):
            os.makedirs(args.NERoutpath)
            
        if args.SAoutpath[-1]!='/':
            args.SAoutpath+='/'
        if not os.path.exists(args.SAoutpath):
            os.makedirs(args.SAoutpath)
        
        #1. gene NER, the results are saved in outpath/ner_tmp/
        geneNER(args.infolder,args.NERoutpath, args.NERmodel)
        
        
        #2. species assignment, the results are saved in outpath/sa_tmp/
        speciesAss(args.NERoutpath,args.SAoutpath, args.SAmodel)
        
    elif args.NERmodel!='' and args.SAmodel=='':
        if args.infolder[-1]!='/':
            args.infolder+='/'
        if not os.path.exists(args.infolder):
            os.makedirs(args.infolder)
            
        # only geneNER
        if args.NERoutpath[-1]!='/':
            args.NERoutpath+='/'
        if not os.path.exists(args.NERoutpath):
            os.makedirs(args.NERoutpath)
            
        print('==============\n| GeneNER |\n==============')
        geneNER(args.infolder,args.NERoutpath,args.NERmodel)
        
    elif args.NERmodel=='' and args.SAmodel!='':
        # only speass
        if args.SAoutpath[-1]!='/':
            args.SAoutpath+='/'
        if not os.path.exists(args.SAoutpath):
            os.makedirs(args.SAoutpath)
            
        print('==============\n| SpeAss |\n==============')
        speciesAss(args.infolder,args.SAoutpath,args.SAmodel)
    else:
        print('Please provide models!')
        
     
    