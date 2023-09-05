# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 11:52:42 2022

@author: luol2
"""

import io
import os
import argparse
import stanza
import sys
import re
import bioc
from src_python.SpeAss.ml_tagging_score_sa import NER_Tag


def ssplit_token(infile,nlp_token):
    fin=io.StringIO(infile.getvalue())
    fout=io.StringIO()
    # fout=open(outfile,'w',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    ori_text_newentity={} #{line[0]+line[1]:[all entity]}
    entity_type=set()
    token_text_new={}#{pmid:token_text}
    for doc_text in all_in:
        lines=doc_text.split('\n')
        ori_text=lines[0].split('|t|')[1]+' '+lines[1].split('|a|')[1]
        pmid=lines[0].split('|t|')[0]
        # print(pmid)
        entity_all=[]   #[[seg0,seg1,...,],[]]
        entity_all_ori=[]
        entity_num=0
        
        #first sort
        doc_result={}
        for i in range(2,len(lines)):
            segs=lines[i].split('\t')
            doc_result[lines[i]]=[int(segs[2]),int(segs[3])]
        doc_result=sorted(doc_result.items(), key=lambda kv:(kv[1]), reverse=False)
        doc_result_sort=[]
        for ele in doc_result:
            doc_result_sort.append(ele[0])
        
        for i in range(0,len(doc_result_sort)):
            seg=doc_result_sort[i].strip().split('\t')
            entity_type.add(seg[5])
            # print(seg)
            if len(seg)<=6:#Gene
                entity_all_ori.append([seg[0],seg[1],'M'+str(entity_num),seg[2],seg[3],seg[4],seg[5],'-'])
                entity_all.append([seg[0],seg[1],'M'+str(entity_num),seg[2],seg[3],seg[4],'Gene','-'])
                entity_num+=1
            elif seg[-1].find('*')>=0:# *Species
                entity_all_ori.append([seg[0],seg[1],'M'+str(entity_num),seg[2],seg[3],seg[4],seg[5],seg[6]])
                entity_all.append([seg[0],seg[1],'M'+str(entity_num),seg[2],seg[3],seg[4],'Species',seg[6]])
                entity_num+=1
        ori_text_newentity[lines[0]+'\n'+lines[1]]=entity_all_ori
        # sys.exit()

        #ssplit token
        doc_stanza = nlp_token(ori_text)
        token_text=''
        sentence_index=[] #[text_offset]
        for sent in doc_stanza.sentences:
            for word in sent.words:
                if word.text.strip()=='':
                    # print('token is blank!')
                    pass
                token_text+=word.text+' '
            token_text=token_text+''  #sentence split 
            sentence_index.append(len(token_text))
        
        #ori_index map token_index
        index_map=[-1]*len(ori_text)
        j=0
        space_list=[' ',chr(160),chr(8201),chr(8194),chr(8197),chr(8202)] #空格有好几种，第一个是常用32,第二个shi 160,8201,8194,8197
        for i in range(0,len(ori_text)):
            if ori_text[i] in space_list:
                pass
            elif ori_text[i]==token_text[j]:
                index_map[i]=j
                j+=1
            else:
                j+=1
                temp_log=j
                try:
                    while(ori_text[i]!=token_text[j]):
                        j+=1
                except:
                    print('doc',doc_text)
                    print('token_text:',token_text)
                    print('error:',ori_text[i-10:i+10],'i:',ori_text[i],'j:',token_text[temp_log],',',token_text[temp_log-10:temp_log+10])
                    print(ord(ori_text[i]),ord(' '))
                    sys.exit()
                index_map[i]=j
                j+=1
        # token_text=token_text.replace('     ','<EOS>')
        # print(token_text)
        fout.write(token_text+'\n')
        token_text_new[pmid]=token_text
        entity_i=0
        cur_sent_i=0
        new_ente=0
        cur_sents=0
        cur_sente=sentence_index[0]
        if entity_all!=[]:
            bug_new_entity=[]
            for entity_i in range(0,len(entity_all)):
                new_ents=index_map[int(entity_all[entity_i][3])]
                new_ente=index_map[int(entity_all[entity_i][4])-1]+1
                new_ent=token_text[new_ents:new_ente]
                old_ent=entity_all[entity_i][5]
                cur_sent_i=0
                cur_sents=0
                cur_sente=sentence_index[0]
                while (not (max(new_ents,cur_sents) <= min(new_ente,cur_sente))) and (cur_sent_i<len(sentence_index)-1):
                    cur_sent_i+=1
                    cur_sents=sentence_index[cur_sent_i-1]
                    cur_sente=sentence_index[cur_sent_i]
                fout.write(entity_all[entity_i][0]+'\t'+entity_all[entity_i][1]+'\t'+entity_all[entity_i][2]+'-'+str(cur_sent_i)+'\t'+str(new_ents)+'\t'+str(new_ente)+'\t'+new_ent+'\t'+entity_all[entity_i][6]+'\t'+entity_all[entity_i][7]+'\n')
                bug_new_entity.append(entity_all[entity_i][0]+'\t'+entity_all[entity_i][1]+'\t'+entity_all[entity_i][2]+'-'+str(cur_sent_i)+'\t'+str(new_ents)+'\t'+str(new_ente)+'\t'+new_ent+'\t'+entity_all[entity_i][6]+'\t'+entity_all[entity_i][7]+'\n')

            """
            new_ents=index_map[int(entity_all[entity_i][3])]
            new_ente=index_map[int(entity_all[entity_i][4])-1]+1
            new_ent=token_text[new_ents:new_ente]
            old_ent=entity_all[entity_i][5]
            while True:
                while new_ents>=cur_sents and new_ente< cur_sente:
    
                    if new_ent.replace(' ','') !=old_ent.replace(' ',''):
                        # print('entity error:',pmid,old_ent,new_ent,entity_all[entity_i][2],entity_all[entity_i][3])
                        pass
                    fout.write(entity_all[entity_i][0]+'\t'+entity_all[entity_i][1]+'\t'+entity_all[entity_i][2]+'-'+str(cur_sent_i)+'\t'+str(new_ents)+'\t'+str(new_ente)+'\t'+new_ent+'\t'+entity_all[entity_i][6]+'\t'+entity_all[entity_i][7]+'\n')
                    entity_i+=1
                    if entity_i>=len(entity_all):
                        break
                    new_ents=index_map[int(entity_all[entity_i][3])]
                    new_ente=index_map[int(entity_all[entity_i][4])-1]+1
                    new_ent=token_text[new_ents:new_ente]
                    old_ent=entity_all[entity_i][5]
                cur_sent_i+=1
                if cur_sent_i >= len(sentence_index):
                    break
                cur_sents=sentence_index[cur_sent_i-1]
                cur_sente=sentence_index[cur_sent_i]
            """
        fout.write('\n')
    # print(entity_type)
    # fout.close()
    return ori_text_newentity,token_text_new,fout

def filter_nest(infile): #nonest

    # fin=open(infile,'r',encoding='utf-8')
    # fout=open(outfile,'w',encoding='utf-8')
    fin=io.StringIO(infile.getvalue())
    fout=io.StringIO()
    
    documents=fin.read().strip().split('\n\n')
    fin.close()
    total_entity=0
    over_entity=0
    nest_entity=0
    for doc in documents:
        lines=doc.split('\n')
        context=lines[0]
        entity_list=[]
        if len(lines)>1:
            first_entity=lines[1].split('\t')
            nest_list=[first_entity]
            max_eid=int(first_entity[4])
            total_entity+=len(lines)-2
            for i in range(2,len(lines)):
                segs=lines[i].split('\t')
                if int(segs[3])> max_eid:
                    if len(nest_list)==1:
                        entity_list.append(nest_list[0])
                        nest_list=[]
                        nest_list.append(segs)
                        if int(segs[4])>max_eid:
                            max_eid=int(segs[4])
                    else:
                        # print(nest_list)
                        nest_entity+=len(nest_list)-1
                        tem=find_max_entity(nest_list)#find max entity
                        # if len(tem)>1:
                        #     print('max nest >1:',tem)
                        entity_list.extend(tem)
                        nest_list=[]
                        nest_list.append(segs)
                        if int(segs[4])>max_eid:
                            max_eid=int(segs[4])
                        
                else:
                    nest_list.append(segs)
                    if int(segs[4])>max_eid:
                        max_eid=int(segs[4])
            if nest_list!=[]:
                if len(nest_list)==1:
                    entity_list.append(nest_list[0])

                else:
                    tem=find_max_entity(nest_list)#find max entity
                    # if len(tem)>1:
                    #     print('max nest >1:',tem)
                    entity_list.extend(tem)
        fout.write(context+'\n')
        for ele in entity_list:
            fout.write('\t'.join(ele)+'\n')
        fout.write('\n')
    # print(total_entity,over_entity, nest_entity)
    return fout
def find_max_entity(nest_list):
    max_len=0
    final_tem=[]
    max_index=0
    for i in range(0, len(nest_list)):
        cur_len=int(nest_list[i][4])-int(nest_list[i][3])
        if cur_len>max_len:
            max_len=cur_len
            max_index=i
        elif cur_len==max_len:
            if nest_list[i][6] =='Gene':
                max_index=i
        # elif nest_list[i][5] =='Species':
        #     final_tem.append(nest_list[i])

    final_tem.append(nest_list[max_index])
    return final_tem


# machine learning species assignment
def ml_tag(infile,nn_model):
    
    #tagging text
    fin=io.StringIO(infile.getvalue())
    fout=io.StringIO()
    # fin=open(infile,'r',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()

    for doc in all_in:            
        pre_result,entity_all=NER_Tag(doc, nn_model)
        for ele in entity_all:
            ent_id=ele[3]+'-'+ele[4]
            if ent_id in pre_result.keys():
                fout.write('\t'.join(ele)+'\t'+','.join(pre_result[ent_id])+'\n')
            else:
                fout.write('\t'.join(ele)+'\t-\n')
        fout.write('\n')
    
    return fout


# details nearest+and
def post_rule1(ori_context,token_text,infile,outfile):
    fin=open(infile,'r',encoding='utf-8')
    fout=open(outfile,'w',encoding='utf-8')
    pred_results={} #{pmid:{'M0':{'sent':'','offset':[sid,eid],'score':[[id,score],[id,score]]}}} #gene
    species_index={} #{pmid:{sentid:[[spe_seg1],[spe_seg]]}}
    mem_sent={} #{pmid:{'M0':sentid}}
    gene_num=0
    gene_none=0
    for line in fin:
        seg=line.strip().split('\t')
        if len(seg)>1:
            if seg[0] not in mem_sent.keys():
                _term_seg=seg[1].split('-')
                mem_sent[seg[0]]={_term_seg[0]:_term_seg[1]}
            else:
                _term_seg=seg[1].split('-')
                mem_sent[seg[0]][_term_seg[0]]=_term_seg[1]
            if seg[5]=='Species':   
                if seg[0] not in species_index.keys():
                    _sent_id=seg[1].split('-')[1]
                    species_index[seg[0]]={_sent_id:[seg]}
                else:
                    _sent_id=seg[1].split('-')[1]
                    if _sent_id in species_index[seg[0]].keys():
                        species_index[seg[0]][_sent_id].append(seg)
                    else:
                        species_index[seg[0]][_sent_id]=[seg]
            else:
                _pred_ids=seg[-1].split(',')
                _temp_id_score=[] #[[spe_id,score]]
                _sent_id=seg[1].split('-')[1]
                for _pred_id in _pred_ids:
                    _temp_id_score.append(_pred_id.split('|'))
                if seg[0] not in pred_results.keys():                   
                    pred_results[seg[0]]={seg[1].split('-')[0]:{'sent':_sent_id,'offset':[seg[2],seg[3]],'score':_temp_id_score}}
                else:
                    pred_results[seg[0]][seg[1].split('-')[0]]={'sent':_sent_id,'offset':[seg[2],seg[3]],'score':_temp_id_score}
    #print(pred_results)
    for pmid_text in ori_context.keys():
        #print(pmid_text)
        lines=pmid_text.split('\n')
        ori_text=lines[0].split('|t|')[1]+' '+lines[1].split('|a|')[1]
        # print(ori_text)
        fout.write(pmid_text+'\n')
        pmid=lines[0].split('|t|')[0]
        before_species=[] #nearest [eid,spe_id]
        after_species=[] #nearest  [sid,spe_id]
        doc_specs=species_index[pmid]
        #mul and spe
        mul_and_spe=[]
        for spe_sent in doc_specs.keys():
            last_id=''
            new_diff_spe=[]
            _temp_speid=set()
            for ele in doc_specs[spe_sent]:
                if ele[-2] !=last_id:
                    new_diff_spe.append(ele)
                    last_id =ele[-2]
                    _temp_speid.add(ele[-2])
                else:
                    
                    new_diff_spe.pop()
                    new_diff_spe.append(ele)
                    _temp_speid.add(ele[-2])
                    last_id =ele[-2]
            if len(new_diff_spe)==2:
                spe_and_text=new_diff_spe[0][4]+' and '+new_diff_spe[1][4]
                if ori_text.find(spe_and_text)>=0:
                    # print('old:',doc_specs[spe_sent])
                    # print('new:',new_diff_spe)
                    # print('\n')
                    mul_and_spe=list(_temp_speid)
                    # print(mul_and_spe)
            elif len(new_diff_spe)>2:
                spe_and_text=''
                for i in range(0,len(new_diff_spe)-1):
                    spe_and_text+=new_diff_spe[i][4]+', '
                spe_and_text1=spe_and_text[0:-2]+' and '+new_diff_spe[-1][4]
                spe_and_text2=spe_and_text+'and '+new_diff_spe[-1][4]
                if ori_text.find(spe_and_text1)>=0 or ori_text.find(spe_and_text2)>=0:
                    # print('old:',doc_specs[spe_sent])
                    # print('new:',new_diff_spe)
                    mul_and_spe=list(_temp_speid)
                    # print(mul_and_spe)
        Gene_type_list=['Gene','FamilyName','DomainMotif']
        for i,ele in enumerate(ori_context[pmid_text]):
            #print(ele)
            
            if ele[5] in Gene_type_list:
                gene_num+=1
                final_preds=set()
                if ele[1] in pred_results[ele[0]].keys():
                    temp_preds=pred_results[ele[0]][ele[1]]['score']
                    
                    if temp_preds!=[['-']]:
                        if len(temp_preds)==1:
                            final_preds.add(temp_preds[0][0])
                        else:
                            max_id=''
                            max_score=0
                            for _temp_pred in temp_preds:
                                _score=float(_temp_pred[1])
                                _id_ass=_temp_pred[0]
                                if len(mul_and_spe)>1:
                                    if _score>0.5 and (_id_ass in mul_and_spe):
                                        # print(_score)
                                        final_preds.add(_id_ass)
                                if _score>max_score:
                                    max_id=_id_ass
                                    max_score=_score
                            if len(final_preds)==0:
                                final_preds.add(max_id) 
                            # final_preds.add(multi_id)
                    else: #'-' nearst rule
                        gene_none+=1
                        # print(mem_sent[ele[0]])
                        _sent_id_gene=mem_sent[ele[0]][ele[1]]
                     
                        for j in range(i+1,len(ori_context[pmid_text])):
                            temp_seg=ori_context[pmid_text][j]
                            if temp_seg[5]=='Species':
                                after_species=[int(temp_seg[2]),temp_seg[6]]
                                break
                        # print(before_species,after_species)
                        # print(seg)
                        if before_species!=[] and after_species!=[]:
                            if len(ori_text[before_species[0]:int(ele[2])].split()) > len(ori_text[int(ele[3]):after_species[0]].split()):
                                final_preds.add(after_species[1])
                            else:
                                final_preds.add(before_species[1])
                        elif before_species==[]:
                            final_preds.add(after_species[1])
                        elif after_species==[]:
                            final_preds.add(before_species[1])
                    if len(final_preds)==0:
                        print('none pred!!!')                             
                    fout.write(ele[0]+'\t'+'\t'.join(ele[2:])+'\t'+','.join(final_preds)+'\n')
                else:
                    # gene_none+=1
                    # print(ele)
                    for j in range(i+1,len(ori_context[pmid_text])):
                        temp_seg=ori_context[pmid_text][j]
                        if temp_seg[5]=='Species':
                            after_species=[int(temp_seg[2]),temp_seg[6]]
                            break
                    # print(before_species,after_species)
                    # print(seg)
                    if before_species!=[] and after_species!=[]:
                        if len(ori_text[before_species[0]:int(ele[2])].split()) > len(ori_text[int(ele[3]):after_species[0]].split()):
                            final_preds.add(after_species[1])
                        else:
                            final_preds.add(before_species[1])
                    elif before_species==[]:
                        final_preds.add(after_species[1])
                    elif after_species==[]:
                        final_preds.add(before_species[1])
                    fout.write(ele[0]+'\t'+'\t'.join(ele[2:])+'\t'+','.join(final_preds)+'\n')
            else:
                fout.write(ele[0]+'\t'+'\t'.join(ele[2:])+'\t-\n')
                before_species=[int(ele[3]),ele[6]]
        fout.write('\n')
    print('gene, none:',gene_num,gene_none)
    fout.close()
    
# major+and
def post_rule2(ori_context,token_text,infile):
    # fin=open(infile,'r',encoding='utf-8')
    # fout=open(outfile,'w',encoding='utf-8')
    fin=io.StringIO(infile.getvalue())
    fout=io.StringIO()
    pred_results={} #{pmid:{'M0':{'sent':'','offset':[sid,eid],'score':[[id,score],[id,score]]}}} #gene
    species_index={} #{pmid:{sentid:[[spe_seg1],[spe_seg]]}}
    species_count={}#{pmid:{speid:num}}
    gene_num=0
    gene_none=0
    for line in fin:
        seg=line.strip().split('\t')
        if len(seg)>1:
            if seg[6]=='Species':
                if seg[0] not in species_count.keys():
                    species_count[seg[0]]={seg[-2]:1}
                else:
                    if seg[-2] not in species_count[seg[0]].keys():
                        species_count[seg[0]][seg[-2]]=1
                    else:
                        species_count[seg[0]][seg[-2]]+=1
                        
                if seg[0] not in species_index.keys():
                    _sent_id=seg[2].split('-')[1]
                    species_index[seg[0]]={_sent_id:[seg]}
                else:
                    _sent_id=seg[2].split('-')[1]
                    if _sent_id in species_index[seg[0]].keys():
                        species_index[seg[0]][_sent_id].append(seg)
                    else:
                        species_index[seg[0]][_sent_id]=[seg]
            else:
                _pred_ids=seg[-1].split(',')
                _temp_id_score=[] #[[spe_id,score]]
                _sent_id=seg[2].split('-')[1]
                for _pred_id in _pred_ids:
                    _temp_id_score.append(_pred_id.split('|'))
                if seg[0] not in pred_results.keys():                   
                    pred_results[seg[0]]={seg[2].split('-')[0]:{'sent':_sent_id,'offset':[seg[3],seg[4]],'score':_temp_id_score}}
                else:
                    pred_results[seg[0]][seg[2].split('-')[0]]={'sent':_sent_id,'offset':[seg[3],seg[4]],'score':_temp_id_score}
    fin.close()
    #print(pred_results)
    for pmid_text in ori_context.keys():
        #print(pmid_text)
        lines=pmid_text.split('\n')
        ori_text=lines[0].split('|t|')[1]+' '+lines[1].split('|a|')[1]
        # print(ori_text)
        fout.write(pmid_text+'\n')
        pmid=lines[0].split('|t|')[0]
        if pmid in species_count.keys():
            marjor_species = max(zip(species_count[pmid].values(), species_count[pmid].keys()))
        else:
            marjor_species = (1000,'*9606')
        
        if pmid in species_index.keys():
            doc_specs=species_index[pmid]
            #mul and spe
            mul_and_spe=[]
            for spe_sent in doc_specs.keys():
                last_id=''
                new_diff_spe=[]
                _temp_speid=set()
                for ele in doc_specs[spe_sent]:
                    if ele[-2] !=last_id:
                        new_diff_spe.append(ele)
                        last_id =ele[-2]
                        _temp_speid.add(ele[-2])
                    else:
                        
                        new_diff_spe.pop()
                        new_diff_spe.append(ele)
                        _temp_speid.add(ele[-2])
                        last_id =ele[-2]
                if len(new_diff_spe)==2:
                    spe_and_text=new_diff_spe[0][5]+' and '+new_diff_spe[1][5]
                    if ori_text.find(spe_and_text)>=0:
                        # print('old:',doc_specs[spe_sent])
                        # print('new:',new_diff_spe)
                        # print('\n')
                        mul_and_spe=list(_temp_speid)
                        # print(mul_and_spe)
                elif len(new_diff_spe)>2:
                    spe_and_text=''
                    for i in range(0,len(new_diff_spe)-1):
                        spe_and_text+=new_diff_spe[i][5]+', '
                    spe_and_text1=spe_and_text[0:-2]+' and '+new_diff_spe[-1][5]
                    spe_and_text2=spe_and_text+'and '+new_diff_spe[-1][5]
                    if ori_text.find(spe_and_text1)>=0 or ori_text.find(spe_and_text2)>=0:
                        # print('old:',doc_specs[spe_sent])
                        # print('new:',new_diff_spe)
                        mul_and_spe=list(_temp_speid)
                        # print(mul_and_spe)
        else:
            mul_and_spe=[]
            
        Gene_type_list=['Gene','FamilyName','DomainMotif']
        for i,ele in enumerate(ori_context[pmid_text]):
            #print(ele)
            
            if ele[6] in Gene_type_list:
                gene_num+=1
                final_preds=set()
                if (ele[0] in pred_results.keys()) and (ele[2] in pred_results[ele[0]].keys()):
                    temp_preds=pred_results[ele[0]][ele[2]]['score']
                    
                    if temp_preds!=[['-']]:
                        if len(temp_preds)==1:
                            final_preds.add(temp_preds[0][0])
                        else:
                            max_id=''
                            max_score=0
                            for _temp_pred in temp_preds:
                                _score=float(_temp_pred[1])
                                _id_ass=_temp_pred[0]
                                if len(mul_and_spe)>1:
                                    if _score>0.5 and (_id_ass in mul_and_spe):
                                        # print(_score)
                                        final_preds.add(_id_ass)
                                if _score>max_score:
                                    max_id=_id_ass
                                    max_score=_score
                            if len(final_preds)==0:
                                final_preds.add(max_id) 
                            # final_preds.add(multi_id)
                    else: #'-' major species
                        gene_none+=1
                        # print(mem_sent[ele[0]])
                        final_preds.add(marjor_species[1])
                    if len(final_preds)==0:
                        print('none pred!!!')                             
                    fout.write(ele[0]+'\t'+ele[1]+'\t'+'\t'.join(ele[3:-1])+'\t'+','.join(final_preds).replace('*','')+'\n')
                else:
                    final_preds.add(marjor_species[1])
                    fout.write(ele[0]+'\t'+ele[1]+'\t'+'\t'.join(ele[3:-1])+'\t'+','.join(final_preds).replace('*','')+'\n')
            else:
                fout.write(ele[0]+'\t'+ele[1]+'\t'+'\t'.join(ele[3:])+'\n')
        fout.write('\n')
    # print('gene, none:',gene_num,gene_none)
    return fout

def ml_tag_main(fin_pubtator,nlp_token, nn_model):
    #print('.......senten split, tokenizer.........')
    #print('...in...\n',fin_pubtator.getvalue())
    ori_text_newentity,token_text,token_out=ssplit_token(fin_pubtator,nlp_token)
    #print('...token....\n',token_out.getvalue())

    #2. filter nest entity
    nonest_out=filter_nest(token_out)
    #print(nonest_out.getvalue())

    #3.ml tag
    #print('.......machine learning-based tagging.........')
    
    ml_out=ml_tag(nonest_out, nn_model)
    #print('.....ml.....\n',ml_out.getvalue())
    
    #4. post processing
    #print('.......post processing.........')
    post_out=post_rule2(ori_text_newentity,token_text,ml_out)
    #print('.........ori_text...............\n', ori_text_newentity)
    #print('.....post.....\n',post_out.getvalue())
    return post_out
    