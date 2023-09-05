# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 09:29:46 2022

@author: luol2

machine learning tagging

"""


import time
import io

from src_python.SpeAss.processing_data_sa import ml_intext_fn,out_BIO_BERT_softmax_score_fn
import tensorflow as tf
gpu = tf.config.list_physical_devices('GPU')
print("Num GPUs Available: ", len(gpu))
if len(gpu) > 0:
    tf.config.experimental.set_memory_growth(gpu[0], True)
#tf.compat.v1.disable_eager_execution()

REL_ENT={'arg1':'Species',
         'arg2':'Gene'}

entity_tag={'arg1':['arg1s','arg1e'],
            'gene':['gene1s','gene1e'],
            'species':['species1s','species1e']
    }

def input_preprocess_notoken(doc_text):
    final_input=[]
    final_id=[]
    
    lines=doc_text.split('\n')
    token_text=lines[0]
    pmid=lines[1].split('\t')[0]
    entity_arg1={} #{species_id:[[spe_sid1,sep_eid1],[...]]}
    entity_all=[]
    for i in range(1,len(lines)):
        seg=lines[i].split('\t')
        if seg[6]==REL_ENT['arg1']:
            if seg[-1] in entity_arg1.keys():
                entity_arg1[seg[-1]].append([seg[3],seg[4]])
            else:
                entity_arg1[seg[-1]]=[[seg[3],seg[4]]]
        entity_all.append(seg)
    
    #print(token_text)
    #print(entity_chemical)
    #generate input instance
    for cur_ele in entity_arg1:
        
        #2. ner label text
        ner_text=''
        text_sid=0
        #print('nonest:',entity_nonest)
        for ele_nonest in entity_all:
            ent_id=[ele_nonest[3],ele_nonest[4]]
            ent_spe_id=ele_nonest[-1]
            ent_sid=int(ele_nonest[3])
            ent_eid=int(ele_nonest[4])
            # print('sid,eid:',ent_sid,ent_eid)
            ent_text=ele_nonest[5]
            ent_type=ele_nonest[6]
            if ent_sid>=text_sid:
                # if token_text[ent_sid:ent_eid]!=ent_text:
                #     print('error!index_text,entext:',token_text[ent_sid:ent_eid],ent_text)
                if ent_id in entity_arg1[cur_ele]: #is species    
                    ner_text+=token_text[text_sid:ent_sid]+' '+ent_spe_id+'|'+entity_tag['arg1'][0]+' '+ent_text+' '+entity_tag['arg1'][1]+' '
                else:
                    ner_text+=token_text[text_sid:ent_sid]+' '+str(ent_sid)+'-'+str(ent_eid)+'|'+entity_tag[ent_type.lower()][0]+' '+ent_text+' '+entity_tag[ent_type.lower()][1]+' '
                text_sid=ent_eid                                      
        ner_text+=token_text[text_sid:]
        sen_tokens=ner_text.split()
        #print('\nner_text:',ner_text)
        
        #3. produce input
        temp_input=[]
        temp_id={'species':'','gene':[]}
        for sen_token in sen_tokens:
            if sen_token.find(entity_tag['arg1'][0])>=0:
                en_id=sen_token.split('|')[0]
                temp_id['species']=en_id
                temp_input.append(entity_tag['arg1'][0]+'\tO')
            elif sen_token.find(entity_tag['gene'][0])>=0:
                en_id=sen_token.split('|')[0]
                temp_id['gene'].append(en_id)
                temp_input.append(entity_tag['gene'][0]+'\tO')
            elif sen_token.find(entity_tag['species'][0])>=0:
                en_id=sen_token.split('|')[0]
                # temp_id.append(en_id)
                temp_input.append(entity_tag['species'][0]+'\tO')
            else:
                if sen_token=='':
                    # print('token is none!error!')
                    pass
                else:
                    temp_input.append(sen_token+'\tO')
        final_input.append('\n'.join(temp_input))
        final_id.append(temp_id)

        # print(entity_nonest)
    return final_input,final_id,entity_all,pmid


def ml_tagging(ml_input,nn_model):

    test_set,test_label = ml_intext_fn(ml_input)
    test_x,test_y, test_bert_text_label=nn_model.rep.load_data_hugface(test_set,test_label,word_max_len=nn_model.maxlen,label_type='softmax')
    test_pre = nn_model.model.predict(test_x)
    ml_out=out_BIO_BERT_softmax_score_fn(test_pre,test_bert_text_label,nn_model.rep.index_2_label)
    return ml_out

def output_rel(ml_output,entity_map,pmid):
    fin=io.StringIO(ml_output)
    alltexts=fin.read().strip().split('\n\n')
    fin.close()
    final_out={} #{'sid-eid':[spechies id]}
    for sen_id,sentence in enumerate(alltexts):
        tokens=sentence.split('\n')
        gene_entity_id=0
        token_id=0
        arg1=''
        arg2_list=[] #[[ID, score],[id,score]]
        while (token_id<len(tokens)):
            seg=tokens[token_id].split('\t')
            if seg[0]==entity_tag['arg1'][0]:
                arg1=entity_map[sen_id]['species']
                token_id+=1
                if token_id >=len(tokens):
                    break
                seg=tokens[token_id].split('\t')
                while seg[0]!=entity_tag['arg1'][1]:
                    token_id+=1
                    if token_id >=len(tokens):
                        break
                    seg=tokens[token_id].split('\t')
            elif seg[0]==entity_tag[REL_ENT['arg2'].lower()][0]:
                temp_rel=seg[-2]
                temp_score=seg[-1]
                arg2_id=entity_map[sen_id]['gene'][gene_entity_id]
                gene_entity_id+=1
                token_id+=1
                if token_id >=len(tokens):
                    break
                seg=tokens[token_id].split('\t')
                while seg[0]!=entity_tag[REL_ENT['arg2'].lower()][1]:
                    token_id+=1
                    if token_id >=len(tokens):
                        break
                    seg=tokens[token_id].split('\t')
                    if seg[-2].find('ARG2')>=0 and temp_rel.find('ARG2')<0:
                        temp_rel=seg[-2]
                        temp_score=seg[-1]
                if temp_rel.find('ARG2')>=0:
                    arg2_list.append([arg2_id,temp_score])
            elif seg[0]==entity_tag[REL_ENT['arg1'].lower()][0]:
                token_id+=1
                if token_id >=len(tokens):
                    break
                seg=tokens[token_id].split('\t')
                while seg[0]!=entity_tag[REL_ENT['arg1'].lower()][1]:
                    token_id+=1
                    if token_id >=len(tokens):
                        break
                    seg=tokens[token_id].split('\t')

            else:
                pass
            token_id+=1
        #print(arg1,arg2_list)
        if arg2_list!=[] and arg1!='':
            for arg2_ele in arg2_list:
                if arg2_ele[0] not in final_out.keys():
                    final_out[arg2_ele[0]]=[arg1+'|'+arg2_ele[1]]
                else:
                    final_out[arg2_ele[0]].append(arg1+'|'+arg2_ele[1])
    return(final_out)

def NER_Tag(doc_in,nn_model):
    
    #1. preprocess input, input_text:conll格式, input_entity：相应的实体列表
    #print(doc_in)
    input_text,entity_index,entity_all,pmid=input_preprocess_notoken(doc_in)
    # print('pmid:',pmid)
    # print('\entity_index:',entity_index)

   
    #2. ml tagging
    if input_text!=[]:
        ml_pre=ml_tagging(input_text,nn_model)
        #print('\noutput:')
        #print(ml_pre)

        #3.generate output
        final_output=output_rel(ml_pre,entity_index,pmid)
    else:
        final_output={}
    return final_output,entity_all
    
    


    
    














