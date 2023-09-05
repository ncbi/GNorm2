# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 08:58:22 2022

@author: luol2
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 11:27:57 2022

@author: luol2
"""


import stanza
import sys
import os
import argparse
import io
import json
import re
#sort entity by position in text
def pubtator_entitysort(infile):
    
    fin=open(infile,'r',encoding='utf-8')    
    # fout=open(path+'LitCoin/sort/Train_sort.PubTator','w',encoding='utf-8')
    fout=io.StringIO()
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    error_dict={} #use to debug error
    for doc in all_in:
        entity_dict={}
        lines=doc.split('\n')
        fout.write(lines[0]+'\n'+lines[1]+'\n')
        for i in range(2,len(lines)):
            segs=lines[i].split('\t')
            if len(segs)>=5:
                if lines[i] not in entity_dict.keys():
                    entity_dict[lines[i]]=int(segs[1])
                else:
                    print('entity have in',lines[i])
                    if segs[0] not in error_dict.keys():
                        error_dict[segs[0]]=[lines[i]]
                    else:
                        if lines[i] not in error_dict[segs[0]]:
                            error_dict[segs[0]].append(lines[i])

        entity_sort=sorted(entity_dict.items(), key=lambda kv:(kv[1]), reverse=False)
        for ele in entity_sort:
            fout.write(ele[0]+'\n')
        fout.write('\n')
    return fout

def filter_overlap(infile): #nonest

    fin=io.StringIO(infile.getvalue())
    fout=io.StringIO()
    
    documents=fin.read().strip().split('\n\n')
    fin.close()
    total_entity=0
    over_entity=0
    nest_entity=0
    for doc in documents:
        lines=doc.split('\n')
        entity_list=[]
        if len(lines)>2:
            first_entity=lines[2].split('\t')
            nest_list=[first_entity]
            max_eid=int(first_entity[2])
            total_entity+=len(lines)-2
            for i in range(3,len(lines)):
                segs=lines[i].split('\t')
                if int(segs[1])> max_eid:
                    if len(nest_list)==1:
                        entity_list.append(nest_list[0])
                        nest_list=[]
                        nest_list.append(segs)
                        if int(segs[2])>max_eid:
                            max_eid=int(segs[2])
                    else:
                        # print(nest_list)
                        nest_entity+=len(nest_list)-1
                        tem=find_max_entity(nest_list)#find max entity
                        # if len(tem)>1:
                        #     print('max nest >1:',tem)
                        entity_list.extend(tem)
                        nest_list=[]
                        nest_list.append(segs)
                        if int(segs[2])>max_eid:
                            max_eid=int(segs[2])
                        
                else:
                    nest_list.append(segs)
                    if int(segs[2])>max_eid:
                        max_eid=int(segs[2])
            if nest_list!=[]:
                if len(nest_list)==1:
                    entity_list.append(nest_list[0])

                else:
                    tem=find_max_entity(nest_list)#find max entity
                    # if len(tem)>1:
                    #     print('max nest >1:',tem)
                    entity_list.extend(tem)
        fout.write(lines[0]+'\n'+lines[1]+'\n')
        for ele in entity_list:
            fout.write('\t'.join(ele)+'\n')
        fout.write('\n')
    # print(total_entity,over_entity, nest_entity)
    return fout
def find_max_entity(nest_list): #longest entity
    max_len=0
    final_tem=[]
    max_index=0
    for i in range(0, len(nest_list)):
        cur_len=int(nest_list[i][2])-int(nest_list[i][1])
        if cur_len>max_len:
            max_len=cur_len
            max_index=i

    final_tem.append(nest_list[max_index])
    return final_tem    
                
# change ori pubtator format to labeled text , entity begin with " ssss", end with 'eeee '
def pubtator_to_labeltext(infile):
    
    fin=io.StringIO(infile.getvalue())
    all_context=fin.read().strip().split('\n\n')
    fin.close()
    fout=io.StringIO()
    label_dic={}
    
    for doc in all_context:
        lines=doc.split('\n')
        ori_text=lines[0].split('|t|')[1]+' '+lines[1].split('|a|')[1]
        pmid=lines[0].split('|t|')[0]
        s_index=0
        e_index=0
        new_text=''
        for i in range(2,len(lines)):
            segs=lines[i].split('\t')
            label_dic[segs[4].lower()]=segs[4]
            if len(segs)==6:
                e_index=int(segs[1])
                new_text+=ori_text[s_index:e_index]+' ssss'+segs[4].lower()+' '+ori_text[int(segs[1]):int(segs[2])]+' eeee'+segs[4].lower()+' '  
                s_index=int(segs[2])
                # if ori_text[int(segs[1]):int(segs[2])]!=segs[3]:
                #     print('error(ori,label):',ori_text[int(segs[1]):int(segs[2])],segs[3])

        new_text+=ori_text[s_index:] 
        fout.write(pmid+'\t'+' '.join(new_text.strip().split())+'\n')
    return fout, label_dic
  

def pre_token(sentence):
    sentence=re.sub("([\=\/\(\)\<\>\+\-\_])"," \\1 ",sentence)
    sentence=re.sub("[ ]+"," ",sentence);
    return sentence

# labeltext to conll format (BIO), a token (including features) per line. sentences are split by '\n', or docs are split by '\n'
def labeltext_to_conll_fasttoken(infile,label_dic):
    
    fin=io.StringIO(infile.getvalue())
    all_context=fin.read().strip().split('\n')
    fin.close()
    fout=io.StringIO()
    
    # nlp = stanza.Pipeline(lang='en', processors='tokenize',package='craft') #package='craft'
    nlp = stanza.Pipeline(lang='en', processors={'tokenize': 'spacy'},package='None') #package='craft'

    doc_i=0
    for doc in all_context:
        doc_text=doc.split('\t')[1]
        doc_text=pre_token(doc_text)
        doc_stanza = nlp(doc_text)
        doc_i+=1
        #print(doc_i)
        inentity_flag=0
        last_label='O'
        for sent in doc_stanza.sentences:
            temp_sent=[]
            word_num=0
            for word in sent.words:
                word_num+=1
                # print(word.text)
                if word.text.strip()=='':
                    continue
                temp_sent.append(word.text)
                if word.text.startswith('ssss')==True:
                    last_label=word.text
                    inentity_flag=1
                elif word.text.startswith('eeee')==True:
                    last_label=word.text
                    inentity_flag=0                    
                else:
                    if last_label=='O':
                        now_label='O'
                    elif last_label.startswith('ssss')==True:
                        now_label='B-'+label_dic[last_label[4:]]
                        
                    elif last_label.startswith('B-')==True:
                        now_label='I-'+last_label[2:]
                    elif last_label.startswith('I-')==True:
                        now_label='I-'+last_label[2:]    
                    elif last_label.startswith('eeee')==True:
                        now_label='O'
                        
                    fout.write(word.text+'\t'+now_label+'\n')
                    last_label=now_label
            if inentity_flag==1: # if entity is split by sentence, will connate the sentence
                # print('sentence error!!!')
                # print(word.text,word_num)
                # print(temp_sent)
                pass
            else:
                fout.write('\n')
    return fout
        
def pubtator_to_conll(infile):
    
    #1.entity sort 
    input_sort=pubtator_entitysort(infile)
    #print(input_sort.getvalue())
    
    #2. no overlap, if overlap get longest entity
    input_nonest=filter_overlap(input_sort)
    # print('......sort.....\n',input_sort.getvalue())
    
    #3. pubtator to label text
    input_labtext,label_dic=pubtator_to_labeltext(input_nonest)
    # print('......label.....\n',input_labtext.getvalue())
    #print(label_dic)
    
    #4. label text to conll
    output = labeltext_to_conll_fasttoken(input_labtext,label_dic)
    # print('......output.....\n',output.getvalue())
    # fout=open(outfile,'w',encoding='utf-8')
    # fout.write(input_nonest.getvalue())
    # fout.close()
    return output

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='python BIO_format.py -i inputfile -o outputfile')
    parser.add_argument('--inputfile', '-i', help="inputfile in PubTator format",default='corpus/NLM-Gene2.2nd/merge.TrainingList.6.GN.PubTator')
    parser.add_argument('--outputfile', '-o', help="outputfile in BIO-conll format",default='corpus/NLM-Gene2.2nd/merge.TrainingList.6.GN.conll')
    args = parser.parse_args()
    
    infile=args.inputfile
    outfile=args.outputfile
    
    output=pubtator_to_conll(infile)
    fout=open(outfile,'w',encoding='utf-8')
    fout.write(output.getvalue())
    fout.close()
    output.close()
    


    
   
    