# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 10:40:08 2021

@author: luol2
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 17:19:02 2020

@author: luol2
"""

import io
import sys

# from BIO format to entity,list line is sentence, follwing the entity(start, end, text, entity, type)
def NN_BIO_tag_entity(pre_BIO):
    sentences=pre_BIO.strip().split('\n\n')

    pre_result=[]
    #print(sentences)
    for sent in sentences:
        tokens=sent.split('\n')
        pre_entity=[]
        pre_start,pre_end=0,0
        sent_text=''
        for i in range(0,len(tokens)):
            segs=tokens[i].split('\t')
            sent_text+=segs[0]+' '
            if len(segs)<3:
                continue
            #print(tokens)
            # generate prediction entity            
            if segs[2].startswith('B-')>0:
                pre_start=i
                pre_type=segs[2][2:]
                if i+1>=len(tokens): # the last word
                    pre_end=i
                    pre_entity.append([pre_start,pre_end,pre_type])
                else: # non last word
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start,pre_end,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2].startswith('I-')>0:
                if i==0 and i+1<len(tokens): # the first word and not only a word
                    pre_start=i
                    pre_type=segs[2][2:]
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start,pre_end,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
                elif i==0 and i+1==len(tokens):# only one word:
                    pre_start=i
                    pre_type=segs[2][2:]
                    pre_end=i
                    pre_entity.append([pre_start,pre_end,pre_type])
                elif i+1>=len(tokens): # the last word
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]                
                    pre_end=i
                    pre_entity.append([pre_start,pre_end,pre_type])
                elif i+1< len(tokens): # non last word
                    next_seg=tokens[i+1].split('\t')
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]  
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start,pre_end,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2]=='O':
                pass        
        pre_result.append([sent_text.rstrip(),pre_entity])


        # print(pre_entity)
    return pre_result

def NN_restore_index_fn(ori_text,file_pre):

    input_result=NN_BIO_tag_entity(file_pre)
    #print(input_result)
    
    
    new_sentence=''
    restore_result=[]
    
    sentence_ori=ori_text.lower()

    for sent_ele in input_result:

        #print(pre_lines)
#        print(sentence_ori)
        if len(sent_ele[1])>0:
            #print(pre_lines)
            sentence_pre=sent_ele[0].lower()
            sentence_pre=sentence_pre.split()
            
            pre_result=sent_ele[1]

            
            restore_sid=0
            restore_eid=0
            each_word_id=[]
            
            for i in range(0,len(sentence_pre)):

                temp_id=sentence_ori.find(sentence_pre[i])
                if temp_id<0:
                        #print('ori:',sentence_ori)
                        print('resotr index error:',sentence_pre[i])
                new_sentence+=sentence_ori[0:temp_id]
                
                restore_sid=len(new_sentence)
                restore_eid=len(new_sentence)+len(sentence_pre[i])
                each_word_id.append([str(restore_sid),str(restore_eid)])
                new_sentence+=sentence_ori[temp_id:temp_id+len(sentence_pre[i])]
                sentence_ori=sentence_ori[temp_id+len(sentence_pre[i]):]
#            print('each_word:',each_word_id)    
            for pre_ele in pre_result:
                temp_pre_result=[each_word_id[int(pre_ele[0])][0],each_word_id[int(pre_ele[1])][1],pre_ele[2]]
                if temp_pre_result not in restore_result:
                    restore_result.append(temp_pre_result)
        else:
            sentence_pre=sent_ele[0].lower()
            sentence_pre=sentence_pre.split()
           
            for i in range(0,len(sentence_pre)):

                temp_id=sentence_ori.find(sentence_pre[i])
                if temp_id<0:
                    print('resotr index error:',sentence_pre[i])
                new_sentence+=sentence_ori[0:temp_id]
                new_sentence+=sentence_ori[temp_id:temp_id+len(sentence_pre[i])]
                sentence_ori=sentence_ori[temp_id+len(sentence_pre[i]):]
    #print('resotre:',restore_result)
    return restore_result

def BERT_BIO_tag_entity(pre_BIO):
    sentences=pre_BIO.strip().split('\n\n')

    pre_result=[]
    for sent in sentences:
        tokens=sent.split('\n')
        pre_entity=[]
        pre_start,pre_end=0,0
        sent_text=''
        for i in range(1,len(tokens)-1):
            segs=tokens[i].split('\t')
            sent_text+=segs[0]+' '
            # generate prediction entity            
            if segs[2].startswith('B-')>0:
                pre_start=i
                pre_type=segs[2][2:]
                if i+1>=len(tokens): # the last word
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                else: # non last word
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2].startswith('I-')>0:
                if i==0 and i+1<len(tokens): # the first word and not only a word
                    pre_start=i
                    pre_type=segs[2][2:]
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
                elif i==0 and i+1==len(tokens):# only one word:
                    pre_start=i
                    pre_type=segs[2][2:]
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                elif i+1>=len(tokens): # the last word
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]                
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                elif i+1< len(tokens): # non last word
                    next_seg=tokens[i+1].split('\t')
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]  
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2]=='O':
                pass        
        pre_result.append([sent_text.rstrip(),pre_entity])


    #print(pre_result)
    return pre_result

def BERT_BIO_tag_entity_revised(pre_BIO):
    print('revised version')
    sentences=pre_BIO.strip().split('\n\n')

    pre_result=[]
    for sent in sentences:
        tokens=sent.split('\n')
        pre_entity=[]
        pre_start,pre_end=0,0
        sent_text=''
        for i in range(1,len(tokens)-1):
            segs=tokens[i].split('\t')
            sent_text+=segs[0]+' '
            # generate prediction entity            
            if segs[2].startswith('B-')>0:
                pre_start=i
                pre_type=segs[2][2:]
                if i+1>=len(tokens)-1: # the last word
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                else: # non last word
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2].startswith('I-')>0:
                if i==1 and i+1<len(tokens)-1: # the first word and not only a word
                    pre_start=i
                    pre_type=segs[2][2:]
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
                elif i==1 and i+1==len(tokens)-1:# only one word:
                    pre_start=i
                    pre_type=segs[2][2:]
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                elif i+1>=len(tokens)-1: # the last word
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]                
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                elif i+1< len(tokens)-1: # non last word
                    next_seg=tokens[i+1].split('\t')
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]  
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2]=='O':
                pass        
        pre_result.append([sent_text.rstrip(),pre_entity])


    #print(pre_result)
    return pre_result

# only predict on the first token of the ori word 
def BERT_BIO_tag_entity_word(pre_BIO):
    sentences=pre_BIO.strip().split('\n\n')

    pre_result=[]
    for sent in sentences:
        tokens=sent.split('\n')
        pre_entity=[]
        pre_start,pre_end=0,0
        sent_text=''
        i=1
        while i< len(tokens)-1:
        # for i in range(1,len(tokens)-1):
            segs=tokens[i].split('\t')
            sent_text+=segs[0]+' '
            # generate prediction entity            
            if segs[2].startswith('B-')>0:
                pre_start=i
                pre_type=segs[2][2:]
                if i+1>=len(tokens)-1: # the last word
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                else: # non last word
                    #pass a word
                    sub_segs=tokens[i+1].split('\t')
                    while(sub_segs[0].find('##')==0):
                        i+=1
                        sent_text+=sub_segs[0]+' '
                        sub_segs=tokens[i+1].split('\t')
                        
                        
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2].startswith('I-')>0:
                if i==1 and i+1<len(tokens)-1: # the first word and not only a word
                    pre_start=i
                    pre_type=segs[2][2:]
                    #pass a word
                    sub_segs=tokens[i+1].split('\t')
                    while(sub_segs[0].find('##')==0):
                        i+=1
                        sent_text+=sub_segs[0]+' '
                        sub_segs=tokens[i+1].split('\t')
                    
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
                elif i==1 and i+1==len(tokens)-1:# only one word:
                    pre_start=i
                    pre_type=segs[2][2:]
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                elif i+1>=len(tokens)-1: # the last word
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]                
                    pre_end=i
                    pre_entity.append([pre_start-1,pre_end-1,pre_type])
                elif i+1< len(tokens)-1: # non last word
                    
                    last_seg=tokens[i-1].split('\t')
                    if last_seg[2]=='O':
                        pre_start=i
                        pre_type=segs[2][2:]
                        #pass a word
                        sub_segs=tokens[i+1].split('\t')
                        while(sub_segs[0].find('##')==0):
                            i+=1
                            sent_text+=sub_segs[0]+' '
                            sub_segs=tokens[i+1].split('\t')
                    next_seg=tokens[i+1].split('\t')
                    if next_seg[2].startswith('B-')>0 or next_seg[2]=='O':
                        pre_end=i
                        pre_entity.append([pre_start-1,pre_end-1,pre_type])
                    elif next_seg[2].startswith('I-')>0:
                        pass
            elif segs[2]=='O':
                pass
            i+=1
        pre_result.append([sent_text.rstrip(),pre_entity])


    #print(pre_result)
    return pre_result


def BERT_restore_index_fn(ori_text,file_pre):

    # input_result=BERT_BIO_tag_entity_revised(file_pre)
    input_result=BERT_BIO_tag_entity_word(file_pre)
    #print(input_result)
    
    
    new_sentence=''
    restore_result=[]
    
    sentence_ori=ori_text.lower()

    for sent_ele in input_result:

        #print(pre_lines)
#        print(sentence_ori)
        if len(sent_ele[1])>0:
            #print(pre_lines)
            sentence_pre=sent_ele[0].lower()
            sentence_pre=sentence_pre.split()
            
            pre_result=sent_ele[1]

            
            restore_sid=0
            restore_eid=0
            each_word_id=[]
            
            
            for i in range(0,len(sentence_pre)):
                if sentence_pre[i][0:2]=="##":
                    sentence_pre[i]=sentence_pre[i][2:]
                temp_id=sentence_ori.find(sentence_pre[i])
                if temp_id<0:
                        #print('ori:',sentence_ori)
                        print('resotr index error:',sentence_pre[i])
                new_sentence+=sentence_ori[0:temp_id]
                
                restore_sid=len(new_sentence)
                restore_eid=len(new_sentence)+len(sentence_pre[i])
                each_word_id.append([str(restore_sid),str(restore_eid)])
                new_sentence+=sentence_ori[temp_id:temp_id+len(sentence_pre[i])]
                sentence_ori=sentence_ori[temp_id+len(sentence_pre[i]):]
#            print('each_word:',each_word_id)    
            for pre_ele in pre_result:
                temp_pre_result=[each_word_id[int(pre_ele[0])][0],each_word_id[int(pre_ele[1])][1],pre_ele[2]]
                if temp_pre_result not in restore_result:
                    restore_result.append(temp_pre_result)
        else:
            sentence_pre=sent_ele[0].lower()
            sentence_pre=sentence_pre.split()
           
            for i in range(0,len(sentence_pre)):
                if sentence_pre[i][0:2]=="##":
                    sentence_pre[i]=sentence_pre[i][2:]
                temp_id=sentence_ori.find(sentence_pre[i])
                if temp_id<0:
                    print('resotr index error:',sentence_pre[i])
                new_sentence+=sentence_ori[0:temp_id]
                new_sentence+=sentence_ori[temp_id:temp_id+len(sentence_pre[i])]
                sentence_ori=sentence_ori[temp_id+len(sentence_pre[i]):]
    #print('resotre:',restore_result)
    return restore_result
if __name__=='__main__':
    path='//panfs/pan1/bionlp/lulab/luoling/OpenBioIE_project/models/'
    fin=open(path+'devout_test.txt','r',encoding='utf-8')
    file_pre=fin.read()
    ori_text="D90A-SOD1 mediated amyotrophic lateral sclerosis: a single founder for all cases with evidence for a Cis-acting disease modifier in the recessive haplotype. More than 100 different heterozygous mutations in copper/zinc superoxide dismutase (SOD1) have been found in patients with amyotrophic lateral sclerosis (ALS), a fatal neurodegenerative disease. Uniquely, D90A-SOD1 has been identified in recessive, dominant and apparently sporadic pedigrees. The phenotype of homozygotes is stereotyped with an extended survival, whereas that of affected heterozygotes varies. The frequency of D90A-SOD1 is 50 times higher in Scandinavia (2.5%) than elsewhere, though ALS prevalence is not raised there. Our earlier study indicated separate founders for recessive and dominant/sporadic ALS and we proposed a disease-modifying factor linked to the recessive mutation. Here we have doubled our sample set and employed novel markers to characterise the mutation's origin and localise any modifying factor. Linkage disequilibrium analysis indicates that D90A homozygotes and heterozygotes share a rare haplotype and are all descended from a single ancient founder (alpha 0.974) c.895 generations ago. Homozygotes arose subsequently only c.63 generations ago (alpha 0.878). Recombination has reduced the region shared by recessive kindreds to 97-265 kb around SOD1, excluding all neighbouring genes. We propose that a cis-acting regulatory polymorphism has arisen close to D90A-SOD1 in the recessive founder, which decreases ALS susceptibility in heterozygotes and slows disease progression."
    NN_restore_index_fn(ori_text,file_pre)
