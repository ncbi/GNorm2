# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 16:34:12 2020

@author: luol2
"""
import numpy as np
import io
import sys
#read ner text (word\tlabel), generate the list[[[w1,label],[w2,label]]]
def ml_intext(file):
    fin=open(file,'r',encoding='utf-8')
    alltexts=fin.read().strip().split('\n\n')
    fin.close()
    data_list=[]
    label_list=[]

    for sents in alltexts:
        lines=sents.split('\n')
        temp_sentece=[]
        for i in range(0,len(lines)):
            seg=lines[i].split('\t')
            temp_sentece.append(seg[:])
            label_list.append(seg[-1])
        
        data_list.append(temp_sentece)
    #print(data_list)
    #print(label_list)
    return data_list,label_list

def ml_intext_fn(alltexts):
    # fin=io.StringIO(ml_input)
    # alltexts=fin.read().strip().split('\n\n')
    # fin.close()
    data_list=[]
    label_list=[]

    for sents in alltexts:
        lines=sents.split('\n')
        temp_sentece=[]
        for i in range(0,len(lines)):
            seg=lines[i].split('\t')
            temp_sentece.append(seg[:])
            label_list.append(seg[-1])
        
        data_list.append(temp_sentece)
    #print(data_list)
    #print(label_list)
    return data_list,label_list

# model predict result to conll evalute format  [token answer predict]        
def out_BIO(file,raw_pre,raw_input,label_set):
    fout=open(file,'w',encoding='utf-8')
    for i in range(len(raw_input)):
        
        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                label_id = raw_pre[i][j]
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\n')
        fout.write('\n')
    fout.close()
    
def out_BIO_softmax(file,raw_pre,raw_input,label_set):
    fout=open(file,'w',encoding='utf-8')
    #print(raw_pre[0:2])
    for i in range(len(raw_input)):
        
        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                label_id = np.argmax(raw_pre[i][j])
                #print(label_id)
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\n')
        fout.write('\n')
    fout.close()
    
def out_BIO_fn(raw_pre,raw_input,label_set):
    fout=io.StringIO()
    for i in range(len(raw_input)):
        
        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                label_id = raw_pre[i][j] 
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\n')
        fout.write('\n')
    return fout.getvalue()  
def out_BIO_BERT_softmax(file,raw_pre,raw_input,label_set):
    fout=open(file,'w',encoding='utf-8')
    for i in range(len(raw_input)):
        
        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                # label_id = raw_pre[i][j]
                label_id = np.argmax(raw_pre[i][j])
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\n')
        fout.write('\n')
    fout.close() 
def out_BIO_BERT(file,raw_pre,raw_input,label_set):
    fout=open(file,'w',encoding='utf-8')
    for i in range(len(raw_input)):
        
        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                label_id = raw_pre[i][j] 
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\n')
        fout.write('\n')
    fout.close() 
def out_BIO_BERT_fn(raw_pre,raw_input,label_set):
    fout=io.StringIO()
    for i in range(len(raw_input)):
        
        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                label_id = raw_pre[i][j] 
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\n')
        fout.write('\n')
    return fout.getvalue()                  
def out_BIO_BERT_softmax_fn(raw_pre,raw_input,label_set):
    fout=io.StringIO()
    for i in range(len(raw_input)):

        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                #label_id = raw_pre[i][j]
                label_id = np.argmax(raw_pre[i][j])
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\n')
        fout.write('\n')
    return fout.getvalue()
def out_BIO_BERT_softmax_score_fn(raw_pre,raw_input,label_set):
    fout=io.StringIO()
    for i in range(len(raw_input)):

        for j in range(len(raw_input[i])):
            if j<len(raw_pre[i]):
                #label_id = raw_pre[i][j]
                label_id = np.argmax(raw_pre[i][j])
                label_score = round(raw_pre[i][j][label_id],4)
                label_tag = label_set[str(label_id)]
            else:
                label_tag='O'
                label_score = 0.0
            fout.write(raw_input[i][j][0]+'\t'+raw_input[i][j][-1]+'\t'+label_tag+'\t'+str(label_score)+'\n')
        fout.write('\n')
    return fout.getvalue()
#generate char vocab
def char_vocab(infile,outfile_char):
    fin=open(infile,'r',encoding='utf-8')
    #fout=open(outfile,'w',encoding='utf-8')
    fout_char=open(outfile_char,'w',encoding='utf-8')
    char_vocab=['oov_char']
    max_len=0
    for line in fin:
        if line.strip()!='':
            seg=line.split('\t')
            word_len=len(seg[0])
            #if word_len<1000:
            #    fout.write(line)
            if word_len>max_len:
                max_len=word_len
                print(seg[0])
            for i in range(word_len):
                if seg[0][i] not in char_vocab:
                    char_vocab.append(seg[0][i])
        #else:
        #    fout.write(line)
    fin.close()
    #fout.close()
    for ele in char_vocab:
        fout_char.write(ele+'\n')
    fout_char.close()
    print('max_len:',max_len)


if __name__=='__main__':
    # infile='//panfs/pan1/bionlp/lulab/luoling/HPO_project/AutoPhe/data/pubmed_unlabel/mutation_disease_1990.ner_BIO'
    # #outfile='//panfs/pan1/bionlp/lulab/luoling/HPO_project/AutoPhe/data/pubmed_unlabel/mutation_disease_1990.ner_BIO_new'
    # outfile_char='//panfs/pan1/bionlp/lulab/luoling/HPO_project/AutoPhe/src/nn_model/vocab/char_vocab'
    # #processing_text(file)
    # char_vocab(infile,outfile_char)
    a=[1,2,3]
    print(a[:-1])
