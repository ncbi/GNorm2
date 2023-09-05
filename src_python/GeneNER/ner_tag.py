# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 11:01:23 2022

@author: luol2
"""



import io
import re
from src_python.GeneNER.processing_data_ner import ml_intext_fn,out_BIO_BERT_softmax_fn
from src_python.GeneNER.restore_index_ner import NN_restore_index_fn
import tensorflow as tf
gpu = tf.config.list_physical_devices('GPU')
print("Num GPUs Available: ", len(gpu))
if len(gpu) > 0:
    tf.config.experimental.set_memory_growth(gpu[0], True)

def pre_token(sentence):
    sentence=re.sub("([\W\-\_])"," \\1 ",sentence)
    sentence=re.sub("[ ]+"," ",sentence);
    return sentence

def ssplit_token_pos_lemma(in_text,text_level,nlp_token, max_len=400):
    #print('max_len:',max_len)
    fout=io.StringIO()

    in_text=in_text.strip()
    in_text=pre_token(in_text)
    doc_stanza = nlp_token(in_text)
    strlen=0
    for sent in doc_stanza.sentences:
        for word in sent.words:
            strlen+=1
            if word.text.strip()=='':
                pass
                #print('!!!!blank token text!!!')
            else:
                fout.write(word.text+'\tO\n')
            if strlen>=max_len:
                #print('long sentence:',strlen)
                fout.write('\n')
                strlen=0
        if text_level=='SENT':
            fout.write('\n')
            strlen=0
    if text_level=='DOC':
        fout.write('\n')
           
    return fout.getvalue()

def ml_tagging(ml_input,nn_model):

    test_list = ml_intext_fn(ml_input)
    test_x,test_y, test_bert_text_label=nn_model.rep.load_data_hugface(test_list,word_max_len=nn_model.maxlen,label_type='softmax')
    test_pre = nn_model.model.predict(test_x,batch_size=64)
    test_decode_temp=out_BIO_BERT_softmax_fn(test_pre,test_bert_text_label,nn_model.rep.index_2_label)
    
    return test_decode_temp
# only machine learning-based method
def ML_Tag(text,ml_model,nlp_token,text_level='SENT'):

#    startTime=time.time()
    ssplit_token=ssplit_token_pos_lemma(text, text_level, nlp_token, max_len=ml_model.maxlen)
    #print(ssplit_token) 
#    print('ssplit token:',time.time()-startTime)
    
#    startTime=time.time()
    ml_tsv=ml_tagging(ssplit_token,ml_model)
    #print(ml_tsv)
#    print('ml ner:',time.time()-startTime)
   
    final_result= NN_restore_index_fn(text,ml_tsv)

    # print('final ner:',time.time()-startTime)
    
    return final_result





    

