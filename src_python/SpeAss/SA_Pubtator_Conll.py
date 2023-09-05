# -*- coding: utf-8 -*-

import sys
import io
import argparse
import stanza
# nlp = stanza.Pipeline(lang='en', processors='tokenize,mwt,pos,lemma',package='craft') #package='craft'
nlp = stanza.Pipeline(lang='en', processors={'tokenize': 'spacy'},package='None') #package='craft'
REL_ENT={'arg1':'Species',
         'arg2':'Gene'}

ENTITY_TAG={'arg1':['arg1s','arg1e'],
            'arg2':['arg2s','arg2e'],
            'gene':['gene1s','gene1e'],
            'species':['species1s','species1e']
            }

# ssplit token and revise index
def ssplit_token(infile):
    fin=open(infile,'r',encoding='utf-8')
    fout=io.StringIO()
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    for doc_text in all_in:
        lines=doc_text.split('\n')
        ori_text=lines[0].split('|t|')[1]+' '+lines[1].split('|a|')[1]
        pmid=lines[0].split('|t|')[0]
        # print(pmid)
        entity_all=[]   #[[seg0,seg1,...,],[]]
        for i in range(2,len(lines)):
            seg=lines[i].split('\t')
            entity_all.append(seg)

        #ssplit token
        doc_stanza = nlp(ori_text)
        token_text=''
        for sent in doc_stanza.sentences:
            for word in sent.words:
                if word.text==' ':
                    pass
                    # print('token is blank!')
                else:
                    token_text+=word.text+' '
            #token_text=token_text+'    '  #sentence split by four blank
        
        #ori_index map token_index
        index_map=[-1]*len(ori_text)
        j=0
        space_list=[' ',chr(160),chr(8201),chr(8194),chr(8197),chr(8202)] #空格有好几种，第一个是常用32,第二个shi 160,8201,8194,8197
        for i in range(0,len(ori_text)):
            if ori_text[i] in space_list:
                pass
            elif ori_text[i]==token_text[j]:
                #if i>0 and i<285:
                #    print('=i,j:',i,j,ori_text[i-1:i+1],token_text[j-1:j+1])
                index_map[i]=j
                j+=1
            else:
                #if i==283:
                #    print('!i,j:',i,j,ori_text[i-1:i+1],token_text[j-1:j+1])
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
        # print(index_map)
        # token_text=token_text.replace('     ','<EOS>')
        # print(token_text)
        fout.write(token_text+'\n')
        for ele in entity_all:
            if index_map[int(ele[1])]==-1:
                new_ents=index_map[int(ele[1])+1]
            else:
                new_ents=index_map[int(ele[1])]
            if index_map[int(ele[2])-1]==-1: 
                new_ente=index_map[int(ele[2])-1-1]+1
            else:
                new_ente=index_map[int(ele[2])-1]+1
            new_ent=token_text[new_ents:new_ente]
            if ele[4]=='Species' or ele[4]=='Gene':
                fout.write(ele[0]+'\t'+str(new_ents)+'\t'+str(new_ente)+'\t'+new_ent+'\t'+ele[4]+'\t'+ele[5]+'\n')
            else:
                # print(ele[4])
                fout.write(ele[0]+'\t'+str(new_ents)+'\t'+str(new_ente)+'\t'+new_ent+'\t'+'Gene'+'\t'+ele[5]+'\n')
        fout.write('\n')
    return fout.getvalue()


def corpus_noNest(token_input):
    
    fin=io.StringIO(token_input)
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
            doc_result={}
            for i in range(1,len(lines)):
                segs=lines[i].split('\t')
                doc_result[lines[i]]=[int(segs[1]),int(segs[2])]
            doc_result=sorted(doc_result.items(), key=lambda kv:(kv[1]), reverse=False)
            doc_result_sort=[]
            for ele in doc_result:
                doc_result_sort.append(ele[0])
            
            first_entity=doc_result_sort[0].split('\t')
            nest_list=[first_entity]
            max_eid=int(first_entity[2])
            total_entity+=len(lines)-2
            for i in range(1,len(doc_result_sort)):
                segs=doc_result_sort[i].split('\t')
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
                        tem=find_max_entity(nest_list,context)#find max entity
                        # if len(tem)>1:
                            # print('max nest >1:',tem)
                        entity_list.extend(tem)
                        nest_list=[]
                        nest_list.append(segs)
                        if int(segs[2])>max_eid:
                            max_eid=int(segs[2])
                        
                else:
                    nest_list.append(segs)
                    over_entity+=1
                    if int(segs[2])>max_eid:
                        max_eid=int(segs[2])
            if nest_list!=[]:
                if len(nest_list)==1:
                    entity_list.append(nest_list[0])

                else:
                    tem=find_max_entity(nest_list,context)#find max entity
                    # if len(tem)>1:
                    #     print('max nest >1:',tem)
                    entity_list.extend(tem)
        fout.write(context+'\n')
        for ele in entity_list:
            if ele[4]=='Gene':
                temp_gene={}
                gene_ids=ele[5].split(',')
                for gene_id in gene_ids:
                    temp_id=gene_id[gene_id.find('Species:'):-1]
                    spe_id=temp_id[len('Species:'):]
                    temp_gene[temp_id]=int(spe_id)
                temp_gene_sort=sorted(temp_gene.items(), key=lambda kv:(kv[1]), reverse=False)
                final_gene_id=''
                for temp_ele in temp_gene_sort:
                    final_gene_id+=temp_ele[0]+','
                fout.write('\t'.join(ele[:-1])+'\t'+final_gene_id[:-1]+'\n')
            else:
                fout.write('\t'.join(ele)+'\n')
        fout.write('\n')
    # print(total_entity,over_entity, nest_entity)
    return fout.getvalue()

def find_max_entity(nest_list,text):
    max_len=0
    final_tem=[]
    max_index=0
    for i in range(0, len(nest_list)):
        if nest_list[i][4] =='Species':
            final_tem.append(nest_list[i])
        else:
            cur_len=int(nest_list[i][2])-int(nest_list[i][1])
            if cur_len>max_len:
                max_len=cur_len
                max_index=i
    final_tem.append(nest_list[max_index])
    return final_tem


def generate_seq_input(nonest_input,outfile):
    
    fin=io.StringIO(nonest_input)
    fout=open(outfile,'w',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()

    final_input=[]

    for doc in all_in:
        lines=doc.split('\n')
        token_text=lines[0]
        pmid=lines[1].split('\t')[0]
        # print(pmid)
        #read entity and relation
        entity_arg1={} #only entity offset
        entity_arg2={} #only entity offset
        entity_all=[] #all entity infor

        for i in range(1,len(lines)):
            seg=lines[i].split('\t')
            if seg[4]==REL_ENT['arg1']:
                if seg[-1] in entity_arg1.keys():
                    entity_arg1[seg[-1]].append([seg[1],seg[2]])
                else:
                    entity_arg1[seg[-1]]=[[seg[1],seg[2]]]
            elif seg[4]==REL_ENT['arg2']:
                temp_spes=seg[-1].split(',')
                for ele in temp_spes:
                    gene_spe_id=ele
                    if gene_spe_id in entity_arg2.keys():
                        entity_arg2[gene_spe_id].append([seg[1],seg[2]])
                    else:
                        entity_arg2[gene_spe_id]=[[seg[1],seg[2]]]
                
            entity_all.append(seg)
        # print('\narg1:',entity_arg1)
        # print('\narg2:',entity_arg2)
        # print('\nall entity:',entity_all)
        # for all arg1 to produce inst
        for cur_ele in entity_arg1.keys():

            #1. ner label text
            #check cur_ele in relation?
            # print(relation_all.keys())
            if cur_ele in entity_arg2.keys(): #pos instance
                rel_ent2=entity_arg2[cur_ele]
                ner_text=''
                text_sid=0
                #print('nonest:',entity_nonest)
                for ele_nonest in entity_all:
                    ent_id=[ele_nonest[1],ele_nonest[2]]
                    ent_sid=int(ele_nonest[1])
                    ent_eid=int(ele_nonest[2])
                    # print('sid,eid:',ent_sid,ent_eid)
                    ent_text=ele_nonest[3]
                    ent_type=ele_nonest[4]
                    if ent_sid>=text_sid:
                        if ent_id in entity_arg1[cur_ele]:    
                            ner_text+=token_text[text_sid:ent_sid]+' '+ENTITY_TAG['arg1'][0]+' '+ent_text+ ' '+ENTITY_TAG['arg1'][1]+' '
                        else:
                            if ent_id in rel_ent2: #arg2 entity
                                if ent_type!=REL_ENT['arg2']:
                                    pass
                                    # print('arg2 is error! not ',REL_ENT['arg2'], ele_nonest)
                                ner_text+=token_text[text_sid:ent_sid]+' '+ENTITY_TAG['arg2'][0]+' '+ent_text+ ' '+ENTITY_TAG['arg2'][1]+' '
                            else:
                                ner_text+=token_text[text_sid:ent_sid]+' '+ENTITY_TAG[ent_type.lower()][0]+' '+ent_text+ ' '+ENTITY_TAG[ent_type.lower()][1]+' '
                        text_sid=ent_eid  
                    else:
                        pass
                        # print('ner entity error!!!',ele_nonest,text_sid)                                    
                ner_text+=token_text[text_sid:]
                sen_tokens=ner_text.split()
                # print('\nner_text:',ner_text)
                
                #3 produce pos input
 
                temp_input=[]
                token_id=0
                while token_id <len(sen_tokens):
                    if sen_tokens[token_id].find(ENTITY_TAG['arg1'][0])>=0:
                        temp_input.append(ENTITY_TAG['arg1'][0]+'\tO')
                        token_id+=1
                        while(sen_tokens[token_id]!=ENTITY_TAG['arg1'][1]):
                            temp_input.append(sen_tokens[token_id]+'\tO')
                            token_id+=1
                        temp_input.append(ENTITY_TAG['arg1'][1]+'\tO')
                    elif sen_tokens[token_id].find(ENTITY_TAG['arg2'][0])>=0:
                        temp_input.append(ENTITY_TAG[REL_ENT['arg2'].lower()][0]+'\tARG2')
                        token_id+=1
                        while(sen_tokens[token_id]!=ENTITY_TAG['arg2'][1]):
                            temp_input.append(sen_tokens[token_id]+'\tARG2')
                            token_id+=1
                        temp_input.append(ENTITY_TAG[REL_ENT['arg2'].lower()][1]+'\tARG2')
                    elif sen_tokens[token_id].find(ENTITY_TAG['gene'][0])>=0:
                        temp_input.append(ENTITY_TAG['gene'][0]+'\tO')
                        token_id+=1
                        while(sen_tokens[token_id]!=ENTITY_TAG['gene'][1]):
                            temp_input.append(sen_tokens[token_id]+'\tO')
                            token_id+=1
                        temp_input.append(ENTITY_TAG['gene'][1]+'\tO')
                    elif sen_tokens[token_id].find(ENTITY_TAG['species'][0])>=0:
                        temp_input.append(ENTITY_TAG['species'][0]+'\tO')
                        token_id+=1
                        while(sen_tokens[token_id]!=ENTITY_TAG['species'][1]):
                            temp_input.append(sen_tokens[token_id]+'\tO')
                            token_id+=1
                        temp_input.append(ENTITY_TAG['species'][1]+'\tO')
                    else:
                        if sen_tokens[token_id]=='':
                            # print('token is none!error!')
                            pass
                        else:
                            temp_input.append(sen_tokens[token_id]+'\tO')
                    token_id+=1
                    
                final_input.append('\n'.join(temp_input))
                
            else:     #neg instance 
                ner_text=''
                text_sid=0
                #print('nonest:',entity_nonest)
                for ele_nonest in entity_all:
                    ent_id=[ele_nonest[1],ele_nonest[2]]
                    ent_sid=int(ele_nonest[1])
                    ent_eid=int(ele_nonest[2])
                    # print('sid,eid:',ent_sid,ent_eid)
                    ent_text=ele_nonest[3]
                    ent_type=ele_nonest[4]
                    if ent_sid>=text_sid:
                        if ent_id in entity_arg1[cur_ele]:    
                            ner_text+=token_text[text_sid:ent_sid]+' '+ENTITY_TAG['arg1'][0]+' '+ent_text+ ' '+ENTITY_TAG['arg1'][1]+' '
                        else:
                            ner_text+=token_text[text_sid:ent_sid]+' '+ENTITY_TAG[ent_type.lower()][0]+' '+ent_text+ ' '+ENTITY_TAG[ent_type.lower()][1]+' '
                        text_sid=ent_eid  
                    else:
                        pass
                        # print('ner entity error!!!')                                    
                ner_text+=token_text[text_sid:]
                sen_tokens=ner_text.split()
                # print('\nner_text:',ner_text)
                # print('ner_Text')
                #3 produce NEG input
                   
                temp_input=[]
                token_id=0
                while token_id <len(sen_tokens):
                    if sen_tokens[token_id].find(ENTITY_TAG['arg1'][0])>=0:
                        temp_input.append(ENTITY_TAG['arg1'][0]+'\tO')
                        token_id+=1
                        while(sen_tokens[token_id]!=ENTITY_TAG['arg1'][1]):
                            temp_input.append(sen_tokens[token_id]+'\tO')
                            token_id+=1
                        temp_input.append(ENTITY_TAG['arg1'][1]+'\tO')
                    elif sen_tokens[token_id].find(ENTITY_TAG['gene'][0])>=0:
                        temp_input.append(ENTITY_TAG['gene'][0]+'\tO')
                        token_id+=1
                        while(sen_tokens[token_id]!=ENTITY_TAG['gene'][1]):
                            temp_input.append(sen_tokens[token_id]+'\tO')
                            token_id+=1
                        temp_input.append(ENTITY_TAG['gene'][1]+'\tO')
                    elif sen_tokens[token_id].find(ENTITY_TAG['species'][0])>=0:
                        temp_input.append(ENTITY_TAG['species'][0]+'\tO')
                        token_id+=1
                        while(sen_tokens[token_id]!=ENTITY_TAG['species'][1]):
                            temp_input.append(sen_tokens[token_id]+'\tO')
                            token_id+=1
                        temp_input.append(ENTITY_TAG['species'][1]+'\tO')
                    else:
                        if sen_tokens[token_id]=='':
                            print('token is none!error!')
                        else:
                            temp_input.append(sen_tokens[token_id]+'\tO')
                    token_id+=1
                    
                final_input.append('\n'.join(temp_input))
                # print(entity_nonest)
        # sys.exit()
    fout.write('\n\n'.join(final_input))
    fout.write('\n')
    fout.close()   

def check_entity_pos(line,relations):
    
    seg=line.split(' ')
    stack_ent=[] 
    # print(seg)
    entity_num={'arg1':0,'arg2':0, 'gene':0,'chemical':0}

    temp_arg2=[]
    for i in range(0,len(seg)):
        if seg[i].find(ENTITY_TAG['gene'][0])>=0:
            entity_num['gene']+=1
            stack_ent.append(seg[i])
        elif seg[i].find(ENTITY_TAG['chemical'][0])>=0:
            entity_num['chemical']+=1
            stack_ent.append(seg[i])
            # print(stack_ent)
        elif seg[i].find(ENTITY_TAG['arg1'][0])>=0:
            entity_num['arg1']+=1
            stack_ent.append(seg[i]) 
        elif seg[i].find(ENTITY_TAG['arg2'][0])>=0:
            entity_num['arg2']+=1
            temp_arg2.append(seg[i].split('|')[0])
            stack_ent.append(seg[i])
        elif seg[i].find(ENTITY_TAG['arg1'][1])>=0 or seg[i].find(ENTITY_TAG['arg2'][1])>=0 or seg[i].find(ENTITY_TAG['gene'][1])>=0 or seg[i].find(ENTITY_TAG['chemical'][1])>=0:
            stack_ent.pop()
    if stack_ent!=[]:
        # print('entity no match!',stack_ent)
        return(-1,seg,entity_num)
    
    else:
        if entity_num['arg1']!=0:
            for arg2_id in relations.keys():
                if arg2_id not in temp_arg2:
                    # print('\ntemp_arg2:',temp_arg2)
                    # print('\narg2_id:',arg2_id)
                    return(0,seg,entity_num) #some arg2 not in sentence
        if entity_num['arg2']!=0 and entity_num['arg1']==0:
            return(0,seg,entity_num) #only arg2, but no arg1
        return(1,seg,entity_num) 

def check_entity_neg(line):
    
    seg=line.split(' ')
    stack_ent=[] 
    # print(seg)
    entity_num={'arg1':0,'gene':0,'chemical':0}
    for i in range(0,len(seg)):
        if seg[i].find(ENTITY_TAG['gene'][0])>=0:
            entity_num['gene']+=1
            stack_ent.append(seg[i])
        elif seg[i].find(ENTITY_TAG['chemical'][0])>=0:
            entity_num['chemical']+=1
            stack_ent.append(seg[i])
            # print(stack_ent)
        elif seg[i].find(ENTITY_TAG['arg1'][0])>=0:
            entity_num['arg1']+=1
            stack_ent.append(seg[i])        
        elif seg[i].find(ENTITY_TAG['arg1'][1])>=0  or seg[i].find(ENTITY_TAG['gene'][1])>=0 or seg[i].find(ENTITY_TAG['chemical'][1])>=0:
            stack_ent.pop()
    if stack_ent!=[]:
        # print('entity no match!',stack_ent)
        return(-1,seg,entity_num)
    
    else:
        return(1,seg,entity_num) 

def get_one_entity(nest_list,cur_ent,rel_entity2_id):
    max_len=0
    max_entity=[]
    final_entity=[]
    for i in range(0, len(nest_list)):
        if nest_list[i][1]==cur_ent:#current entity
            final_entity=[]
            max_entity=nest_list[i]
            final_entity.append(nest_list[i])
            return(final_entity)
        if nest_list[i][1] in rel_entity2_id: #invole rel
            final_entity.append(nest_list[i])
            continue
        length=int(nest_list[i][4])-int(nest_list[i][3])
        if max_entity==[]: #first entity
            max_len=length
            max_entity=nest_list[i]
        else:
            if length>max_len:
                if max_entity[2]==REL_ENT['arg1']:
                    max_len=length
                    max_entity=nest_list[i]
                else:
                    if nest_list[i][2]==REL_ENT['arg2'] and max_entity[1] not in rel_entity2_id:
                        max_len=length
                        max_entity=nest_list[i] 
                        
            else:
                if nest_list[i][1] in rel_entity2_id:
                    max_len=length
                    max_entity=nest_list[i] 
                elif max_entity[2]==REL_ENT['arg1'] and nest_list[i][2]==REL_ENT['arg2']:
                    max_len=length
                    max_entity=nest_list[i]
    if final_entity==[]:
        final_entity.append(max_entity)
    return final_entity

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='python SA_Pubtator_Conll.py -i inputfile -o outputfile')
    parser.add_argument('--inputfile', '-i', help="inputfile in PubTator format",default='corpus/NLM-Gene2.2nd/merge.TrainingList.6.Sp.PubTator')
    parser.add_argument('--outputfile', '-o', help="outputfile in BIO-conll format",default='corpus/NLM-Gene2.2nd/merge.TrainingList.6.Sp.conll')
    args = parser.parse_args()
    
    infile=args.inputfile
    outfile=args.outputfile
    
    #tokenizer
    token_input=ssplit_token(infile)

    #filter nest entity
    nonest_input=corpus_noNest(token_input)

    # to conll
    generate_seq_input(nonest_input,outfile)