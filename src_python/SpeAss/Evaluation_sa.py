# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 15:33:54 2021

@author: luol2
"""
# compute metrics using IO prefile
#ignore arg1
def Rel_Evaluation(prefile):
    fin=open(prefile,'r',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    TP=0 #gold=pre=pos
    FP=0 #gold=neg, pre=pos
    FN=0 #gold=pos, pre=Neg
    for sentence in all_in:
        tokens=sentence.split('\n')
        entity_id=0
        token_id=0
        temp_gold='O'
        temp_pre='O'
        while (token_id<len(tokens)):
            seg=tokens[token_id].split('\t')
            if seg[0]=='<GENE>':
                if seg[1]=='O':
                    temp_gold=seg[1]
                else:
                    temp_gold=seg[1][2:]
                if seg[2]=='O':
                    temp_pre=seg[2]
                else:
                    temp_pre=seg[2][2:]  
                token_id+=1
                seg=tokens[token_id].split('\t')
                while seg[0]!='</GENE>':
                    token_id+=1
                    seg=tokens[token_id].split('\t')
                    if seg[1]!='O' and temp_gold=='O':
                        temp_gold=seg[1][2:]
                    if seg[2]!='O' and temp_pre=='O':
                        temp_pre=seg[2][2:]
                if temp_pre!='O' and temp_gold!='O' and temp_pre==temp_gold:
                    TP+=1
                elif temp_pre!='O' and temp_gold!='O' and temp_pre!=temp_gold:
                    FP+=1
                    FN+=1
                elif temp_pre!='O' and temp_gold=='O' :
                    FP+=1
                elif temp_pre=='O' and temp_gold!='O' :
                    FN+=1
                temp_pre='O'
                temp_gold='O'
        
            else:
                pass
            token_id+=1
    # print('TP,FP,FN:',TP,FP,FN)
    if TP+FP==0:
        P=0
    else:
        P=TP/(TP+FP)
    if TP+FN==0:
        R=0
    else:
        R=TP/(TP+FN)
    if P+R==0:
        F1=0
    else:
        F1=2*P*R/(P+R)
    print('TP,FP,FN:',TP,FP,FN)
    print('P,R,F1:',P,R,F1)


def Rel_Evaluation_fn(prefile):
    fin=open(prefile,'r',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    TP=0 #gold=pre=pos
    FP=0 #gold=neg, pre=pos
    FN=0 #gold=pos, pre=Neg
    for sentence in all_in:
        tokens=sentence.split('\n')
        entity_id=0
        token_id=0
        temp_gold='O'
        temp_pre='O'
        while (token_id<len(tokens)):
            seg=tokens[token_id].split('\t')
            if seg[0]=='<GENE>':
                if seg[1]=='O':
                    temp_gold=seg[1]
                else:
                    temp_gold=seg[1][2:]
                if seg[2]=='O':
                    temp_pre=seg[2]
                else:
                    temp_pre=seg[2][2:]  
                token_id+=1
                seg=tokens[token_id].split('\t')
                while seg[0]!='</GENE>':
                    token_id+=1
                    seg=tokens[token_id].split('\t')
                    if seg[1]!='O' and temp_gold=='O':
                        temp_gold=seg[1][2:]
                    if seg[2]!='O' and temp_pre=='O':
                        temp_pre=seg[2][2:]
                if temp_pre!='O' and temp_gold!='O' and temp_pre==temp_gold:
                    TP+=1
                elif temp_pre!='O' and temp_gold!='O' and temp_pre!=temp_gold:
                    FP+=1
                elif temp_pre!='O' and temp_gold=='O' :
                    FP+=1
                elif temp_pre=='O' and temp_gold!='O' :
                    FN+=1
                temp_pre='O'
                temp_gold='O'
        
            else:
                pass
            token_id+=1
    print('TP,FP,FN:',TP,FP,FN)
    if TP+FP==0:
        P=0
    else:
        P=TP/(TP+FP)
    if TP+FN==0:
        R=0
    else:
        R=TP/(TP+FN)
    if P+R==0:
        F1=0
    else:
        F1=2*P*R/(P+R)
    # print('TP,FP,FN:',TP,FP,FN)
    print('P,R,F1:',P,R,F1)
    return F1

def Rel_Evaluation_Hugface_fn(prefile,ARG2_label='gene1s'):
    fin=open(prefile,'r',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    TP=0 #gold=pre=pos
    FP=0 #gold=neg, pre=pos
    FN=0 #gold=pos, pre=Neg
    result_dict={}#{'rel type':[TP,FP,FN],...,}
    for sentence in all_in:
        tokens=sentence.split('\n')
        for token in tokens:
            seg=token.split('\t')
            if seg[0]==ARG2_label:
                if seg[1].find('ARG2')>=0:
                    if seg[2]==seg[1]:
                        if seg[1] not in result_dict.keys():
                            result_dict[seg[1]]=[1,0,0]
                        else:
                            result_dict[seg[1]][0]+=1
                        TP+=1
                    elif seg[2].find('ARG2')>=0:
                        if seg[1] not in result_dict.keys():
                            result_dict[seg[1]]=[0,0,1]
                        else:
                            result_dict[seg[1]][2]+=1
                        if seg[2] not in result_dict.keys():
                            result_dict[seg[2]]=[0,1,0]
                        else:
                            result_dict[seg[2]][1]+=1
                        FP+=1
                        FN+=1
                    else:
                        if seg[1] not in result_dict.keys():
                            result_dict[seg[1]]=[0,0,1]
                        else:
                            result_dict[seg[1]][2]+=1
                        FN+=1

                else:
                    if seg[2].find('ARG2')>=0:
                        if seg[2] not in result_dict.keys():
                            result_dict[seg[2]]=[0,1,0]
                        else:
                            result_dict[seg[2]][1]+=1
                        FP+=1
    # print('TP,FP,FN:',TP,FP,FN)
    rel_metrics={}
    for rel_type in result_dict.keys():
        if result_dict[rel_type][0]+result_dict[rel_type][1]==0:
            p=0
        else:
            p=result_dict[rel_type][0]/(result_dict[rel_type][0]+result_dict[rel_type][1])
        if result_dict[rel_type][0]+result_dict[rel_type][2]==0:
            r=0
        else:
            r=result_dict[rel_type][0]/(result_dict[rel_type][0]+result_dict[rel_type][2])
        if p+r==0:
            f1=0
        else:
            f1=2*p*r/(p+r)
        rel_metrics[rel_type]=[round(p,4),round(r,4),round(f1,4)]
    if TP+FP==0:
        P=0
    else:
        P=TP/(TP+FP)
    if TP+FN==0:
        R=0
    else:
        R=TP/(TP+FN)
    if P+R==0:
        F1=0
    else:
        F1=2*P*R/(P+R)
    P=round(P,4)
    R=round(R,4)
    F1=round(F1,4)
    print('mertics:\n',rel_metrics)
    print('\nTP,FP,FN:',TP,FP,FN)
    print('Overall P,R,F1:',P,R,F1)   
    return [P,R,F1],rel_metrics

def Rel_Evaluation_AIO_fn(prefile):
    fin=open(prefile,'r',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    TP=0 #gold=pre=pos
    FP=0 #gold=neg, pre=pos
    FN=0 #gold=pos, pre=Neg
    for sentence in all_in:
        tokens=sentence.split('\n')
        for token in tokens:
            seg=token.split('\t')
            if seg[0]=='<GENE>':
                if seg[1].find('ARG2-')>=0:
                    if seg[2]==seg[1]:
                        TP+=1
                    elif seg[2].find('ARG2-')>=0:
                        FP+=1
                        FN+=1
                    else:
                        FN+=1

                else:
                    if seg[2].find('ARG2-')>=0:
                        FP+=1
    # print('TP,FP,FN:',TP,FP,FN)
    if TP+FP==0:
        P=0
    else:
        P=TP/(TP+FP)
    if TP+FN==0:
        R=0
    else:
        R=TP/(TP+FN)
    if P+R==0:
        F1=0
    else:
        F1=2*P*R/(P+R)
    P=round(P,4)
    R=round(R,4)
    F1=round(F1,4)
    print('TP,FP,FN:',TP,FP,FN)
    print('P,R,F1:',P,R,F1)   
    return [P,R,F1]

def Rel_Evaluation_AIO_GC_fn(prefile):
    fin=open(prefile,'r',encoding='utf-8')
    all_in=fin.read().strip().split('\n\n')
    fin.close()
    TP=0 #gold=pre=pos
    FP=0 #gold=neg, pre=pos
    FN=0 #gold=pos, pre=Neg
    for sentence in all_in:
        tokens=sentence.split('\n')
        for token in tokens:
            seg=token.split('\t')
            if seg[0]=='<CHEMICAL>':
                if seg[1].find('ARG2-')>=0:
                    if seg[2]==seg[1]:
                        TP+=1
                    elif seg[2].find('ARG2-')>=0:
                        FP+=1
                        FN+=1
                    else:
                        FN+=1

                else:
                    if seg[2].find('ARG2-')>=0:
                        FP+=1
    # print('TP,FP,FN:',TP,FP,FN)
    if TP+FP==0:
        P=0
    else:
        P=TP/(TP+FP)
    if TP+FN==0:
        R=0
    else:
        R=TP/(TP+FN)
    if P+R==0:
        F1=0
    else:
        F1=2*P*R/(P+R)
    P=round(P,4)
    R=round(R,4)
    F1=round(F1,4)
    print('TP,FP,FN:',TP,FP,FN)
    print('P,R,F1:',P,R,F1)   
    return [P,R,F1]
    
def office_evaluation(goldfile,prefile):
    fin_gold=open(goldfile,'r',encoding='utf-8')
    all_gold=fin_gold.read().strip().split('\n')
    fin_gold.close()
    fin_pre=open(prefile,'r',encoding='utf-8')
    all_pre=fin_pre.read().strip().split('\n')
    fin_pre.close()
    
    gold_result={}#{'relation type':set(line)}
    pre_result={}
    all_result={} #{'relation type':[tp,fp,fn]}
    for line in all_gold:
        seg=line.split('\t')
        if seg[1] not in all_result.keys():
            all_result[seg[1]]=[0,0,0]
        if seg[1] not in gold_result.keys(): 
            gold_result[seg[1]]=set()
            gold_result[seg[1]].add(line)         
        else:
            gold_result[seg[1]].add(line)
      
    for line in all_pre:
        seg=line.split('\t')
        if seg[1] not in pre_result.keys():
            pre_result[seg[1]]=set()
            pre_result[seg[1]].add(line)
        else:
            pre_result[seg[1]].add(line)
            
    for rel_type in gold_result.keys():
        for gold_ele in gold_result[rel_type]:
            if rel_type not in pre_result.keys():
                all_result[rel_type][2]+=1
            else:
                if gold_ele in pre_result[rel_type]:
                    all_result[rel_type][0]+=1
                else:
                    all_result[rel_type][2]+=1
        if rel_type in pre_result.keys():
            for pre_ele in pre_result[rel_type]:
                if pre_ele not in gold_result[rel_type]:
                    all_result[rel_type][1]+=1
    ave_f=0
    TP,FP,FN=0,0,0
    print(all_result)
    for rel_type in all_result.keys():
        TP+=all_result[rel_type][0]
        FP+=all_result[rel_type][1]
        FN+=all_result[rel_type][2]
        tem_p,tem_r,tem_f=0,0,0
        if all_result[rel_type][0]+all_result[rel_type][1]==0:
            tem_p=0
        else:
            tem_p=all_result[rel_type][0]/(all_result[rel_type][0]+all_result[rel_type][1])
        if all_result[rel_type][0]+all_result[rel_type][2]==0:
            tem_r=0
        else:
            tem_r=all_result[rel_type][0]/(all_result[rel_type][0]+all_result[rel_type][2])
        if tem_p+tem_r==0:
            tem_f=0
        else:
            tem_f=2*tem_p*tem_r/(tem_p+tem_r)
        ave_f+=tem_f
        print('%s:p=%.4f,r=%.4f,f=%.4f' % (rel_type,tem_p,tem_r,tem_f))
    
    if TP+FP==0:
        P=0
    else:
        P=TP/(TP+FP)
    if TP+FN==0:
        R=0
    else:
        R=TP/(TP+FN)
    if P+R==0:
        F1=0
    else:
        F1=2*P*R/(P+R)
    ave_f+=tem_f
    
    print('Overall:')
    print('ave_f1:',ave_f/len(all_result))
    print('TP=%d, FP=%d, FN=%d'%(TP,FP,FN))
    print('P=%.4f, R=%.4f, F1=%.4f'%(P,R,F1))
    

if __name__=='__main__':
    path='//panfs/pan1/bionlplab/luol2/BC7DrugProt/results/'
    office_evaluation(path+'dev/dev_gold_relations.tsv',path+'drugprot_dev_LSTM-CRF-ES_pre.tsv')
    print('............')
    Rel_Evaluation_check('//panfs/pan1/bionlplab/luol2/BC7DrugProt/check/dev_pre_temp.conll')
