/**
 * Project: GNormPlus
 * Function: Gene Normalization
 */

package GNormPluslib;

import bioc.BioCAnnotation;
import bioc.BioCCollection;
import bioc.BioCDocument;
import bioc.BioCLocation;
import bioc.BioCPassage;

import bioc.io.BioCDocumentWriter;
import bioc.io.BioCFactory;
import bioc.io.woodstox.ConnectorWoodstox;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.BreakIterator;
import java.time.LocalDate;
import java.time.ZoneId;
import java.text.DecimalFormat;
import java.math.RoundingMode;

import javax.xml.stream.XMLStreamException;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

public class GN 
{
	public static HashMap<String, String> MatchedTokens_hash = new HashMap<String, String>();
	private double ScoringFunction(String geneid,HashMap<String,String> Mention_hash,String LF)
	{
		/*
		 * define gene/homo id
		 */
		
		//LF
		LF = LF.toLowerCase();
		LF = LF.replaceAll("([0-9])([a-z])", "$1 $2");
		LF = LF.replaceAll("([a-z])([0-9])", "$1 $2");
		LF = LF.replaceAll("([\\W\\-\\_])", " ");
		LF = LF.replaceAll("[ ]+", " ");
		String LF_tkn[]=LF.split(" ");
		int LF_ParticalMatch = 0;
		
		Pattern ptmp = Pattern.compile("[0-9]+\\-([0-9]+)");
		Matcher mtmp = ptmp.matcher(geneid);
		Pattern ptmp2 = Pattern.compile("([0-9]+)");
		Matcher mtmp2 = ptmp.matcher(geneid);
		if(mtmp.find())
		{
			geneid = "Homo:"+mtmp.group(1);
		}
		else
		{
			geneid = "Gene:"+geneid;
		}
		
		if(GNormPlus.GeneScoring_hash.containsKey(geneid))
		{
			HashMap<String,Double> TF = new HashMap<String,Double>(); // token i in gene j
			HashMap<String,Double> TermFrequency = new HashMap<String,Double>();
			
			/*
			 * Tokens in Query (Gene id lexicon)
			 */
			String l[]=GNormPlus.GeneScoring_hash.get(geneid).split("\t"); // Gene:2664293	cmk-1,cytidylate-1,kinase-1,mssa-1	0.4096	4	0.0625	1	2.0
			String tkns_Gene[] = l[0].split(",");
			for(int i=0;i<tkns_Gene.length;i++)
			{
				String Tkn_Freq[] = tkns_Gene[i].split("-");
				TermFrequency.put(Tkn_Freq[0], Double.parseDouble(Tkn_Freq[1]));
			}
			Double Cj =  Double.parseDouble(l[1]);
			Double AllTknNum = Double.parseDouble(l[2]);
			//Double Cj_max =  Double.parseDouble(l[3]);
			//Double MaxTknNum = Double.parseDouble(l[4]);
			Double Norm = Double.parseDouble(l[5]);
			if(Norm == 0.0){Norm=1.0;}
			
			/*
			 * Tokens in Document (recognized mentions)
			 */
			for(String Mention : Mention_hash.keySet())
			{
				Mention = Mention.toLowerCase();
				Mention = Mention.replaceAll("([0-9])([a-z])", "$1 $2");
				Mention = Mention.replaceAll("([a-z])([0-9])", "$1 $2");
				Mention = Mention.replaceAll("([\\W\\-\\_])", " ");
				Mention = Mention.replaceAll("[ ]+", " ");
				String tkns_Mention[]=Mention.split(" ");
				for(int i=0;i<tkns_Mention.length;i++)
				{
					if(TermFrequency.containsKey(tkns_Mention[i]))
					{
						TF.put(tkns_Mention[i], TermFrequency.get(tkns_Mention[i]));
					}
				}
			}
			
			Double score=0.0;
			for(String Tkn : TF.keySet())
			{
				//LF
				for(int t=0;t<LF_tkn.length;t++)
				{
					if(LF_tkn[t].equals(Tkn))
					{
						LF_ParticalMatch++;
					}
				}
				
				double TFij = TF.get(Tkn)/AllTknNum;
				double IDFi=GNormPlus.GeneScoringDF_hash.get(Tkn);
				score=score+TFij*IDFi*(1/(1-TFij));
			}
			//score = Cj * (1/Norm) *score;
			if(LF_ParticalMatch>0){score = score + LF_ParticalMatch;/*System.out.println(geneid+"\t"+LF+"\t"+score);*/}
			return score;
		}
		else
		{
			//System.out.println("Error: cannot find geneid: "+geneid+" in GeneScoring_hash");
			return 0.0;
		}
	}
	
	public void PreProcessing4GN(String Filename,String FilenameBioC) throws IOException, XMLStreamException
	{
		for (int i = 0; i < GNormPlus.BioCDocobj.Annotations.size(); i++) 
		{
			for (int j = 0; j < GNormPlus.BioCDocobj.Annotations.get(i).size(); j++) 
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++)
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
    				String start=anno[0];
					String last=anno[1];
					String mentions=anno[2];
					String type=anno[3];
					String id="";
					if(anno.length>=5)
					{
						id=anno[4];	
					}
					
					if(type.matches("Gene|GENERIF"))
					{
						String mentionArr[] = mentions.split("\\|");
						boolean update=false;
						for(int m=0;m<mentionArr.length;m++)
						{
							Pattern ptmp = Pattern.compile("^(.*[0-9A-Z])[ ]*p$");
							Matcher mtmp = ptmp.matcher(mentionArr[m]);
							Pattern ptmp2 = Pattern.compile("^(.+)nu$");
							Matcher mtmp2 = ptmp2.matcher(mentionArr[m]);
							Pattern ptmp3 = Pattern.compile("^(.*)alpha(.*)$");
							Matcher mtmp3 = ptmp3.matcher(mentionArr[m]);
							Pattern ptmp4 = Pattern.compile("^(.*)beta(.*)$");
							Matcher mtmp4 = ptmp4.matcher(mentionArr[m]);
							Pattern ptmp5 = Pattern.compile("^(.+[0-9])a$");
							Matcher mtmp5 = ptmp5.matcher(mentionArr[m]);
							Pattern ptmp6 = Pattern.compile("^(.+[0-9])b$");
							Matcher mtmp6 = ptmp6.matcher(mentionArr[m]);
							Pattern ptmp7 = Pattern.compile("^(.+)II([a-z])$");
							Matcher mtmp7 = ptmp7.matcher(mentionArr[m]);
							Pattern ptmp8 = Pattern.compile("^(.+)III([a-z])$");
							Matcher mtmp8 = ptmp8.matcher(mentionArr[m]);
							if(mtmp.find())
							{
								mentions=mentions+"|"+mtmp.group(1);
								update=true;
							}
							if(mtmp2.find())
							{
								mentions=mentions+"|"+mtmp2.group(1);
								update=true;
							}
							if(mtmp3.find())
							{
								mentions=mentions+"|"+mtmp3.group(1)+"a"+mtmp3.group(2);
								update=true;
							}
							if(mtmp4.find())
							{
								mentions=mentions+"|"+mtmp4.group(1)+"b"+mtmp4.group(2);
								update=true;
							}
							if(mtmp5.find())
							{
								mentions=mentions+"|"+mtmp5.group(1)+"alpha";
								update=true;
							}
							if(mtmp6.find())
							{
								mentions=mentions+"|"+mtmp6.group(1)+"beta";
								update=true;
							}
							if(mtmp7.find())
							{
								mentions=mentions+"|"+mtmp7.group(1)+"2"+mtmp7.group(2);
								update=true;
							}
							if(mtmp8.find())
							{
								mentions=mentions+"|"+mtmp8.group(1)+"3"+mtmp8.group(2);
								update=true;
							}
						}
						if(update == true)
						{
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start + "\t" + last + "\t" + mentions + "\t" + type + "\t" + id );
						}
					}
				}
			}
		}
		//GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,true);
	}
	
	public void ChromosomeRecognition(String Filename,String FilenameBioC) throws IOException, XMLStreamException
	{
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++) /** PMIDs : i */
		{
			String Pmid = GNormPlus.BioCDocobj.PMIDs.get(i);
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
				
				/** Chromosome recognition */
				ArrayList<String> locations = GNormPlus.PT_GeneChromosome.SearchMentionLocation(PassageContext,"ChromosomeLocation");
				for (int k = 0 ; k < locations.size() ; k++)
				{
					String anno[]=locations.get(k).split("\t");
					//int start= Integer.parseInt(anno[0]);
	        		//int last= Integer.parseInt(anno[1]);
	        		//String mention = anno[2];
	        		String ids = anno[3];
	        		//GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tChromosomeLocation\t"+ids); //paragraph
	        		String IDs[] = ids.split("[\\|,]");
	        		for(int idcount=0;idcount<IDs.length;idcount++)
	        		{
	        			//IDs[idcount] = IDs[idcount].replaceAll("\\-[0-9]+", "");
	        			GNormPlus.Pmid2ChromosomeGene_hash.put(Pmid+"\t"+IDs[idcount],"");
	        		}
				}
			}
		}
		//GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,true);
	}
	
	public void GeneNormalization(String Filename,String FilenameBioC,boolean GeneIDMatch) throws IOException, XMLStreamException
	{
		final DecimalFormat df = new DecimalFormat("0.####");
        df.setRoundingMode(RoundingMode.HALF_UP);
		
		//Tokenization
		for (int i = 0; i < GNormPlus.BioCDocobj.Annotations.size(); i++) /** PMIDs : i */
		{
			String Pmid = GNormPlus.BioCDocobj.PMIDs.get(i);
			
			/** Species */
			HashMap<String,String> Species_hash = new HashMap<String,String>();
			for (int j = 0; j < GNormPlus.BioCDocobj.Annotations.get(i).size(); j++) /** Paragraphs : j */
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) /** Annotation : k */
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
    				String mentions=anno[2];
					String type=anno[3];
					if(type.matches("(Species|Genus|Strain|CellLine|Cell)"))
					{
						Species_hash.put(mentions,"");
					}
				}
			}
			
			
			/*
			 * Collect Gene mentions :
			 * 
			 *  GeneMention-taxid	->	"ID" : geneid
			 *  					->	"type" : "Gene"
			 *  					->	start1-last1 : ""
			 *  					->	start2-last2 : ""
			 *  					->	start3-last3 : ""
			 */

			String tiabs="";
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageContexts.get(i).size(); j++) /** Paragraphs : j */
			{
				tiabs=tiabs+GNormPlus.BioCDocobj.PassageContexts.get(i).get(j).toLowerCase();
			}
			HashMap<String,HashMap<String,String>> GeneMention_hash = new HashMap<String,HashMap<String,String>>();
			HashMap<String,String> Mention_hash = new HashMap<String,String>();
			for (int j = 0; j < GNormPlus.BioCDocobj.Annotations.get(i).size(); j++) /** Paragraphs : j */
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) /** Annotation : k */
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
    				String start=anno[0];
					String last=anno[1];
					String mentions=anno[2];
					String type=anno[3];
					String taxids="Tax:9606";
					
					if(anno.length>=5)
					{
						taxids=anno[4];
					}
					String mentions_tmp=mentions.toLowerCase();
					mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
					mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
					taxids=taxids.replaceAll("(Focus|Right|Left|Prefix|Tax):","");
					if(taxids.equals(""))
					{
						taxids="9606";
					}
					/** Filtering */
					boolean found_filter = false;
					if(GNormPlus.Filtering_hash.containsKey(mentions_tmp)) // filtering
					{
						found_filter=true;
					}
					
					if(found_filter==false) //abbreviation
					{
						for(String f : GNormPlus.Filtering_WithLongForm_hash.keySet())
						{
							if( GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).matches(".*[\\t\\|]"+f+"\tGene.*") ||
								GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).matches(".*\\t"+f+"\\|[^\t]+\tGene.*")
									)
							{
								String lf=GNormPlus.Filtering_WithLongForm_hash.get(f);
								if(tiabs.matches(".*"+lf+".*"))
								{
									found_filter=true;
									break;
								}
							}
						}
					}
					
					if(found_filter==false)
					{
						if( GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).matches(".*[\\t\\|][a-z]\tGene.*") ||
								GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).matches(".*\\t[a-z]\\|[^\t]+\tGene.*") //32171191	Wuhan's
									)
						{
							found_filter=true;
			 
						}
					}
					
					if(found_filter == false)
					{
						if(type.matches("(Gene|GENERIF)"))
	    				{
							if(GeneMention_hash.containsKey(mentions+"\t"+taxids))
	    					{
	    						GeneMention_hash.get(mentions+"\t"+taxids).put(start+"\t"+last,"");
	    					}
	    					else 
	    					{
	    						HashMap<String,String> offset_hash = new HashMap<String,String>();
	    						offset_hash.put(start+"\t"+last,"");
	    						GeneMention_hash.put(mentions+"\t"+taxids, offset_hash);
	    						GeneMention_hash.get(mentions+"\t"+taxids).put("type", type);
	    						Mention_hash.put(mentions,type);
	    					}
	    				}
	    				else if(type.matches("(FamilyName|DomainMotif)"))
	    				{
	    					String GMs[]=mentions.split("\\|");
	    					for(int g=0;g<GMs.length;g++)
							{
								String mention = GMs[g];
								Mention_hash.put(mention,"FamilyDomain");
							}
	    				}
					}
					
				}
			}
			
			/*
			 * Gene id refinement:
			 *  1. Official name
			 *  2. only one gene
			 */
			HashMap<String,String> GuaranteedGene2ID = new HashMap<String,String>();
			HashMap<String,String> MultiGene2ID = new HashMap<String,String>();
			for(String GeneMentionTax : GeneMention_hash.keySet())
			{
				String GT[]=GeneMentionTax.split("\\t");
				String mentions=GT[0];
				String taxids=GT[1];
				String GMs[]=mentions.split("\\|");
				
				HashMap<String,String> taxids_hash = new HashMap<String,String>();
				String taxids_arr[]=taxids.split(",");
				for(int t=0;t<taxids_arr.length;t++)
				{
					taxids_hash.put(taxids_arr[t], "");
				}
				
				for(int ms=0;ms<GMs.length;ms++)
				{
					String mention = GMs[ms];
					String IDstr = GNormPlus.PT_Gene.MentionMatch(mention); /** searched by PT_Gene */
					String IDs[]=IDstr.split("\\|");
					
					/*
					 * printing the ambiguous gene mentions and candidates
					 */
					//String IDs_s[]=IDstr.split(",");
					//if(IDs_s.length>1)
					//{
					//	System.out.println(Pmid+"\t"+mention+"\t"+mentions+"\t"+IDstr);
					//}
					
					for(int c=0;c<IDs.length;c++)
					{
						String tax2ID[]=IDs[c].split(":"); // tax2ID[0] = taxid ; tax2ID[1] = geneids
						if(taxids_hash.containsKey(tax2ID[0]))
						{
							String geneid=tax2ID[1];
							String TargetTax=tax2ID[0];
							GeneMention_hash.get(GeneMentionTax).put("ID", geneid);
							GeneMention_hash.get(GeneMentionTax).put("TargetTax", TargetTax);
							break;
						}
					}
					
					//geneid refinement
					if(GeneMention_hash.get(GeneMentionTax).containsKey("ID"))
					{
						Pattern ptmp = Pattern.compile("\\*([0-9]+(\\-[0-9]+|))");
						Matcher mtmp = ptmp.matcher(GeneMention_hash.get(GeneMentionTax).get("ID"));
						
						if(mtmp.find()) // 1. Official Name
						{
							GeneMention_hash.get(GeneMentionTax).put("ID",mtmp.group(1));
							GuaranteedGene2ID.put(GeneMentionTax,mtmp.group(1));
						}
						else if(GeneMention_hash.get(GeneMentionTax).get("ID").matches("[0-9]+(\\-[0-9]+|)")) // 2. only one gene
						{
							GuaranteedGene2ID.put(GeneMentionTax,GeneMention_hash.get(GeneMentionTax).get("ID"));
						}
						else
						{
							String ID[] = GeneMention_hash.get(GeneMentionTax).get("ID").split(",");
							boolean FoundByChroLoca=false;
							for(int idcount=0;idcount<ID.length;idcount++)
							{
								if(GNormPlus.Pmid2ChromosomeGene_hash.containsKey(Pmid+"\t"+ID[idcount])) // 3. Chromosome location
								{
									GuaranteedGene2ID.put(GeneMentionTax,ID[idcount]);
									FoundByChroLoca=true;
									break;
								}
							}
							if(FoundByChroLoca == false)
							{
								MultiGene2ID.put(GeneMentionTax, GeneMention_hash.get(GeneMentionTax).get("ID"));
							}
						}
					}
					if(GNormPlus.suffixprefix_orig2modified.containsKey(mention) && (!IDstr.equals("-1")) && (!IDstr.equals("-2")) && (!IDstr.equals("-3")))
					{
						break;
					}
				}
			}
			
			/*
			 * Gene id refinement:
			 *  3. multiple genes but can be inferred by 1. and 2.
			 */
			for(String GeneMentionTax_M : MultiGene2ID.keySet())
			{
				for(String GeneMentionTax_G : GuaranteedGene2ID.keySet())
				{
					String MG[] = MultiGene2ID.get(GeneMentionTax_M).split(",");
					for(int m=0;m<MG.length;m++)
					{
						if(MG[m].equals(GuaranteedGene2ID.get(GeneMentionTax_G)))
						{
							GeneMention_hash.get(GeneMentionTax_M).put("ID",MG[m]);
						}
					}
				}
			}
			
			/*
			 * Gene id refinement:
			 *  4. FullName -> Abbreviation
			 */
			for(String GeneMentionTax : GeneMention_hash.keySet())
			{
				String MT[] = GeneMentionTax.split("\\t");
				if(GNormPlus.PmidLF2Abb_hash.containsKey(Pmid+"\t"+MT[0]))
				{
					String GeneMentionTax_Abb = GNormPlus.PmidLF2Abb_hash.get(Pmid+"\t"+MT[0]) + "\t" + MT[1];
					if(GeneMention_hash.containsKey(GeneMentionTax_Abb) && GeneMention_hash.get(GeneMentionTax).containsKey("ID"))
					{
						GeneMention_hash.get(GeneMentionTax_Abb).put("ID", GeneMention_hash.get(GeneMentionTax).get("ID"));
					}
				}
			}
			
			/*
			 * Gene id refinement:
			 *  5. Ranking by scoring function (inference network)
			 */
			for(String GeneMentionTax : GeneMention_hash.keySet())
			{
				if(GeneMention_hash.get(GeneMentionTax).containsKey("ID") && GeneMention_hash.get(GeneMentionTax).get("ID").matches(".+,.+"))
				{
					String geneids=GeneMention_hash.get(GeneMentionTax).get("ID");
					String geneid[] = geneids.split(",");
					
					String OutputStyle="Top1";
					if(OutputStyle.equals("Top1"))
					{
						//only return the best one
						double max_score=0.0;
						String target_geneid="";
						for(int g=0;g<geneid.length;g++)
						{
							String MT[] = GeneMentionTax.split("\\t");
							String LF="";
							if(GNormPlus.PmidAbb2LF_hash.containsKey(Pmid+"\t"+MT[0]))
							{
								LF = GNormPlus.PmidAbb2LF_hash.get(Pmid+"\t"+MT[0]);
							}
							double score = ScoringFunction(geneid[g],Mention_hash,LF);
							if(score>max_score)
							{
								max_score=score;
								target_geneid=geneid[g];
							}
							else if(score == 0.0)
							{
								//System.out.println(GeneMentionTax);
							}
						}
						GeneMention_hash.get(GeneMentionTax).put("ID", target_geneid);
					}
					else // "All"
					{
						//return all geneids
						String geneSTR="";
						for(int g=0;g<geneid.length;g++)
						{
							String MT[] = GeneMentionTax.split("\\t");
							String LF="";
							if(GNormPlus.PmidAbb2LF_hash.containsKey(Pmid+"\t"+MT[0]))
							{
								LF = GNormPlus.PmidAbb2LF_hash.get(Pmid+"\t"+MT[0]);
							}
							double score = ScoringFunction(geneid[g],Mention_hash,LF);
							String hoge = df.format(score);
							score=Double.parseDouble(hoge);
							
							if(geneSTR.equals(""))
							{
								geneSTR=geneid[g]+"-"+score;
							}
							else
							{
								geneSTR=geneSTR+","+geneid[g]+"-"+score;
							}
						}
						GeneMention_hash.get(GeneMentionTax).put("ID", geneSTR);
					}
				}
			}
			
			/*
			 * Gene id refinement: - removed (Reason: cause too much False Positive)
			 *  6. Abbreviation -> FullName
			 *  
			 */
			for(String GeneMentionTax : GeneMention_hash.keySet())
			{
				String MT[] = GeneMentionTax.split("\\t");
				if(GNormPlus.PmidAbb2LF_hash.containsKey(Pmid+"\t"+MT[0]))
				{
					String GeneMentionTax_LF = GNormPlus.PmidAbb2LF_hash.get(Pmid+"\t"+MT[0]) + "\t" + MT[1];
					if(GeneMention_hash.containsKey(GeneMentionTax_LF) && GeneMention_hash.get(GeneMentionTax).containsKey("ID"))
					{
						GeneMention_hash.get(GeneMentionTax_LF).put("ID", GeneMention_hash.get(GeneMentionTax).get("ID"));
					}
				}
			}
			
			/*
			 * Gene id refinement:
			 *  7. The inference network tokens of Abbreviation.ID should contain at least LF tokens
			 *  8. The short mention should be filtered if not long form support
			 */
			ArrayList<String> removeGMT = new ArrayList<String>();
			for(String GeneMentionTax : GeneMention_hash.keySet())
			{
				String GT[]=GeneMentionTax.split("\\t");
				String mentions=GT[0];
				String tax=GT[1];
				if(GeneMention_hash.get(GeneMentionTax).containsKey("type") && GeneMention_hash.get(GeneMentionTax).get("type").matches("Gene|GENERIF") && GeneMention_hash.get(GeneMentionTax).containsKey("ID"))
				{
					String type = GeneMention_hash.get(GeneMentionTax).get("type");
					String id = GeneMention_hash.get(GeneMentionTax).get("ID");
					String geneid="";
					Pattern ptmp1 = Pattern.compile("^([0-9]+)\\-([0-9]+)$");
					Pattern ptmp2 = Pattern.compile("^([0-9]+)$");
					Matcher mtmp1 = ptmp1.matcher(id);
					Matcher mtmp2 = ptmp2.matcher(id);
					//System.out.println(id);
					if(mtmp1.find())
					{
						geneid = "Homo:"+mtmp1.group(2);
					}
					else if(mtmp2.find())
					{
						geneid = "Gene:"+mtmp2.group(1);
					}
					
					boolean LongFormTknMatch= false;
					boolean LongFormExist= true;
					if(GNormPlus.GeneScoring_hash.containsKey(geneid))
					{
						if(GNormPlus.PmidAbb2LF_lc_hash.containsKey(Pmid+"\t"+mentions.toLowerCase()))
						{
							/*
							 * token in lexicon : tkn_lexicon
							 * token in mention : tkn_mention
							 */
							String l[]=GNormPlus.GeneScoring_hash.get(geneid).split("\t"); // Gene:2664293	cmk-1,cytidylate-1,kinase-1,mssa-1	0.4096	4	0.0625	1	2.0
							String tkns_Gene[] = l[0].split(",");
							ArrayList<String> tkn_lexicon = new ArrayList<String>();
							for(int ti=0;ti<tkns_Gene.length;ti++)
							{
								String Tkn_Freq[] = tkns_Gene[ti].split("-");
								tkn_lexicon.add(Tkn_Freq[0]);
							}
							
							String LF_lc=GNormPlus.PmidAbb2LF_lc_hash.get(Pmid+"\t"+mentions.toLowerCase());
							LF_lc = LF_lc.replaceAll("([0-9])([A-Za-z])", "$1 $2");
							LF_lc = LF_lc.replaceAll("([A-Za-z])([0-9])", "$1 $2");
							String tkn_mention[] = LF_lc.split("[\\W\\-\\_]");
							for(int tl=0;tl<tkn_lexicon.size();tl++)
							{
								for(int tm=0;tm<tkn_mention.length;tm++)
								{
									if(tkn_lexicon.get(tl).equals(tkn_mention[tm]) && (!tkn_mention[tm].matches("[0-9]+")))
									{
										LongFormTknMatch = true;
									}
								}	
							}
						}
						else{LongFormExist = false;}
					}
					else{LongFormTknMatch = true;} // exception
					
					if(LongFormTknMatch == false && LongFormExist == true) // 7.
					{
						removeGMT.add(GeneMentionTax); //remove short form
						removeGMT.add(GNormPlus.PmidAbb2LF_hash.get(Pmid+"\t"+mentions)+"\t"+tax); //remove long form
					}
					else if(mentions.length()<=2 && LongFormExist == false) // 8.
					{
						removeGMT.add(GeneMentionTax);
					}
				}
			}
			
			for(int gmti=0;gmti<removeGMT.size();gmti++) // remove
			{
				GeneMention_hash.remove(removeGMT.get(gmti));
			}
						
			// Append gene ids
			for (int j = 0; j < GNormPlus.BioCDocobj.Annotations.get(i).size(); j++) // Paragraphs : j
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
    				String start=anno[0];
					String last=anno[1];
					String mentions=anno[2];
					String type=anno[3];
					String taxid_org="Tax:9606";
					if(anno.length>=5)
					{
						taxid_org=anno[4];
					}
					String taxids=taxid_org.replaceAll("(Focus|Right|Left|Prefix|Tax):","");
					String GMs[]=mentions.split("\\|");
					
					if(GeneMention_hash.containsKey(mentions+"\t"+taxids) && GeneMention_hash.get(mentions+"\t"+taxids).containsKey("TargetTax"))
					{
						String taxtype=taxid_org.replaceAll(":([0-9,]+)","");
						String taxid=GeneMention_hash.get(mentions+"\t"+taxids).get("TargetTax");
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mentions+"\t"+type+"\t"+taxtype+":"+taxid);
					}
					
					if(type.matches("Gene|GENERIF"))
					{
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k) + "|");
						
					
						if(GeneMention_hash.containsKey(mentions+"\t"+taxids) && GeneMention_hash.get(mentions+"\t"+taxids).containsKey("ID"))
						{
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k) + GeneMention_hash.get(mentions+"\t"+taxids).get("ID") + "," );
						}
						else // cannot find appropriate species
						{
							//System.out.println(mention+"\t"+taxid);
						}
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).substring(0, GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).length()-1)); // remove ",$"
					}
				}
			}
			
			//Extend to all gene mentions
			HashMap<String,String> GeneMentions = new HashMap<String,String>(); // Extending Gene mentions
			HashMap<String,String> GeneMentionLocation = new HashMap<String,String>(); // Extending Gene mentions
			for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					int start = Integer.parseInt(anno[0]);
					int last = Integer.parseInt(anno[1]);
					String mentions=anno[2];
					String type=anno[3];
					String id="Tax:9606";
					if(anno.length>=5)
					{
						id=anno[4];
					}
					if(type.matches("Gene|GENERIF") && id.matches("(Focus|Right|Left|Prefix|Tax)\\:([0-9]+)\\|([0-9]+)\\-([0-9]+)"))
					{
						GeneMentions.put(mentions.toLowerCase(), id);
						for (int s=start ;s<=last;s++)
						{
							GeneMentionLocation.put(j+"\t"+s,"");
						}
					}
					else if(type.matches("Gene|GENERIF") && id.matches("(Focus|Right|Left|Prefix|Tax)\\:([0-9]+)\\|([0-9]+)"))
					{
						GeneMentions.put(mentions.toLowerCase(), id);
						for (int s=start ;s<=last;s++)
						{
							GeneMentionLocation.put(j+"\t"+s,"");
						}
					}
				}
			}
			for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>i && GNormPlus.BioCDocobj.PassageContexts.get(i).size()>j)
				{
					String PassageContexts = " " + GNormPlus.BioCDocobj.PassageContexts.get(i).get(j) + " ";
					String PassageContexts_tmp = PassageContexts.toLowerCase();
					for(String gm : GeneMentions.keySet())
					{
						String id = GeneMentions.get(gm);
						if(gm.length()>=3)
						{
							gm = gm.replaceAll("[ ]*[\\|]*$", "");
							gm = gm.replaceAll("^[\\|]*[ ]*", "");
							gm = gm.replaceAll("[\\|][\\|]+", "\\|");
							if(!gm.matches("[\\W\\-\\_]*"))
							{
								gm = gm.replaceAll("([^A-Za-z0-9\\| ])", "\\\\$1");
								Pattern ptmp = Pattern.compile("^(.*[\\W\\-\\_])("+gm+")([\\W\\-\\_].*)$");
								Matcher mtmp = ptmp.matcher(PassageContexts_tmp);
								while(mtmp.find())
								{
									String pre = mtmp.group(1);
									String gmtmp = mtmp.group(2);
									String post = mtmp.group(3);
			
									int start = pre.length()-1;
									int last = start+gmtmp.length();
									if(PassageContexts.length()>=last+1)
									{
										String mention = PassageContexts.substring(start+1,last+1);
										if(!GeneMentionLocation.containsKey(j+"\t"+start) && !GeneMentionLocation.containsKey(j+"\t"+last))
										{
											GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tGene\t"+id);
										}
									}
									gmtmp = gmtmp.replaceAll(".", "\\@");
									PassageContexts_tmp=pre+""+gmtmp+""+post;
									mtmp = ptmp.matcher(PassageContexts_tmp);
								}
							}
						}
					}
				}
			}
			
			//Apply to FamilyNames
			HashMap<String,String> geneids = new HashMap<String,String>(); // Extending Gene mentions
			for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					String type=anno[3];
					if(type.matches("(Gene|GENERIF)"))
					{
						String id="Tax:9606";
						if(anno.length>=5)
						{
							id=anno[4];
						}
						Pattern ptmp0 = Pattern.compile("^(Focus|Right|Left|Prefix|GeneID|Tax)\\:([0-9]+)\\|([0-9]+)$");
						Matcher mtmp0 = ptmp0.matcher(id);
						Pattern ptmp1 = Pattern.compile("^(Focus|Right|Left|Prefix|GeneID|Tax)\\:([0-9]+)\\|([0-9]+)\\-([0-9]+)$");
						Matcher mtmp1 = ptmp1.matcher(id);
						if(mtmp0.find())
						{
							geneids.put(mtmp0.group(3), "");
						}
						if(mtmp1.find())
						{
							geneids.put(mtmp1.group(3), "");
						}
					}
				}
			}
			for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
			{
				for (int k = GNormPlus.BioCDocobj.Annotations.get(i).get(j).size()-1; k >=0 ; k--) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					String mention=anno[2];
					String type=anno[3];
					if(type.matches("(FamilyName|DomainMotif)"))
					{
						String id="Tax:9606";
						if(anno.length>=5)
						{
							id=anno[4];
						}
						String IDstrs = GNormPlus.PT_FamilyName.MentionMatch(mention);
						String IDstr[]=IDstrs.split("\\|");
						String ids="";
						for(int id_i=0;id_i<IDstr.length;id_i++)
						{
							if(geneids.containsKey(IDstr[id_i]))
							{
								 if(ids.equals(""))
								 {
									 ids=IDstr[id_i];
								 }
								 else
								 {
									 ids=ids+";"+IDstr[id_i];
								 }
							}
						}
						if(!ids.equals(""))
						{
							if(type.equals("FamilyName")){type="Gene";}
							String Annotation_k=anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+type+"\tTax:9606";
							if(anno.length>=5)
							{
								Annotation_k=anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+type+"\t"+anno[4];
							}
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k,Annotation_k+"|"+ids);
						}
						else
						{
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).remove(k);
						}
					}
				}
			}
			//Species "*" and "(anti)" removed.
			for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
			{
				for (int k = GNormPlus.BioCDocobj.Annotations.get(i).get(j).size()-1; k >=0 ; k--) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					String type=anno[3];
					if(type.equals("Species") || type.equals("Genus") || type.equals("Strain") || type.equals("CellLine") || type.equals("Cell"))
					{
						String id=anno[4];
						id=id.replaceAll("\\*", "");
						id=id.replaceAll("\\(anti\\)", "");
						String Annotation_k=anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+type+"\t"+id;
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k,Annotation_k);
					}
				}
			}
			
			for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
			{
				
				for (int k = GNormPlus.BioCDocobj.Annotations.get(i).get(j).size()-1; k >=0 ; k--) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					int start = Integer.parseInt(anno[0]);
					int last = Integer.parseInt(anno[1]);
					String mention = anno[2];
					String type = anno[3];
					String id = anno[4];
					if(type.matches("(Gene|GENERIF)") && Species_hash.containsKey(mention))
					{
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).remove(k);
					}
					else if(type.matches("(Gene|GENERIF)") && id.equals(""))
					{
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).remove(k);
					}
					else
					{
						for (int k1 = GNormPlus.BioCDocobj.Annotations.get(i).get(j).size()-1; k1 >=0 ; k1--) // Annotation : k
						{
							if(k1 != k)
							{
								String anno1[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k1).split("\t");
								int start1 = Integer.parseInt(anno1[0]);
								int last1 = Integer.parseInt(anno1[1]);
								if((start1<start && last1>=last) || (start1<=start && last1>last))
								{
									GNormPlus.BioCDocobj.Annotations.get(i).get(j).remove(k);
									break;
								}
							}
						}
					}
				}
			}
		}
		if(GeneIDMatch == true)
		{
			//GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,true);
		}
		else
		{
			GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,true,true);
		}
	}
	/*
	 * Search Potential GeneID in the Prefix Tree
	 */
	public ArrayList<String> SearchGeneIDLocation(String Doc)
	{
		ArrayList<String> location = new ArrayList<String>();
		
		String Doc_tmp=" "+Doc+" ";
		Pattern ptmp = Pattern.compile("^(.*[^A-Za-z0-9]+)([0-9]+\\S*[A-Za-z]+|[A-Za-z]+\\S*[0-9]+|[0-9]+\\S*[A-Za-z]+\\S*[0-9]+|[A-Za-z]+\\S*[0-9]+\\S*[A-Za-z]+)([^A-Za-z0-9]+.*)$");
		Matcher mtmp = ptmp.matcher(Doc_tmp);
		while(mtmp.find())
		{
			String str1=mtmp.group(1);
			String str2=mtmp.group(2);
			String str3=mtmp.group(3);
			for(int m=str1.length();m<=(str1.length()+str2.length());m++)
			{
				int start = str1.length()-1;
				int last = start+str2.length();
				String mention = Doc.substring(start, last);
				if(!mention.matches(".*[\\'\\;\\[\\]\\+\\*\\\\].*"))
				{
					if(last-start>6 && (mention.matches(".*\\(.*\\).*") || mention.matches("[^\\(\\)]+")) )
					{
						Pattern ptmp1 = Pattern.compile("^(.+[^0-9])([0-9]+)\\-([0-9]+)$");
						Matcher mtmp1 = ptmp1.matcher(mention);
						Pattern ptmp2 = Pattern.compile("^(.+[^0-9])([0-9]+)\\-(.+[^0-9])([0-9]+)$");
						Matcher mtmp2 = ptmp2.matcher(mention);
						if(mtmp1.find())
						{
							String S1 = mtmp1.group(1);
							if(mtmp1.group(2).length()<=6 && mtmp1.group(3).length()<=6)
							{
								int Num1 = Integer.parseInt(mtmp1.group(2));
								int Num2 = Integer.parseInt(mtmp1.group(3));
								String prefix = "";
								Pattern ptmp3 = Pattern.compile("^([0]+)");
								Matcher mtmp3 = ptmp3.matcher(mtmp1.group(2));
								if(mtmp3.find())
								{
									prefix = mtmp3.group(1);
								}
								if(Num2-Num1>0 && (Num2-Num1<=20))
								{
									for(int n=Num1;n<=Num2;n++)
									{
										String StrNum=S1+prefix+n;
										if(StrNum.length()>=5)
										{
											location.add(start+"\t"+last+"\t"+StrNum+"\tGeneID");
										}
									}
								}
							}
						}
						else if(mtmp2.find())
						{
							if(mtmp2.group(2).length()<=6 && mtmp2.group(4).length()<=6)
							{
								String S1 = mtmp2.group(1);
								int Num1 = Integer.parseInt(mtmp2.group(2));
								String S2 = mtmp2.group(3);
								int Num2 = Integer.parseInt(mtmp2.group(4));
								if(S1.equals(S2))
								{
									String prefix = "";
									Pattern ptmp3 = Pattern.compile("^([0]+)");
									Matcher mtmp3 = ptmp3.matcher(mtmp2.group(2));
									if(mtmp3.find())
									{
										prefix = mtmp3.group(1);
									}
									if(Num2-Num1>0 && (Num2-Num1<=20))
									{
										for(int n=Num1;n<=Num2;n++)
										{
											String StrNum=S1+prefix+n;
											if(StrNum.length()>=5)
											{
												location.add(start+"\t"+last+"\t"+StrNum+"\tGeneID");
											}
										}
									}
								}
							}
						}
					}
					location.add(start+"\t"+last+"\t"+mention+"\tGeneID");
				}
			}
			String men="";
			for(int m=0;m<str2.length();m++){men=men+"@";}
			Doc_tmp=str1+men+str3;
			mtmp = ptmp.matcher(Doc_tmp);
		}
		return location;
	}
	public void GeneIDRecognition(String Filename,String FilenameBioC) throws IOException, XMLStreamException
	{
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++) /** PMIDs : i */
		{
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
				/** GeneID recognition by pattern match */
				ArrayList<String> locations = SearchGeneIDLocation(PassageContext);
				for (int k = 0 ; k < locations.size() ; k++)
				{
					String anno[]=locations.get(k).split("\t");
					String mention = anno[2].toLowerCase();
	        		mention = mention.replaceAll("[\\W\\-\\_]+", "");
	        		if(GNormPlus.GeneIDs_hash.containsKey(mention))
	        		{
	        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(locations.get(k)+"\tGeneID:"+GNormPlus.GeneIDs_hash.get(mention)); //paragraph
	        		}
				}
			}
		}
		GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,true,true);
	}
}