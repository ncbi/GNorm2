/**
 * Project: GNormPlus
 * Function: Species recognition and Species assignment
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

import javax.xml.stream.XMLStreamException;

import org.tartarus.snowball.SnowballStemmer;
import org.tartarus.snowball.ext.englishStemmer;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Collections;

public class SR 
{
	@SuppressWarnings("null")
	public void SpeciesRecognition(String Filename,String FilenameBioC,String StrainFilename,String FilterAntibody) throws IOException, XMLStreamException
	{
		/** Recognizing Species Names: SP */
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++) /** PMIDs : i */
		{
			String Pmid = GNormPlus.BioCDocobj.PMIDs.get(i);
			PrefixTree PT_Genus = new PrefixTree();
			HashMap<String, String> SPID_hash = new HashMap<String, String>();
			ArrayList<String> TargetedLocation = new ArrayList<String>();
			HashMap<String, String> GenusNames = new HashMap<String, String>();
			HashMap<String, String> Mention2ID_lc = new HashMap<String, String>();
			ArrayList<String> IDset = new ArrayList<String>();
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
				
				/** Species recognition */
				ArrayList<String> locations = GNormPlus.PT_Species.SearchMentionLocation(PassageContext,"Species"); /** PT_Species */
				for (int k = 0 ; k < locations.size() ; k++)
				{
					String anno[]=locations.get(k).split("\t");
					int start= Integer.parseInt(anno[0]);
	        		int last= Integer.parseInt(anno[1]);
	        		
	        		// For anti-serum filtering
	        		String ForwardSTR="";
	        		String BackwardSTR="";
	        		if(start>21)
	        		{
	        			ForwardSTR = (PassageContext+"ZZZZZZZZZZZZZZZZZZZZZZZZZZZ").substring(start-21,last);
	        		}
	        		else
	        		{
	        			ForwardSTR = (PassageContext+"ZZZZZZZZZZZZZZZZZZZZZZZZZZZ").substring(0,last);
	        		}
	        		if(PassageContext.length()>last+21)
	        		{
	        			BackwardSTR = PassageContext.substring(start,last+21);
	        		}
	        		else
	        		{
	        			BackwardSTR = PassageContext.substring(start,PassageContext.length());
	        		}
	        			        		
	        		String mention = anno[2];
	        		String id = anno[3];
	        		String mention_tmp=mention.toLowerCase();
	        		mention_tmp = mention_tmp.replaceAll("([^A-Za-z0-9@ ])", "\\\\$1");
	        		String antibody="";
	        		if(ForwardSTR.toLowerCase().matches(".*(anti|antibody|antibodies|serum|polyclonal|monoclonal|igg)[\\W\\-\\_]+"+mention_tmp)) {antibody="(anti)";}//filtering : antibody
	        		else if(BackwardSTR.toLowerCase().matches(mention_tmp+"[\\W\\-\\_]+(anti|antibody|antibodies|serum|polyclonal|monoclonal|igg).*")){antibody="(anti)";} //filtering : antibody
	        		else if(BackwardSTR.toLowerCase().matches(mention_tmp+"[\\W\\-\\_]+[A-Za-z0-9]+[\\W\\-\\_]+(anti|antibody|antibodies|serum|polyclonal|monoclonal|igg).*")){antibody="(anti)";} //filtering : antibody
	        		
					if(mention.matches(".*[\\(\\[\\{].*") && BackwardSTR.toLowerCase().matches(mention_tmp+"\\).*") )
    				{
	        			last=last+1;
	        			mention=mention+")";
    				}
	        		
	        		if(BackwardSTR.toLowerCase().matches(mention_tmp+"[0-9].*")){} // filtered: Bee1p
	        		else if((mention.matches(".*[;:,].*")) && mention.length()<=10){} // filtered : x, XXX
	        		else if(mention.matches("to[\\W\\-\\_]+[0-9]+")){} // to 7
	        		else if(mention.matches("[a-z][\\)\\]\\}].*") && (!mention.matches(".*[\\(\\[\\{].*")) && mention.length()<=10){} // s). Major
	        		else if(mention.matches(".*[\\(\\[\\{].*") && (!mention.matches(".*[\\)\\]\\}].*")) && mention.length()<=10){} // s). Major
	        		else if(!id.equals("NA"))
	        		{
	        			if(GNormPlus.BioCDocobj.Annotations.size()>i && GNormPlus.BioCDocobj.Annotations.get(i).size()>j)
						{
	        				if((!mention.matches("^[A-Za-z] [A-Za-z0-9]+$")) && (mention.length()>=3)) // invalid species: "a group/a GAL4/a strain"
	        				{
	        					if(FilterAntibody.equals("False") || (!antibody.equals("(anti)")))
	        					{
	        						String patt="^(.+?) [sS]train";
									Pattern ptmp = Pattern.compile(patt);
									Matcher mtmp = ptmp.matcher(mention);
									if(mtmp.find())
									{
										mention=mtmp.group(1);
										last=last-7;
									}
			        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tSpecies\t"+id); //+antibody
			        				String mentions_tmp=mention.toLowerCase();
		    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
		    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
		    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
			        				Mention2ID_lc.put(mention.toLowerCase(), id); //+antibody
			        				
			         				String mention_genus = "";
			        				patt="^([A-Za-z]+) ";
									ptmp = Pattern.compile(patt);
									mtmp = ptmp.matcher(mention);
									if(mtmp.find())
									{
										mention_genus=mtmp.group(1); // get genus
									}
									
			        				IDset.add(id);
			        				for(int s=start;s<last;s++)
									{
										TargetedLocation.add(j+"\t"+s);
									}
			        				String ids[]=id.split(";");
			        				for(int x=0;x<ids.length;x++)
			        				{
				        				patt="^\\**([0-9]+)";
										ptmp = Pattern.compile(patt);
										mtmp = ptmp.matcher(ids[x]);
										if(mtmp.find())
										{
											SPID_hash.put(mtmp.group(1), mention_genus);
										}
			        				}
	        					}
							}
						}
	        		}
				}
				
				/** Cell Line recognition */
				locations = GNormPlus.PT_Cell.SearchMentionLocation(PassageContext,"Cell"); /** PT_Cell */
				for (int k = 0 ; k < locations.size() ; k++)
				{
					String anno[]=locations.get(k).split("\t");
					int start= Integer.parseInt(anno[0]);
	        		int last= Integer.parseInt(anno[1]);
	        		String mention = anno[2];
	        		String id = anno[3];
	        		if(GNormPlus.BioCDocobj.Annotations.size()>i && GNormPlus.BioCDocobj.Annotations.get(i).size()>j)
					{
	        			if(!TargetedLocation.contains(j+"\t"+start)) //already exists
	        			{
	        				int last40=0;
		        			if(PassageContext.length()>=last+40)
		        			{
		        				last40=last+40;
		        			}
		        			else
		        			{
		        				last40=PassageContext.length();
		        			}
		        			
		        			// For anti-serum filtering
			        		String ForwardSTR="";
			        		String BackwardSTR="";
			        		if(start>21)
			        		{
			        			ForwardSTR = PassageContext.substring(start-21,last);
			        		}
			        		else
			        		{
			        			ForwardSTR = PassageContext.substring(0,last);
			        		}
			        		if(PassageContext.length()>last+21)
			        		{
			        			BackwardSTR = PassageContext.substring(start,last+21);
			        		}
			        		else
			        		{
			        			BackwardSTR = PassageContext.substring(start,PassageContext.length());
			        		}
			        		String mention_tmp=mention.toLowerCase();
			        		mention_tmp = mention_tmp.replaceAll("([^A-Za-z0-9@ ])", "\\\\$1");
			        		if(mention_tmp.matches(".*[\\[\\]\\(\\)\\{\\}].*")){}
			        		else if(BackwardSTR.toLowerCase().matches(mention_tmp+"[0-9\\-\\_].*")){} // filtered: Bee1p
			        		else if(ForwardSTR.toLowerCase().matches(".*[0-9\\-\\_]"+mention_tmp)){} // filtered: IL-22RA1
			        		else
			        		{
			        			String patt="[\\W\\-]cell([\\- ]*line|)[s]*[\\W\\-]";
			    				Pattern ptmp = Pattern.compile(patt);
			    				Matcher mtmp = ptmp.matcher(PassageContext.substring(last, last40).toLowerCase());
			    				if(mtmp.find())
			    				{
			    					if(GNormPlus.taxid4gene.contains(id)) // for gene
				        			{
				        				id="*"+id;
					        		}
			    					GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tCell\t"+id);
									String mentions_tmp=mention.toLowerCase();
		    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
		    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
		    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
			    					IDset.add(id);
			    					for(int s=start;s<last;s++)
									{
										TargetedLocation.add(j+"\t"+s);
									}
			    				}
			        		}
		        		}
					}
				}
				
				/** Genus names*/
				for(String ID: SPID_hash.keySet())
				{
					if(GNormPlus.GenusID_hash.containsKey(ID))
					{
						GenusNames.put(ID,GNormPlus.GenusID_hash.get(ID));
					}
					if(SPID_hash.get(ID).length()>=7)
					{
						GenusNames.put(ID,SPID_hash.get(ID));
					}
				}
			}
			
			GenusNames.put("3702", "arabidopsis");
			GenusNames.put("4932", "saccharomyces");
			GenusNames.put("562", "escherichia");
			GenusNames.put("7227", "drosophila");
			GenusNames.put("8355", "xenopus");
			
			PT_Genus.Hash2Tree(GenusNames);
			
			/** Genus recognition */
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>i && 
					GNormPlus.BioCDocobj.PassageContexts.get(i).size()>j &&  
					GNormPlus.BioCDocobj.Annotations.size()>i &&  
					GNormPlus.BioCDocobj.Annotations.get(i).size()>j 
					)
				{
					String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j);
					ArrayList<String> locations_Genus = PT_Genus.SearchMentionLocation(PassageContext,"Genus"); /** PT_Genus*/
					for (int k = 0 ; k < locations_Genus.size() ; k++)
					{
						String anno[]=locations_Genus.get(k).split("\t");
						String start= anno[0];
		        		String last= anno[1];
		        		String mention = anno[2];
		        		String id = anno[3];
		        		if(!TargetedLocation.contains(j+"\t"+start)) //already exists
	        			{
		        			String patt="^\\**([0-9]+)$";
							Pattern ptmp = Pattern.compile(patt);
							Matcher mtmp = ptmp.matcher(id);
							if(mtmp.find())
							{
								id = mtmp.group(1);
							}
							
							if(GNormPlus.taxid4gene.contains(id)) // for gene
		        			{
		        				id="*"+id;
			        		}
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tGenus\t"+id);
							String mentions_tmp=mention.toLowerCase();
    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
							IDset.add(id);
							for(int s=Integer.parseInt(start);s<Integer.parseInt(last);s++)
							{
								TargetedLocation.add(j+"\t"+s);
							}
		        		}
		        	}
				}
			}
			
			/** Strain Tree */
			PrefixTree PT_Strain = new PrefixTree();
			HashMap<String, String> StrainID_hash = new HashMap<String, String>();
			BufferedReader br = new BufferedReader(new FileReader(StrainFilename));
			String line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				String ancestor = l[0];
				String tax_id = l[1];
				String tax_names = l[2];
				if(SPID_hash.containsKey(ancestor))
				{
					StrainID_hash.put(tax_id, tax_names); // tax id -> strain
				}
				else if(SPID_hash.containsKey(tax_id))
				{
					StrainID_hash.put(tax_id, tax_names); // tax id -> strain
				}
			}
			br.close();
			HashMap<String, String> StrainNames = new HashMap<String, String>();
			for(String ID: StrainID_hash.keySet())
			{
				StrainNames.put(ID,StrainID_hash.get(ID));
			}
			
			PT_Strain.Hash2Tree(StrainNames);
			
			/** Strain recognition */
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>i && 
					GNormPlus.BioCDocobj.PassageContexts.get(i).size()>j &&  
					GNormPlus.BioCDocobj.Annotations.size()>i &&  
					GNormPlus.BioCDocobj.Annotations.get(i).size()>j 
					)
				{
					String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
					ArrayList<String> locations_Strain = PT_Strain.SearchMentionLocation(PassageContext,"Strain"); /** PT_Strain*/
					for (int k = 0 ; k < locations_Strain.size() ; k++)
					{
						String anno[]=locations_Strain.get(k).split("\t");
						String start= anno[0];
		        		String last= anno[1];
		        		String mention = anno[2];
		        		String id = anno[3];
		        		if(!TargetedLocation.contains(j+"\t"+start)) //already exists
	        			{
		        			if((!mention.matches(".*[;,\\{\\}\\(\\)\\[\\]].*")) && !mention.matches("[a-z]{1,4} [0-9]{1,3}"))
		        			{
			        			if(GNormPlus.taxid4gene.contains(id)) // for gene
			        			{
			        				id="*"+id;
				        		}
			        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tStrain\t"+id);
								String mentions_tmp=mention.toLowerCase();
	    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
	    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
	    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
			        			IDset.add(id);
			        			for(int s=Integer.parseInt(start);s<Integer.parseInt(last);s++)
								{
									TargetedLocation.add(j+"\t"+s);
								}
		        			}
		        		}
					}
				}
			}
			
			HashMap<String, String> OtherNames = new HashMap<String, String>();
			for(String men : Mention2ID_lc.keySet())
			{
				String men_id= Mention2ID_lc.get(men);
				if(GNormPlus.PmidLF2Abb_lc_hash.containsKey(Pmid+"\t"+men))
				{
					String Abb = GNormPlus.PmidLF2Abb_lc_hash.get(Pmid+"\t"+men);
					// Abbreviation
					if(OtherNames.containsKey(men_id))
					{
						OtherNames.put(men_id, OtherNames.get(men_id)+"|"+Abb);
					}
					else
					{
						OtherNames.put(men_id,Abb);
					}
				}
				String men_nospace=men.replaceAll(" ", "");
				// no space
				if(OtherNames.containsKey(men_id))
				{
					OtherNames.put(men_id, OtherNames.get(men_id)+"|"+men_nospace);
				}
				else
				{
					OtherNames.put(men_id,men_nospace);
				}
			}
			PrefixTree PT_Others = new PrefixTree();
			PT_Others.Hash2Tree(OtherNames);
			
			/** 
			 *
			 * Others: 
			 * 1) Abbreviation 
			 * 2) no space
			 * 
			 * */
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>i && 
					GNormPlus.BioCDocobj.PassageContexts.get(i).size()>j &&  
					GNormPlus.BioCDocobj.Annotations.size()>i &&  
					GNormPlus.BioCDocobj.Annotations.get(i).size()>j 
					)
				{
					String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
					ArrayList<String> locations_Abb = PT_Others.SearchMentionLocation(PassageContext,"Species"); /** PT_Abb*/
					for (int k = 0 ; k < locations_Abb.size() ; k++)
					{
						String anno[]=locations_Abb.get(k).split("\t");
						String start= anno[0];
		        		String last= anno[1];
		        		String mention = anno[2];
		        		String id = anno[3];
		        		if(!TargetedLocation.contains(j+"\t"+start)) //already exists
	        			{
		        			if(GNormPlus.taxid4gene.contains(id)) // for gene
		        			{
		        				id="*"+id;
			        		}
		        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tSpecies\t"+id);
							String mentions_tmp=mention.toLowerCase();
    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
		        			Mention2ID_lc.put(mention.toLowerCase(), id);
		        			IDset.add(id);
		        			for(int s=Integer.parseInt(start);s<Integer.parseInt(last);s++)
							{
								TargetedLocation.add(j+"\t"+s);
							}
		        		}
					}
				}
			}
			
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>i && GNormPlus.BioCDocobj.PassageContexts.get(i).size()>j && GNormPlus.BioCDocobj.Annotations.size()>i && GNormPlus.BioCDocobj.Annotations.get(i).size()>j)
				{
					ArrayList <Integer> remove_anno = new ArrayList <Integer>();
					for (int a = 0; a < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); a++) /** Annotations : a */
					{
						String SpAnno[]=GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(a).split("\t");
						String start= SpAnno[0];
		        		String last= SpAnno[1];
		        		String mention = SpAnno[2];
		        		String type = SpAnno[3];
		        		
		        		if(type.matches("Gene|FamilyName|GENERIF"))
		        		{
		        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(a,start+"\t"+last+"\t"+mention+"\t"+type);
			        	}
		        		else if(type.matches("Species|Genus|Strain|Cell") && SpAnno.length==5)
		        		{
		        			//System.out.println(GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(a));
			        		/** Abbreviation solution */
			    			if(GNormPlus.PmidAbb2LF_lc_hash.containsKey(Pmid+"\t"+mention.toLowerCase()) && Mention2ID_lc.containsKey(GNormPlus.PmidAbb2LF_lc_hash.containsKey(Pmid+"\t"+mention.toLowerCase())))
							{
								String LF_lc=GNormPlus.PmidAbb2LF_lc_hash.get(Pmid+"\t"+mention.toLowerCase());
								if(Mention2ID_lc.containsKey(LF_lc))
								{
									String LF_ID=Mention2ID_lc.get(LF_lc);
									GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(a, start+"\t"+last+"\t"+mention+"\t"+type+"\t"+LF_ID);
									String mentions_tmp=mention.toLowerCase();
		    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
		    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
		    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
								}
							}
			    			else if (SpAnno.length>4) 
			    			{
			    				String id = SpAnno[4];
			    				String id_split[]=id.split(";");
		    					if(id_split.length>=2)
		    					{
		    						/** Smallest set of tax ids */
				    				boolean found=false;
		    						for(int x=0;x<IDset.size();x++)
				    				{
				    					String id_tmp= IDset.get(x);
				    					for(int y=0;y<id_split.length;y++) // if any other id is a component of the target id
				    					{
				    						if(id_split[y].equals(id_tmp))
				    						{
				    							found=true;
				    						}
				    					}
				    					if(found == true)
				    					{
				    						GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(a, start+"\t"+last+"\t"+mention+"\t"+type+"\t"+id_tmp);
											String mentions_tmp=mention.toLowerCase();
				    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
				    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
				    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
				    						x=1000000;
				    					}
			    					}
		    						
		    						/** smallest tax id number */
		    						if(found == false)
			    					{
		    							int min=10000000;
				    					String min_id="";
			    						for(int y=0;y<id_split.length;y++) // if any other id is a component of the target id
				    					{
			    							String id_tmp = id_split[y];
			    							String patt="^\\**([0-9]+)";
			    							Pattern ptmp = Pattern.compile(patt);
			    							Matcher mtmp = ptmp.matcher(id_tmp);
			    							if(mtmp.find())
			    							{
			    								id_tmp = mtmp.group(1);
			    							}
			    								
			    							if(y==0)
			    							{
			    								min_id=id_split[y];
			    								min=Integer.parseInt(id_tmp);
			    							}
			    							else if(Integer.parseInt(id_tmp)<min)
			    							{
			    								min=Integer.parseInt(id_tmp);
			    								min_id=id_tmp;
			    							}
				    					}
			    						if(GNormPlus.taxid4gene.contains(min_id)) // for gene
			    						{
			    							min_id="*"+min_id;
			    						}
					        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(a,start+"\t"+last+"\t"+mention+"\tSpecies\t"+min_id);
										String mentions_tmp=mention.toLowerCase();
			    						mentions_tmp=mentions_tmp.replaceAll("[\\W\\-\\_]","");
			    						mentions_tmp=mentions_tmp.replaceAll("[0-9]","0");
			    						GNormPlus.Filtering_hash.put(mentions_tmp,"");
					        		}
			    				}
			    			}
		        		}
		        		else //disease, and other concepts
		        		{
		        			remove_anno.add(a);
		        		}
					}
					
					Collections.sort(remove_anno);
					for (int counter = remove_anno.size()-1; counter >= 0 ; counter--) 
					{ 		      
						int ai=remove_anno.get(counter);
						//System.out.println("\n"+ai+"\t"+GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(ai));
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).remove(ai);
					}
				}
			}	
		}
		GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,true); //save in BioC file
	}
	public void SpeciesAssignment(String Filename,String FilenameBioC) throws IOException, XMLStreamException
	{
		GNormPlus.BioCDocobj.Annotations = new ArrayList();
		GNormPlus.BioCDocobj.BioCReaderWithAnnotation(Filename);
		
		BreakIterator iterator = BreakIterator.getSentenceInstance(Locale.US);	
		for (int i = 0; i < GNormPlus.BioCDocobj.Annotations.size(); i++) /** PMIDs : i */
		{
			HashMap<String, String> PrefixIDTarget_hash = new HashMap<String, String>();
			PrefixIDTarget_hash.put("9606", "h");
			PrefixIDTarget_hash.put("10090", "m");
			PrefixIDTarget_hash.put("10116", "r");
			PrefixIDTarget_hash.put("4932", "y");
			PrefixIDTarget_hash.put("7227", "d");
			PrefixIDTarget_hash.put("7955", "z|zf|Zf|dr|Dr");
			PrefixIDTarget_hash.put("3702", "at|At");
			
			HashMap<String, Double> SP2Num_hash = new HashMap<String, Double>();
			for (int j = 0; j < GNormPlus.BioCDocobj.Annotations.get(i).size(); j++) /** Paragraphs : j */
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					if(anno.length==5) //Species
	        		{
						String patt="^\\**([0-9]+)$";
						Pattern ptmp = Pattern.compile(patt);
						Matcher mtmp = ptmp.matcher(anno[4]);
						if(mtmp.find())
						{
							String id = mtmp.group(1);
							
							if(!PrefixIDTarget_hash.containsKey(id))
							{
								PrefixIDTarget_hash.put(id,GNormPlus.PrefixID_hash.get(id)); // taxid -> prefix
							}
							if(j == 0)//title
			        		{
			        			if(SP2Num_hash.containsKey(id))
			        			{
			        				SP2Num_hash.put(id, SP2Num_hash.get(id)+2);
			        			}
			        			else
			        			{
			        				if(GNormPlus.TaxFreq_hash.containsKey(id))
			        				{
			        					SP2Num_hash.put(id, GNormPlus.TaxFreq_hash.get(id)+2);
			        				}
			        				else
			        				{
			        					SP2Num_hash.put(id, 2.0);
			        				}
			        			}
			        			// Virus -> Human (not to double weight human to virus)
			        			/*if(GNormPlus.SP_Virus2Human_hash.containsKey(id))
		        				{
			        				if(SP2Num_hash.containsKey("9606"))
				        			{
				        				SP2Num_hash.put("9606", SP2Num_hash.get("9606")+2);
				        			}
				        			else
				        			{
				        				SP2Num_hash.put("9606", 2 + GNormPlus.TaxFreq_hash.get("9606")+1);
				        			}
		        				}*/
			        		}
			        		else
			        		{
			        			if(SP2Num_hash.containsKey(id))
			        			{
			        				SP2Num_hash.put(id, SP2Num_hash.get(id)+1);
			        			}
			        			else
			        			{
			        				if(GNormPlus.TaxFreq_hash.containsKey(id))
			        				{
			        					SP2Num_hash.put(id, 1 + GNormPlus.TaxFreq_hash.get(id));
			        				}
			        				else
			        				{
			        					SP2Num_hash.put(id, 1.0);
			        				}
			        			}
			        			// Virus -> Human
			        			/*if(GNormPlus.SP_Virus2Human_hash.containsKey(id))
		        				{
			        				if(SP2Num_hash.containsKey("9606"))
				        			{
				        				SP2Num_hash.put("9606", SP2Num_hash.get("9606")+1);
				        			}
				        			else
				        			{
				        				SP2Num_hash.put("9606", GNormPlus.TaxFreq_hash.get("9606")+1);
				        			}
		        				}*/
			        		}
						}
	        		}
				}
			}
			String MajorSP="9606";
			double MaxSP=0;
			for(String tid : SP2Num_hash.keySet())
			{
				if(SP2Num_hash.get(tid)>MaxSP)
				{
					MajorSP=tid;
					MaxSP=SP2Num_hash.get(tid);
				}
			}
			
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageContexts.get(i).size(); j++) /** Paragraphs : j */
			{
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
				//int PassageOffset = GNormPlus.BioCDocobj.PassageOffsets.get(i).get(j); // Passage offset
				iterator.setText(PassageContext);
				ArrayList<Integer> Sentence_offsets = new ArrayList<Integer>();
				int Sent_start = iterator.first();
				for (int Sent_last = iterator.next(); Sent_last != BreakIterator.DONE; Sent_start = Sent_last, Sent_last = iterator.next()) 
				{
					Sentence_offsets.add(Sent_start);
				}
				
				HashMap<Integer,String> Annotations_Gene_hash = new HashMap<Integer,String>();
				ArrayList<String> Annotations_Species = new ArrayList<String>();
				if(GNormPlus.BioCDocobj.Annotations.get(i).size()>j)
				{
					for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
					{
						String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
						if(anno.length==5) //Species
		        		{
							Annotations_Species.add(GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k));
		        		}
		        		else //Gene : if(anno.length==3)
		        		{
		        			//String mention = PassageContext.substring(Integer.parseInt(anno[0]), Integer.parseInt(anno[1]));
		        			Annotations_Gene_hash.put(k,GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k)); // k -> Gene Annotation
		        		}
					}
	
					//Gene --> Species Inference (PMID:28777492)
					HashMap<String,HashMap<Integer,String>> mention2Location2Species_hash = new HashMap<String,HashMap<Integer,String>>(); 
					HashMap<Integer,String> Location2Species_hash = new HashMap<Integer,String>(); 
					for (int k : Annotations_Gene_hash.keySet()) // k is the index of GNormPlus.BioCDocobj.Annotations.get(i).get(j) 
	    			{
						boolean SPfound = false;
						String anno[] = Annotations_Gene_hash.get(k).split("\t");
	    				int G_Start= Integer.parseInt(anno[0]);
		        		int G_Last= Integer.parseInt(anno[1]);
		        		String G_mentions = anno[2];
		        		/**
		        		 *  2. Co-occurring word
		        		 *  boundary : 
		        		 *  Sentence Start: Sentence_offsets.get(Target_Sentence)
		        		 *  Sentence Last: Sentence_offsets.get(Target_Sentence+1)
		        		 */
		        		//Find the target sentence
		        		int Target_Sentence=0;
		        		if(SPfound == false) // 1. left : Closed to start of the gene mention 
		        		{
		        			for(int s=0;s<Sentence_offsets.size();s++)
		        		
			        		{
			        			int Sentence_last=1000000;
			        			if(s<Sentence_offsets.size()-1)
			        			{
			        				Sentence_last=Sentence_offsets.get(s+1);
			        			}
			        			if(G_Start<Sentence_last)
			        			{
			        				Target_Sentence=s;
			        				break;
			        			}
			        		}
		        		}
		        		int Sentence_Start = Sentence_offsets.get(Target_Sentence);
		        		int Sentence_Last = 1000000;
		        		if(Sentence_offsets.size() > Target_Sentence+1){ Sentence_Last = Sentence_offsets.get(Target_Sentence+1); }
		        		if(SPfound == false) // 1. left : Closed to start of the gene mention 
		        		{
		        			int closet_Sp_Start=0;
		        			for(int sp=0;sp<Annotations_Species.size();sp++) // Find the closet species
			        		{
			        			String AnnoSp[]=Annotations_Species.get(sp).split("\t");
			        			int Sp_Start = Integer.parseInt(AnnoSp[0]);
				        		String patt="^\\**([0-9]+)$";
								Pattern ptmp = Pattern.compile(patt);
								Matcher mtmp = ptmp.matcher(AnnoSp[4]);
								if(mtmp.find())
								{
									String taxid = mtmp.group(1);
									Location2Species_hash.put(Sp_Start,taxid);
				        			if(Sp_Start <= G_Start && Sp_Start >= Sentence_Start && Sp_Start >closet_Sp_Start)
					        		{
					        			closet_Sp_Start=Sp_Start;
					        			Location2Species_hash.put(Integer.parseInt(anno[0]), taxid);

					        			if(mention2Location2Species_hash.containsKey(G_mentions.toLowerCase()))
					        			{
					        				mention2Location2Species_hash.get(G_mentions.toLowerCase()).put(Integer.parseInt(anno[0]), taxid);
					        			}
					        			else
					        			{
					        				mention2Location2Species_hash.put(G_mentions.toLowerCase(),Location2Species_hash);
					        			}
					        			
					        			SPfound=true;
					        		}
								}
				        	}
			        	}
		        		if(SPfound == false) // 2. right : Closed to last of the gene mention
		        		{
		        			int closet_Sp_Last=1000000;
		        			for(int sp=0;sp<Annotations_Species.size();sp++) // Find the closet species
			        		{
			        			String AnnoSp[]=Annotations_Species.get(sp).split("\t");
			        			int Sp_Last = Integer.parseInt(AnnoSp[1]);
				        		String patt="^\\**([0-9]+)$";
								Pattern ptmp = Pattern.compile(patt);
								Matcher mtmp = ptmp.matcher(AnnoSp[4]);
								if(mtmp.find())
								{
									String taxid = mtmp.group(1);
				        			if(Sp_Last >= G_Last && Sp_Last <= Sentence_Last && Sp_Last < closet_Sp_Last)
					        		{
					        			closet_Sp_Last=Sp_Last;
					        			Location2Species_hash.put(Integer.parseInt(anno[0]), taxid);
					        			
					        			if(mention2Location2Species_hash.containsKey(G_mentions.toLowerCase()))
					        			{
					        				mention2Location2Species_hash.get(G_mentions.toLowerCase()).put(Integer.parseInt(anno[0]), taxid);
					        			}
					        			else
					        			{
					        				mention2Location2Species_hash.put(G_mentions.toLowerCase(),Location2Species_hash);
					        			}
					        			
					        			SPfound=true;
					        		}
								}
				        	}
		        		}
	    			}
					
					for (int k : Annotations_Gene_hash.keySet()) // k is the index of GNormPlus.BioCDocobj.Annotations.get(i).get(j) 
	    			{
						String anno[] = Annotations_Gene_hash.get(k).split("\t");
	    				int G_Start= Integer.parseInt(anno[0]);
		        		int G_Last= Integer.parseInt(anno[1]);
		        		String G_mentions = anno[2];
		        		String G_type = anno[3];
		        		String G_mention_list[]=G_mentions.split("\\|");
		        		String G_mention=G_mention_list[0]; // only use the first term to detect species ; should be updated after SimConcept
		        		
		        		/** 1. prefix */
		        		boolean SPfound = false;
		        		for(String taxid: PrefixIDTarget_hash.keySet())
		        		{
							if(GNormPlus.GeneWithoutSPPrefix_hash.containsKey(G_mention.toLowerCase()))
							{
								//special case, and no need for prefix - SA
							}
							else
							{
								Pattern ptmp = Pattern.compile("^("+PrefixIDTarget_hash.get(taxid)+")([A-Z].*)$");
								Matcher mtmp = ptmp.matcher(G_mention);
								if(mtmp.find())
								{
									String MentionWoPrefix=mtmp.group(2);
									GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, anno[0]+"\t"+anno[1]+"\t"+anno[2]+"|"+MentionWoPrefix+"\t"+anno[3]+"\tPrefix:"+taxid);
									SPfound=true;
									break;
								}
							}
		        		}
		        		
		        		/**
		        		 *  2. Co-occurring word
		        		 *  boundary : 
		        		 *  Sentence Start: Sentence_offsets.get(Target_Sentence)
		        		 *  Sentence Last: Sentence_offsets.get(Target_Sentence+1)
		        		 */
		        		//Find the target sentence
		        		int Target_Sentence=0;
		        		if(SPfound == false) // 1. left : Closed to start of the gene mention 
		        		{
		        			for(int s=0;s<Sentence_offsets.size();s++)
		        		
			        		{
			        			int Sentence_last=1000000;
			        			if(s<Sentence_offsets.size()-1)
			        			{
			        				Sentence_last=Sentence_offsets.get(s+1);
			        			}
			        			if(G_Start<Sentence_last)
			        			{
			        				Target_Sentence=s;
			        				break;
			        			}
			        		}
		        		}
		        		int Sentence_Start = Sentence_offsets.get(Target_Sentence);
		        		int Sentence_Last = 1000000;
		        		if(Sentence_offsets.size() > Target_Sentence+1){ Sentence_Last = Sentence_offsets.get(Target_Sentence+1); }
		        		if(SPfound == false) // 1. left : Closed to start of the gene mention 
		        		{
		        			int closet_Sp_Start=0;
		        			for(int sp=0;sp<Annotations_Species.size();sp++) // Find the closet species
			        		{
			        			String AnnoSp[]=Annotations_Species.get(sp).split("\t");
			        			int Sp_Start = Integer.parseInt(AnnoSp[0]);
				        		String patt="^\\**([0-9]+)$";
								Pattern ptmp = Pattern.compile(patt);
								Matcher mtmp = ptmp.matcher(AnnoSp[4]);
								if(mtmp.find())
								{
									String taxid = mtmp.group(1);
				        			if(Sp_Start <= G_Start && Sp_Start >= Sentence_Start && Sp_Start >closet_Sp_Start)
					        		{
					        			closet_Sp_Start=Sp_Start;
										if(GNormPlus.SP_Virus2Human_hash.containsKey(taxid))
				        				{
					        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tLeft:"+taxid+"&9606");
				        				}
					        			else
					        			{
					        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tLeft:"+taxid);
					        			}
					        			SPfound=true;
					        		}
								}
				        	}
			        	}
		        		if(SPfound == false) // 2. right : Closed to last of the gene mention
		        		{
		        			int closet_Sp_Last=1000000;
		        			for(int sp=0;sp<Annotations_Species.size();sp++) // Find the closet species
			        		{
			        			String AnnoSp[]=Annotations_Species.get(sp).split("\t");
			        			int Sp_Last = Integer.parseInt(AnnoSp[1]);
				        		String patt="^\\**([0-9]+)$";
								Pattern ptmp = Pattern.compile(patt);
								Matcher mtmp = ptmp.matcher(AnnoSp[4]);
								if(mtmp.find())
								{
									String taxid = mtmp.group(1);
				        			if(Sp_Last >= G_Last && Sp_Last <= Sentence_Last && Sp_Last < closet_Sp_Last)
					        		{
					        			closet_Sp_Last=Sp_Last;
										if(GNormPlus.SP_Virus2Human_hash.containsKey(taxid))
				        				{
					        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tRight:"+taxid+"&9606");
				        				}
					        			else
					        			{
					        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tRight:"+taxid);
					        			}
					        			SPfound=true;
					        		}
								}
				        	}
		        		}
		        		
		    			/** 3. Focus species */
		        		if(SPfound == false) // 2. right : Closed to last of the gene mention
		        		{
		        			// 1. only the mentions appeared earlier are inferred
		        			//
		        			if(mention2Location2Species_hash.containsKey(G_mentions.toLowerCase()))
		        			{
		        				int closed_loca=0;
		        				for (int loca_start : mention2Location2Species_hash.get(G_mentions.toLowerCase()).keySet())
			        			{
		        					if(loca_start<G_Start)
									{
										if(loca_start>closed_loca)
										{
											closed_loca=loca_start;
										}
									}
			        			}
		        				if(closed_loca>0)
								{
		        					if(GNormPlus.SP_Virus2Human_hash.containsKey(Location2Species_hash.get(closed_loca)))
			        				{
				        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tFocus:"+Location2Species_hash.get(closed_loca)+"&9606");
			        				}
				        			else
				        			{
					        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tFocus:"+Location2Species_hash.get(closed_loca));
				        			}
								}
								else
								{
									if(GNormPlus.SP_Virus2Human_hash.containsKey(MajorSP))
			        				{
				        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tFocus:"+MajorSP+"&9606");
			        				}
				        			else
				        			{
					        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tFocus:"+MajorSP);
				        			}
								}
		        			}
		        			else
							{
		        				if(GNormPlus.SP_Virus2Human_hash.containsKey(MajorSP))
		        				{
			        				GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tFocus:"+MajorSP+"&9606");
		        				}
			        			else
			        			{
				        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tFocus:"+MajorSP);
			        			}
							}
		        		}
	    			}
				}
			}
		}
		GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,true);
	}
	public void SpeciesAssignment(String Filename,String FilenameBioC,String FocusSpecies) throws IOException, XMLStreamException
	{
		for (int i = 0; i < GNormPlus.BioCDocobj.Annotations.size(); i++) /** PMIDs : i */
		{
			for (int j = 0; j < GNormPlus.BioCDocobj.Annotations.get(i).size(); j++) /** Paragraphs : j */
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					if(anno.length==5) //Species
	        		{
						String id=anno[4].replaceAll("\\*", "");
						GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, anno[0]+"\t"+anno[1]+"\t"+anno[2]+"\t"+anno[3]+"\t"+id);
	        		}
	        		else //Gene : if(anno.length==3)
	        		{
	        			/** 1. prefix */
		        		boolean SPfound = false;
		        		if(GNormPlus.GeneWithoutSPPrefix_hash.containsKey(anno[2].toLowerCase()))
						{
							//special case, and no need for prefix - SA
						}
						else
						{
							Pattern ptmp = Pattern.compile("^("+GNormPlus.PrefixID_hash.get(FocusSpecies)+")([A-Z].*)$");
							Matcher mtmp = ptmp.matcher(anno[2]);
							if(mtmp.find())
							{
								String MentionWoPrefix=mtmp.group(2);
								GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, anno[0]+"\t"+anno[1]+"\t"+anno[2]+"|"+MentionWoPrefix+"\t"+anno[3]+"\tPrefix:"+FocusSpecies);
								SPfound=true;
							}
						}
		        		if(SPfound == false)
		        		{
		        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k,  GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k)+"\tFocus:"+FocusSpecies);
		        		}
	        		}
				}
			}
		}
		GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,true);
	}
}