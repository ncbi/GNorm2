/**
 * Project: GNormPlus
 * Function: Gene Name Recognition
 */

package GNormPluslib;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.xml.stream.XMLStreamException;

import org.tartarus.snowball.SnowballStemmer;
import org.tartarus.snowball.ext.englishStemmer;

import GNormPluslib.GNormPlus;
import GNormPluslib.BioCDoc;

public class GNR 
{
	/*
	 * Read BioC files
	 */
	public void Ab3P(String Filename,String FilenameAbb,String TrainTest) throws XMLStreamException,IOException
	{
		/** Abbreviation*/
		//BioC -> Abb input
		String line="";
		BufferedWriter FileAbb = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenameAbb), "UTF-8"));
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++)
		{
			String Pmid = GNormPlus.BioCDocobj.PMIDs.get(i);
			String Context="";
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++)
			{
				String PassageContext=GNormPlus.BioCDocobj.PassageContexts.get(i).get(j);
				if(PassageContext.matches(".*\\([^\\(\\)]+,[^\\(\\)]+\\).*"))
				{
					PassageContext=PassageContext.replaceAll("\\([^\\(\\)]+,[^\\(\\)]+\\)", "");
				}
				if(PassageContext.contains("\\("))
				{
					Context = Context+PassageContext+" ";
				}
			}
			FileAbb.write(Pmid+"\n"+Context+"\n\n");
		}
		FileAbb.close();
		//Abb
		File f = new File(FilenameAbb+".out");
        BufferedWriter fr = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(f), "UTF-8"));
		Runtime runtime = Runtime.getRuntime();
	    String cmd ="./Ab3P "+FilenameAbb+".Abb "+FilenameAbb+".out";
	    
	    String OS=System.getProperty("os.name").toLowerCase();
	    if(OS.contains("windows"))
	    {
	    	cmd ="java -jar bioadi.jar "+FilenameAbb;
	    }
	    else //if(OS.contains("nux")||OS.contains("nix"))
	    {
	    	cmd ="./Ab3P "+FilenameAbb+" "+FilenameAbb+".out";
	    	//cmd ="java -jar bioadi.jar "+FilenameAbb+" > "+FilenameAbb+".out";
	    }
	    
    	Process process = runtime.exec(cmd);
    	InputStream is = process.getInputStream();
    	InputStreamReader isr = new InputStreamReader(is, "UTF-8");
    	BufferedReader br = new BufferedReader(isr);
    	line="";
	    while ( (line = br.readLine()) != null) 
	    {
	    	fr.write(line);
	    	fr.newLine();
	        fr.flush();
	    }
	    is.close();
	    isr.close();
	    br.close();
	    fr.close();
	    //Abb output -> Hash
		BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameAbb+".out"), "UTF-8"));
		line="";
		String pmid="";
		while ((line = inputfile.readLine()) != null)  
		{
			String patt="^  (.+)\\|(.+)\\|([0-9\\.]+)$";
			Pattern ptmp = Pattern.compile(patt);
			Matcher mtmp = ptmp.matcher(line);
			if(line.matches("^[0-9]+$"))
			{
				pmid=line;
			}
			if(mtmp.find())
			{
				String SF = mtmp.group(1);
				String LF = mtmp.group(2);
				double weight=  Double.parseDouble(mtmp.group(3));
				GNormPlus.Pmid2Abb_hash.put(pmid+"\t"+SF, "Abb:SF");
				GNormPlus.Pmid2Abb_hash.put(pmid+"\t"+LF, "Abb:LF");
				GNormPlus.PmidLF2Abb_lc_hash.put(pmid+"\t"+LF.toLowerCase(), SF.toLowerCase());
				GNormPlus.PmidAbb2LF_lc_hash.put(pmid+"\t"+SF.toLowerCase(), LF.toLowerCase());
				GNormPlus.PmidAbb2LF_hash.put(pmid+"\t"+SF, LF);
				if(weight >= 0.9)
				{
					GNormPlus.PmidLF2Abb_hash.put(pmid+"\t"+LF, SF);
				}
			}
		}
		inputfile.close();	
	}
	
	public void LoadInputFile(String Filename,String FilenameAbb,String TrainTest) throws XMLStreamException,IOException
	{
		/** Read BioC file */
		//if(TrainTest.equals("Train"))
		//{
			GNormPlus.BioCDocobj.BioCReaderWithAnnotation(Filename);
		//}
		//else
		//{
		//	GNormPlus.BioCDocobj.BioCReader(Filename);
		//}
		
		
		/** Abbreviation*/
		//BioC -> Abb input
		String line="";
		BufferedWriter FileAbb = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenameAbb), "UTF-8"));
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++)
		{
			String Pmid = GNormPlus.BioCDocobj.PMIDs.get(i);
			String Context="Text:";
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++)
			{
				String PassageContext=GNormPlus.BioCDocobj.PassageContexts.get(i).get(j);
				if(PassageContext.matches(".*\\([^\\(\\)]+,[^\\(\\)]+\\).*"))
				{
					PassageContext=PassageContext.replaceAll("\\([^\\(\\)]+,[^\\(\\)]+\\)", "");
				}
				if(PassageContext.contains("("))
				{
					Context = Context+PassageContext+" ";
				}
			}
			FileAbb.write(Pmid+"\n"+Context+"\n\n");
		}
		FileAbb.close();
		//Abb
		File f = new File(FilenameAbb+".out");
        BufferedWriter fr = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(f), "UTF-8"));
		Runtime runtime = Runtime.getRuntime();
	    String cmd ="./Ab3P "+FilenameAbb+".Abb "+FilenameAbb+".out";
	    
	    String OS=System.getProperty("os.name").toLowerCase();
	    if(OS.contains("windows"))
	    {
	    	cmd ="java -jar bioadi.jar "+FilenameAbb;
	    }
	    else //if(OS.contains("nux")||OS.contains("nix"))
	    {
	    	cmd ="./Ab3P "+FilenameAbb+" "+FilenameAbb+".out";
	    	//cmd ="java -jar bioadi.jar "+FilenameAbb+" > "+FilenameAbb+".out";
	    }
	    
    	Process process = runtime.exec(cmd);
    	InputStream is = process.getInputStream();
    	InputStreamReader isr = new InputStreamReader(is, "UTF-8");
    	BufferedReader br = new BufferedReader(isr);
    	line="";
	    while ( (line = br.readLine()) != null) 
	    {
	    	fr.write(line);
	    	fr.newLine();
	        fr.flush();
	    }
	    is.close();
	    isr.close();
	    br.close();
	    fr.close();
	    //Abb output -> Hash
		BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameAbb+".out"), "UTF-8"));
		line="";
		String pmid="";
		while ((line = inputfile.readLine()) != null)  
		{
			String patt="^  (.+)\\|(.+)\\|([0-9\\.]+)$";
			Pattern ptmp = Pattern.compile(patt);
			Matcher mtmp = ptmp.matcher(line);
			if(line.matches("^[0-9]+$"))
			{
				pmid=line;
			}
			if(mtmp.find())
			{
				String SF = mtmp.group(1);
				String LF = mtmp.group(2);
				double weight=  Double.parseDouble(mtmp.group(3));
				GNormPlus.Pmid2Abb_hash.put(pmid+"\t"+SF, "Abb:SF");
				GNormPlus.Pmid2Abb_hash.put(pmid+"\t"+LF, "Abb:LF");
				GNormPlus.PmidLF2Abb_lc_hash.put(pmid+"\t"+LF.toLowerCase(), SF.toLowerCase());
				GNormPlus.PmidAbb2LF_lc_hash.put(pmid+"\t"+SF.toLowerCase(), LF.toLowerCase());
				GNormPlus.PmidAbb2LF_hash.put(pmid+"\t"+SF, LF);
				if(weight >= 0.9)
				{
					GNormPlus.PmidLF2Abb_hash.put(pmid+"\t"+LF, SF);
				}
			}
		}
		inputfile.close();	
	}
	
	/*
	 * Feature Extraction
	 */
	public void FeatureExtraction(String FilenameData,String FilenameLoca,String TrainTest) throws XMLStreamException
	{
		try 
		{
			/** output files */ 
			BufferedWriter FileLocation = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenameLoca), "UTF-8")); // .location
			BufferedWriter FileData = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(FilenameData), "UTF-8")); // .data
			//NLP modules
			SnowballStemmer stemmer = new englishStemmer();
			/** PMIDs : i */
			for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++)
			{
				String Pmid = GNormPlus.BioCDocobj.PMIDs.get(i);
				
				/** Paragraphs : j */
				for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++)
				{
					String PassageName= GNormPlus.BioCDocobj.PassageNames.get(i).get(j); // Passage name
					int PassageOffset = GNormPlus.BioCDocobj.PassageOffsets.get(i).get(j); // Passage offset
					String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
					ArrayList<String> Annotation = GNormPlus.BioCDocobj.Annotations.get(i).get(j); // Annotation
					HashMap<Integer, String> CTDGene_hash = new HashMap<Integer, String>();
					HashMap<Integer, String> FamilyName_hash = new HashMap<Integer, String>();
					HashMap<Integer, String> character_hash = new HashMap<Integer, String>();
					HashMap<Integer, String> Abbreviation_hash = new HashMap<Integer, String>();
					String PassageContext_tmp=" "+PassageContext+" ";
					
					/** Abbreviation */
					HashMap<Integer, String> Abb_sortebylength = new HashMap<Integer, String>();
					ArrayList<Integer> length_list = new ArrayList<Integer>();
					int countn=0;
					for (Object key : GNormPlus.Pmid2Abb_hash.keySet())
					{
						String pmid2abb[]=key.toString().split("\t");
						if(Pmid.equals(pmid2abb[0]))
						{
							Abb_sortebylength.put(pmid2abb[1].length()*100+countn, pmid2abb[1]);
							length_list.add(pmid2abb[1].length()*100+countn);
							countn++;
						}
					}
					Collections.sort(length_list);
					for (int l=length_list.size()-1;l>=0;l--)
					{
						String AbbLF = Abb_sortebylength.get(length_list.get(l));
						AbbLF=AbbLF.replaceAll("([^A-Za-z0-9@ ])", "\\\\$1");
						AbbLF=AbbLF.replaceAll(" ", "\\[ \\]\\+");
						Pattern ptmp = Pattern.compile("^(.*[^A-Za-z0-9]+)("+AbbLF+")([^A-Za-z0-9]+.*)$");
						Matcher mtmp = ptmp.matcher(PassageContext_tmp);
						while(mtmp.find())
						{
							String str1=mtmp.group(1);
							String str2=mtmp.group(2);
							String str3=mtmp.group(3);
							for(int m=str1.length();m<=(str1.length()+str2.length());m++)
    						{
								Abbreviation_hash.put((m-1),GNormPlus.Pmid2Abb_hash.get(Pmid+"\t"+Abb_sortebylength.get(length_list.get(l))));
							}
							String men="";
							for(int m=0;m<str2.length();m++){men=men+"@";}
    						PassageContext_tmp=str1+men+str3;
    						mtmp = ptmp.matcher(PassageContext_tmp);
						}
					}
					
					/** PT_CTDGene */
					ArrayList<String> locations = GNormPlus.PT_CTDGene.SearchMentionLocation(PassageContext,"CTDGene");
					for (int k = 0 ; k < locations.size() ; k++)
					{
						String anno[]=locations.get(k).split("\t");
						int start= Integer.parseInt(anno[0]) + PassageOffset;
		        		int last= Integer.parseInt(anno[1]) + PassageOffset;
		        		String mention = anno[2];
		        		String id = anno[3];
		        		
		        		CTDGene_hash.put(start,"CTDGene_B");
		        		CTDGene_hash.put(last,"CTDGene_E");
		        		for(int s=start+1;s<last;s++)
	        			{
		        			CTDGene_hash.put(s,"CTDGene_I");
	        			}
					}
					
					/** FamilyName */
					ArrayList<String> locations_Fname = GNormPlus.PT_FamilyName.SearchMentionLocation(PassageContext,"FamilyName");
					for (int k = 0 ; k < locations_Fname.size() ; k++)
					{
						String anno[]=locations_Fname.get(k).split("\t");
						int start= Integer.parseInt(anno[0]) + PassageOffset;
		        		int last= Integer.parseInt(anno[1]) + PassageOffset;
		        		String mention = anno[2];
		        		String id = anno[3];
		        		
		        		if(!CTDGene_hash.containsKey(start))
		        		{
			        		FamilyName_hash.put(start,"famplex_B");
			        		FamilyName_hash.put(last,"famplex_E");
			        		for(int s=start+1;s<last;s++)
		        			{
			        			FamilyName_hash.put(s,"famplex_I");
		        			}
		        		}
					}
					
					/** Annotations : k 
					 * 0 start
					 * 1 last
					 * 2 mention
					 * 3 type
					 * 4 id
					 */
					for (int k = 0; k < Annotation.size() ; k++)   // k : Annotations
					{
						String anno[]=Annotation.get(k).split("\t");
						int start= Integer.parseInt(anno[0]);
		        		int last= Integer.parseInt(anno[1]);
		        		String type = anno[3];
		        		
		        		character_hash.put(start,type+"_B");
		        		character_hash.put(last,type+"_E");
		        		for(int s=start+1;s<last;s++)
	        			{
	        				character_hash.put(s,type+"_I");
	        			}
					}
					
					String PassageContext_rev=PassageContext;
					PassageContext_rev = PassageContext_rev.replaceAll("([A-Z][A-Z])([A-Z][0-9][0-9]+[A-Z][\\W\\-\\_])", "$1 $2"); //PTENK289E
					PassageContext_rev = PassageContext_rev.replaceAll("([0-9])([A-Za-z])", "$1 $2");
					PassageContext_rev = PassageContext_rev.replaceAll("([A-Za-z])([0-9])", "$1 $2");
					//PassageContext_rev = PassageContext_rev.replaceAll("([A-Z])([a-z])", "$1 $2");
					//PassageContext_rev = PassageContext_rev.replaceAll("([a-z])([A-Z])", "$1 $2");
					PassageContext_rev = PassageContext_rev.replaceAll("([\\W])", " $1 ");
					PassageContext_rev = PassageContext_rev.replaceAll("[ ]+", " ");
					PassageContext_tmp = PassageContext;
					int Offset=0;
					String tokens[]=PassageContext_rev.split(" ");
	        		for(int p=0;p<tokens.length;p++)
	                {
	        			String WSB="WSB:NoGap"; //no gap to backward
	        			String WSF="WSF:NoGap"; //no gap to forward
	        			while(PassageContext_tmp.substring(0,1).equals(" "))
						{
	        				PassageContext_tmp=PassageContext_tmp.substring(1);
							Offset++;
							WSB="WSB:Gap";
						}
	        			if(PassageContext_tmp.length()>tokens[p].length() && PassageContext_tmp.substring(tokens[p].length(),tokens[p].length()+1).equals(" "))
	        			{
	        				WSF="WSF:Gap";
	        			}
	        			if(p==0)
	        			{
	        				WSB="WSB:1st";
	        			}
	        			else if(p==tokens.length-1)
	        			{
	        				WSF="WSF:last";
	        			}
	        			
	        			if(PassageContext_tmp.substring(0,tokens[p].length()).equals(tokens[p]))
						{
	        				if(tokens[p].length()>0)
							{
	        					/*
	        					 * .loca
	        					 */
	        					int start=Offset;
	        					int last=Offset+tokens[p].length();
	        					String State="";
	        					if(!character_hash.containsKey(start) || !character_hash.containsKey(last)){}
	        					else if(character_hash.get(start).matches(".*B$"))
	        					{
	        						State=character_hash.get(start);
	        					}
	        					else if(character_hash.get(last).matches(".*E$"))
	        					{
	        						State=character_hash.get(last);
	        					}
	        					else if(character_hash.get(start).matches(".*I$"))
	        					{
	        						State=character_hash.get(start);
	        					}
	        					
	        					if((!tokens[p].equals("\t")))
	        					{
	        						FileLocation.write(Pmid+"\t"+PassageName+"\t"+j+"\t"+tokens[p]+"\t"+(Offset+1)+"\t"+(Offset+tokens[p].length())+"\t"+State+"\n");
	        					}
	        					
	        					/*
	        					 * .data
	        					 */
	        					
	        					//Abbreviation
	        					String Abb_State="__nil__";
	        					if(!Abbreviation_hash.containsKey(start) || !Abbreviation_hash.containsKey(last)){Abb_State="__nil__";}
	        					else if(Abbreviation_hash.containsKey(start))
	        					{
	        						Abb_State=Abbreviation_hash.get(start);
	        					}
	        					
	        					//CTDGene
	        					start=PassageOffset+Offset;
	        					last=PassageOffset+Offset+tokens[p].length();
	        					String CTDGene_State="__nil__";
	        					if(!CTDGene_hash.containsKey(start) || !CTDGene_hash.containsKey(last)){CTDGene_State="__nil__";}
	        					else if(CTDGene_hash.get(start).matches(".*B$"))
	        					{
	        						CTDGene_State=CTDGene_hash.get(start);
	        					}
	        					else if(CTDGene_hash.get(last).matches(".*E$"))
	        					{
	        						CTDGene_State=CTDGene_hash.get(last);
	        					}
	        					else if(CTDGene_hash.get(start).matches(".*I$"))
	        					{
	        						CTDGene_State=CTDGene_hash.get(start);
	        					}
	        					
	        					//FamilyName
	        					if(CTDGene_State.equals("__nil__"))
	        					{
		        					start=PassageOffset+Offset;
		        					last=PassageOffset+Offset+tokens[p].length();
		        					if(!FamilyName_hash.containsKey(start) || !FamilyName_hash.containsKey(last)){}
		        					else if(FamilyName_hash.get(start).matches(".*B$"))
		        					{
		        						CTDGene_State=FamilyName_hash.get(start);
		        					}
		        					else if(FamilyName_hash.get(last).matches(".*E$"))
		        					{
		        						CTDGene_State=FamilyName_hash.get(last);
		        					}
		        					else if(FamilyName_hash.get(start).matches(".*I$"))
		        					{
		        						CTDGene_State=FamilyName_hash.get(start);
		        					}
	        					}
	        					
	        					//stemming
								stemmer.setCurrent(tokens[p].toLowerCase());
								stemmer.stem();
								String stem=stemmer.getCurrent();
								
								//Number of Numbers [0-9]
								String Num_num="";
								String tmp=tokens[p];
								tmp=tmp.replaceAll("[^0-9]","");
								if(tmp.length()>3){Num_num="N:4+";}else{Num_num="N:"+ tmp.length();}
								
								//Number of Uppercase [A-Z]
								String Num_Uc="";
								tmp=tokens[p];
								tmp=tmp.replaceAll("[^A-Z]","");
								if(tmp.length()>3){Num_Uc="U:4+";}else{Num_Uc="U:"+ tmp.length();}
								
								//Number of Lowercase [a-z]
								String Num_lc="";
								tmp=tokens[p];
								tmp=tmp.replaceAll("[^a-z]","");
								if(tmp.length()>3){Num_lc="L:4+";}else{Num_lc="L:"+ tmp.length();}
								
								//Number of ALL char
								String Num_All="";
								if(tokens[p].length()>3){Num_All="A:4+";}else{Num_All="A:"+ tokens[p].length();}
								
								//specific character (;:,.->+_)
								String SpecificC="__nil__";
								if(tokens[p].equals(";") || tokens[p].equals(":") || tokens[p].equals(",") || tokens[p].equals(".") || tokens[p].equals("-") || tokens[p].equals(">") || tokens[p].equals("+") || tokens[p].equals("_"))
								{
									SpecificC="-SpecificC1-";
								}
								else if(tokens[p].equals("(") || tokens[p].equals(")"))
								{
									SpecificC="-SpecificC2-";
								}
								else if(tokens[p].equals("{") || tokens[p].equals("}"))
								{
									SpecificC="-SpecificC3-";
								}
								else if(tokens[p].equals("[") || tokens[p].equals("]"))
								{
									SpecificC="-SpecificC4-";
								}
								else if(tokens[p].equals("\\") || tokens[p].equals("/"))
								{
									SpecificC="-SpecificC5-";
								}
								
								//Chemical Prefix/Suffix
								String ChemPreSuf="__nil__";
								if(tokens[p].matches(".*(yl|ylidyne|oyl|sulfonyl)")){ChemPreSuf="-CHEMinlineSuffix-";}
								else if(tokens[p].matches("(meth|eth|prop|tetracos).*")){ChemPreSuf="-CHEMalkaneStem-";}
								else if(tokens[p].matches("(di|tri|tetra).*")){ChemPreSuf="-CHEMsimpleMultiplier-";}
								else if(tokens[p].matches("(benzen|pyridin|toluen).*")){ChemPreSuf="-CHEMtrivialRing-";}
								else if(tokens[p].matches(".*(one|ol|carboxylic|amide|ate|acid|ium|ylium|ide|uide|iran|olan|inan|pyrid|acrid|amid|keten|formazan|fydrazin)(s|)")){ChemPreSuf="-CHEMsuffix-";}
								
								
								//Mention Type
								String MentionType="__nil__";
								/*
								if($tmp eq "to" && $CTD_result_hash{$count_token-1} eq "CTD_gene" && $CTD_result_hash{$count_token+1} eq "CTD_gene"){$CTD_result_hash{$count_token}="CTD_gene";}
								if($tmp=~/^(or|and|,)$/ && $CTD_result_hash{$count_token-1} eq "CTD_gene" && $CTD_result_hash{$count_token+1} eq "CTD_gene"){$MentionType="-Type_GeneConjunction-";}
								elsif($tmp=~/^(or|and|,)$/ && $last_token=~/^(or|and|,)$/ && $CTD_result_hash{$count_token-2} eq "CTD_gene" && $CTD_result_hash{$count_token+1} eq "CTD_gene"){$MentionType="-Type_GeneConjunction-";}
								elsif($tmp=~/^(or|and|,)$/ && $next_token=~/^(or|and|,)$/ && $CTD_result_hash{$count_token-1} eq "CTD_gene" && $CTD_result_hash{$count_token+2} eq "CTD_gene"){$MentionType="-Type_GeneConjunction-";}
								*/
								if(tokens[p].matches("(ytochrome|cytochrome)")){MentionType="-Type_cytochrome-";}
								else if(tokens[p].matches(".*target") ){MentionType="-Type_target-";}
								else if(tokens[p].matches(".*(irradiation|hybrid|fusion|experiment|gst|est|gap|antigen)") ){MentionType="-Type_ExperimentNoun-";}
								else if(tokens[p].matches(".*(disease|disorder|dystrophy|deficiency|syndrome|dysgenesis|cancer|injury|neoplasm|diabetes|diabete)") ){MentionType="-Type_Disease-";}
								else if(tokens[p].matches(".*(motif|domain|omain|binding|site|region|sequence|frameshift|finger|box).*") ){MentionType="-Type_DomainMotif-";}
								else if(tokens[p].equals("-") && (p<tokens.length-1 && tokens[p+1].matches(".*(motif|domain|omain|binding|site|region|sequence|frameshift|finger|box).*")) ){MentionType="-Type_DomainMotif-";}
								else if(tokens[p].matches("[rmc]") && (p<tokens.length-1 && (tokens[p+1].equals("DNA") || tokens[p+1].equals("RNA")) ) ){MentionType="-Type_DomainMotif-";}
								else if(tokens[p].matches(".*(famil|complex|cluster|proteins|genes|factors|transporter|proteinase|membrane|ligand|enzyme|channels|tors$|ase$|ases$)") ){MentionType="-Type_Family-";}
								else if(tokens[p].toLowerCase().matches("^marker") ){MentionType="-Type_Marker-";}
								else if(tokens[p].equals(".*cell.*") || (p<tokens.length-1 && tokens[p+1].equals("cell") && tokens[p].matches("^(T|B|monocytic|cancer|tumor|myeloma|epithelial|crypt)$") ) ){MentionType="-Type_Cell-";}
								else if(tokens[p].equals(".*chromosome.*") ){MentionType="-Type_Chromosome-";}
								else if(tokens[p].matches("[pq]") && ( (p<tokens.length-1 && tokens[p+1].matches("^[0-9]+$")) || (p>0 && tokens[p-1].matches("^[0-9]+$")) ) ){MentionType="-Type_ChromosomeStrain-";}
								else if(tokens[p].matches(".*(related|regulated|associated|correlated|reactive).*")){MentionType="-Type_relation-";}
								else if(tokens[p].toLowerCase().matches(".*(polymorphism|mutation|deletion|insertion|duplication|genotype|genotypes).*") ){MentionType="-Type_VariationTerms-";}
								else if(tokens[p].matches(".*(oxidase|transferase|transferases|kinase|kinese|subunit|unit|receptor|adrenoceptor|transporter|regulator|transcription|antigen|protein|gene|factor|member|molecule|channel|deaminase|spectrin).*") ){MentionType="-Type_suffix-";}
								else if(tokens[p].matches("[\\(\\-\\_]") && (p<tokens.length-1 && tokens[p+1].toLowerCase().matches(".*(alpha|beta|gamma|delta|theta|kappa|zeta|sigma|omega|i|ii|iii|iv|v|vi|[abcdefgyr])")) ){MentionType="-Type_strain-";}
								else if(tokens[p].matches("(alpha|beta|gamma|delta|theta|kappa|zeta|sigma|omega|i|ii|iii|iv|v|vi|[abcdefgyr])") ){MentionType="-Type_strain-";}
								
								//Protein symbols
								String ProteinSym="__nil__";
								if(tokens[p].matches(".*(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|glycine|aspartate|methionine|threonine|histidine|aspartic|asparticacid|arginine|asparagine|tryptophan|proline|phenylalanine|cysteine|serine|glutamate|tyrosine|stop|frameshift).*")){ProteinSym="-ProteinSymFull-";}
								else if(tokens[p].matches("(cys|ile|ser|gln|met|asn|pro|lys|asp|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr|fs|fsx)")){ProteinSym="-ProteinSymTri-";}
								else if(tokens[p].matches("[CISQMNPKDTFAGHLRWVEYX]")){ProteinSym="-ProteinSymChar-";}
								
								//prefix
								String prefix="";
								tmp=tokens[p];
								if(tmp.length()>=1){ prefix=tmp.substring(0, 1);}else{prefix="__nil__";}
								if(tmp.length()>=2){ prefix=prefix+" "+tmp.substring(0, 2);}else{prefix=prefix+" __nil__";}
								if(tmp.length()>=3){ prefix=prefix+" "+tmp.substring(0, 3);}else{prefix=prefix+" __nil__";}
								if(tmp.length()>=4){ prefix=prefix+" "+tmp.substring(0, 4);}else{prefix=prefix+" __nil__";}
								if(tmp.length()>=5){ prefix=prefix+" "+tmp.substring(0, 5);}else{prefix=prefix+" __nil__";}
								
								
								//suffix
								String suffix="";
								tmp=tokens[p];
								if(tmp.length()>=1){ suffix=tmp.substring(tmp.length()-1, tmp.length());}else{suffix="__nil__";}
								if(tmp.length()>=2){ suffix=suffix+" "+tmp.substring(tmp.length()-2, tmp.length());}else{suffix=suffix+" __nil__";}
								if(tmp.length()>=3){ suffix=suffix+" "+tmp.substring(tmp.length()-3, tmp.length());}else{suffix=suffix+" __nil__";}
								if(tmp.length()>=4){ suffix=suffix+" "+tmp.substring(tmp.length()-4, tmp.length());}else{suffix=suffix+" __nil__";}
								if(tmp.length()>=5){ suffix=suffix+" "+tmp.substring(tmp.length()-5, tmp.length());}else{suffix=suffix+" __nil__";}
								
								if(State.equals(""))
								{
									State="O";
								}
								
								if((!tokens[p].equals("\t")))
	        					{
									if(TrainTest.equals("Train"))
		        					{
										FileData.write(tokens[p]+" "+stem+" "+WSB+" "+WSF+" "+Num_num+" "+Num_Uc+" "+Num_lc+" "+Num_All+" "+SpecificC+" "+ChemPreSuf+" "+MentionType+" "+ProteinSym+" "+prefix+" "+suffix+" "+CTDGene_State+" "+Abb_State+" "+State+"\n");
									}
									else
									{	
										FileData.write(tokens[p]+" "+stem+" "+WSB+" "+WSF+" "+Num_num+" "+Num_Uc+" "+Num_lc+" "+Num_All+" "+SpecificC+" "+ChemPreSuf+" "+MentionType+" "+ProteinSym+" "+prefix+" "+suffix+" "+CTDGene_State+" "+Abb_State+"\n");
									}
								}
			        			PassageContext_tmp=PassageContext_tmp.substring(tokens[p].length()); // remove the token for the context
	        					Offset=Offset+tokens[p].length();
							}
						}
	        		}
	        		if(tokens.length>0)
	        		{
		        		FileLocation.write("\n");
						FileData.write("\n");
	        		}
				}
			}
			FileLocation.close();
			FileData.close();
		}
		catch(IOException e1){ System.out.println("[MR]: Input file is not exist.");}
	}
	/*
	 * Testing by CRF++
	 */
	public void CRF_test(String model, String FilenameData, String FilenameOutput) throws IOException 
	{
		File f = new File(FilenameOutput);
        BufferedWriter fr = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(f), "UTF-8"));
		
		Runtime runtime = Runtime.getRuntime();
	    
		String OS=System.getProperty("os.name").toLowerCase();
		
		String cmd="./CRF/crf_test -m "+model+" -o "+FilenameOutput+" "+FilenameData;
	    if(OS.contains("windows"))
	    {
	    	cmd ="CRF/crf_test -m "+model+" -o "+FilenameOutput+" "+FilenameData;
	    }
	    else //if(OS.contains("nux")||OS.contains("nix"))
	    {
	    	cmd ="./CRF/crf_test -m "+model+" -o "+FilenameOutput+" "+FilenameData;
	    }
	    
	    try {
	    	Process process = runtime.exec(cmd);
	    	InputStream is = process.getInputStream();
	    	InputStreamReader isr = new InputStreamReader(is, "UTF-8");
	    	BufferedReader br = new BufferedReader(isr);
	    	String line="";
		    while ( (line = br.readLine()) != null) 
		    {
		    	fr.write(line);
		    	fr.newLine();
		        fr.flush();
		    }
		    is.close();
		    isr.close();
		    br.close();
		    fr.close();
	    }
	    catch (IOException e) {
	    	System.out.println(e);
	    	runtime.exit(0);
	    }
	}
	
	public void CRF_test(String model,String FilenameData,String FilenameOutput,String top3) throws IOException 
	{
		File f = new File(FilenameOutput);
        BufferedWriter fr = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(f), "UTF-8"));
		
		Runtime runtime = Runtime.getRuntime();
	    
		String OS=System.getProperty("os.name").toLowerCase();
		
		String cmd="./CRF/crf_test -n 3 -m "+model+" -o "+FilenameOutput+" "+FilenameData;
	    if(OS.contains("windows"))
	    {
	    	cmd ="CRF/crf_test -n 3 -m "+model+" -o "+FilenameOutput+" "+FilenameData;
	    }
	    else //if(OS.contains("nux")||OS.contains("nix"))
	    {
	    	cmd ="./CRF/crf_test -n 3 -m "+model+" -o "+FilenameOutput+" "+FilenameData;
	    }
	    
		try {
	    	Process process = runtime.exec(cmd);
	    	InputStream is = process.getInputStream();
	    	InputStreamReader isr = new InputStreamReader(is, "UTF-8");
	    	BufferedReader br = new BufferedReader(isr);
	    	String line="";
		    while ( (line = br.readLine()) != null) 
		    {
		    	fr.write(line);
		    	fr.newLine();
		        fr.flush();
		    }
		    is.close();
		    isr.close();
		    br.close();
		    fr.close();
	    }
	    catch (IOException e) {
	    	System.out.println(e);
	    	runtime.exit(0);
	    }
	}
	
	/*
	 * Learning model by CRF++
	 */
	public void CRF_learn(String model, String FilenameData) throws IOException 
	{
		Runtime runtime = Runtime.getRuntime();
	    
	    Process process = null;
	    String line = null;
	    InputStream is = null;
	    InputStreamReader isr = null;
	    BufferedReader br = null;
	    
	    String OS=System.getProperty("os.name").toLowerCase();
		
		String cmd="./CRF/crf_learn -f 3 -c 4.0 CRF/template_UB "+FilenameData+" "+model; 
	    if(OS.contains("windows"))
	    {
	    	cmd ="CRF/crf_learn -f 3 -c 4.0 CRF/template_UB "+FilenameData+" "+model; 
	    }
	    else //if(OS.contains("nux")||OS.contains("nix"))
	    {
	    	cmd ="./CRF/crf_learn -f 3 -c 4.0 CRF/template_UB "+FilenameData+" "+model; 
	    }
	    
	    try {
	    	process = runtime.exec(cmd);
		    is = process.getInputStream();
		    isr = new InputStreamReader(is, "UTF-8");
		    br = new BufferedReader(isr);
		    while ( (line = br.readLine()) != null) 
		    {
		    	System.out.println(line);
		        System.out.flush();
		    }
		    is.close();
		    isr.close();
		    br.close();
	    }
	    catch (IOException e) {
	    	System.out.println(e);
	    	runtime.exit(0);
	    }
	}
	
	public void ReadCRFresult(String Filename,String FilenameLoca,String FilenameOutput,String FilenameBioC) throws XMLStreamException, IOException
	{
		/** load CRF output */
		ArrayList<String> outputArr = new ArrayList<String>();
		BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameOutput), "UTF-8"));
		String line;
		while ((line = inputfile.readLine()) != null)  
		{	
			outputArr.add(line);
		}
		inputfile.close();
		
		/** load location */
		ArrayList<String> locationArr = new ArrayList<String>(); 
		inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameLoca), "UTF-8"));
		while ((line = inputfile.readLine()) != null)  
		{
			locationArr.add(line);
		}
		inputfile.close();
		
		/** output -> mentions */
		String pmid_last="";
		String paragraph_num_last="";
		String pmid="";
		String paragraph="";
		String paragraph_num="";
		Pattern pat_B = Pattern.compile("((FamilyName|DomainMotif|Gene)_[B])$");
		Pattern pat_IE = Pattern.compile("((FamilyName|DomainMotif|Gene)_[IE])$");
		ArrayList<ArrayList<String>> AnnotationInPMID = new ArrayList(); // array of Annotations in the PMIDs
		ArrayList<String> AnnotationInPassage= new ArrayList<String>(); // array of Annotations in the Passage
		GNormPlus.BioCDocobj.Annotations = new ArrayList();
		int countPMID=0;
		int countPassage=0;
		/** outputArr */
		for(int i=0;i<outputArr.size();i++)
		{
			String outputsRow[]=outputArr.get(i).split("\\t");
			String locationRow[]=locationArr.get(i).split("\\t");
			int start=100000;
			int last=0;
			String MentionType="";
			
			if(locationRow.length>3)
			{
				pmid=locationRow[0];
				paragraph=locationRow[1];
				paragraph_num=locationRow[2];
			}
			
			if( (!paragraph_num_last.equals("")) && (!paragraph_num.equals(paragraph_num_last)) )
			{
				AnnotationInPMID.add(AnnotationInPassage);
				AnnotationInPassage = new ArrayList<String>();
				countPassage++;
			}
			if( (!pmid_last.equals("")) && (!pmid.equals(pmid_last)) )
			{
				GNormPlus.BioCDocobj.Annotations.add(AnnotationInPMID);
				AnnotationInPMID = new ArrayList();
				countPMID++;
				countPassage=0;
			}
			
			boolean F = false; //Flag of Finding
			if(locationRow.length>2)
			{
				Matcher mat = pat_B.matcher(outputsRow[outputsRow.length-1]); // last column : Status
				while(mat.find() && locationRow.length==6)
				{
					MentionType=mat.group(2);	
					pmid=locationRow[0];
					paragraph_num=locationRow[2];
					int start_tmp=Integer.parseInt(locationRow[4])-1;
					int last_tmp=Integer.parseInt(locationRow[5]);
					if(start_tmp<start){start=start_tmp;}
					if(last_tmp>last){last=last_tmp;}
					i++;
					F = true;
					if(locationArr.get(i).length()>0)
					{
						outputsRow=outputArr.get(i).split("\\t");
						locationRow=locationArr.get(i).split("\\t");
						mat = pat_IE.matcher(outputsRow[outputsRow.length-1]);
					}
					else
					{
						break;
					}
				}
			}
			
			if(F == true)
			{
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(countPMID).get(countPassage); // Passage context
				String Mention = PassageContext.substring(start, last); 
				String Mention_nospace = Mention.replaceAll("[\\W\\-\\_]", ""); 
				if(Mention.toLowerCase().matches("(figure|tables|fig|tab|exp\\. [0-9]+).*")){}
				else if(Mention.matches("[A-Z][A-Z]s")){}
				else if(Mention.matches(".*\\|.*")){} 
				else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\;\\,\\'\\/\\\\].*")){}
				else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\(].*") && !Mention.matches(".*[\\)].*")){}
				else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\[].*") && !Mention.matches(".*[\\]].*")){}
				else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\)].*") && !Mention.matches(".*[\\(].*")){}
				else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\]].*") && !Mention.matches(".*[\\[].*")){}
				else
				{
					AnnotationInPassage.add(start+"\t"+last+"\t"+Mention+"\t"+MentionType);
				}
				i--;
			}
			
			paragraph_num_last=paragraph_num;
			pmid_last=pmid;
		}// outputArr1
		AnnotationInPMID.add(AnnotationInPassage);
		GNormPlus.BioCDocobj.Annotations.add(AnnotationInPMID);
		
		//GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,false); //save in BioC file
	}
	
	public void ReadCRFresult(String Filename,String FilenameLoca,String FilenameOutput,String FilenameBioC,double threshold,double threshold_GeneType) throws XMLStreamException, IOException
	{
		/** load CRF output */
		ArrayList<String> outputArr1 = new ArrayList<String>();
		ArrayList<String> outputArr2 = new ArrayList<String>(); 
		ArrayList<String> outputArr3 = new ArrayList<String>();
		ArrayList<String> outputArr1_score = new ArrayList<String>(); 
		ArrayList<String> outputArr2_score = new ArrayList<String>(); 
		ArrayList<String> outputArr3_score = new ArrayList<String>(); 
		BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameOutput), "UTF-8"));
		String line;
		int rank=0;
		String score="";
		Pattern pat_Rank = Pattern.compile("^# ([0-2]) ([0-9\\.]+)$");
		while ((line = inputfile.readLine()) != null)  
		{
			Matcher mat = pat_Rank.matcher(line); // last column : Status
			if(mat.find())
			{
				rank = Integer.parseInt(mat.group(1));
				score = mat.group(2);
			}
			else if(rank == 0)
			{
				outputArr1.add(line);
				outputArr1_score.add(score);
			}
			else if(rank == 1)
			{
				outputArr2.add(line);
				outputArr2_score.add(score);
			}
			else if(rank == 2)
			{
				outputArr3.add(line);
				outputArr3_score.add(score);
			}
		}
		inputfile.close();
		
		/** load location */
		ArrayList<String> locationArr = new ArrayList<String>(); 
		inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(FilenameLoca), "UTF-8"));
		while ((line = inputfile.readLine()) != null)  
		{
			locationArr.add(line);
		}
		inputfile.close();
		
		/** output -> mentions */
		String pmid_last="";
		String paragraph_num_last="";
		String pmid="";
		String paragraph="";
		String paragraph_num="";
		Pattern pat_B = Pattern.compile("((FamilyName|DomainMotif|Gene)_[B])$");
		Pattern pat_IE = Pattern.compile("((FamilyName|DomainMotif|Gene)_[IE])$");
		ArrayList<ArrayList<String>> AnnotationInPMID = new ArrayList(); // array of Annotations in the PMIDs
		ArrayList<String> AnnotationInPassage= new ArrayList<String>(); // array of Annotations in the Passage
		GNormPlus.BioCDocobj.Annotations = new ArrayList();
		int countPMID=0;
		int countPassage=0;
		/** outputArr1 */
		int size_Arr=outputArr1.size();
		if(locationArr.size()<outputArr1.size())
		{
			size_Arr=locationArr.size();
		}
		for(int i=0;i<size_Arr;i++)
		{
			String outputsRow[]=outputArr1.get(i).split("\\t");
			String locationRow[]=locationArr.get(i).split("\\t");
			int start=100000;
			int last=0;
			String MentionType="";
			if(locationRow.length>3)
			{
				pmid=locationRow[0];
				paragraph=locationRow[1];
				paragraph_num=locationRow[2];
			}
			
			boolean F = false; //Flag of Finding
			if(outputsRow.length>=1)
			{
				Matcher mat = pat_B.matcher(outputsRow[outputsRow.length-1]); // last column : Status
				while(mat.find() && locationRow.length==6)
				{
					MentionType=mat.group(2);	
					pmid=locationRow[0];
					int start_tmp=Integer.parseInt(locationRow[4])-1;
					int last_tmp=Integer.parseInt(locationRow[5]);
					if(start_tmp<start){start=start_tmp;}
					if(last_tmp>last){last=last_tmp;}
					i++;
					outputsRow=outputArr1.get(i).split("\\t");
					locationRow=locationArr.get(i).split("\\t");
					mat = pat_IE.matcher(outputsRow[outputsRow.length-1]);
					F = true;
				}
			}
			
			if( (!paragraph_num_last.equals("")) && (!paragraph_num.equals(paragraph_num_last)) ) // paragraph change
			{
				AnnotationInPMID.add(AnnotationInPassage);
				AnnotationInPassage = new ArrayList<String>();
				countPassage++;
			}
			
			if( !pmid.equals(pmid_last) && paragraph_num.equals("0") && paragraph_num_last.equals("0") ) // pmid change (special case : the article only has one paragrpah)
			{
				AnnotationInPMID.add(AnnotationInPassage);
				AnnotationInPassage = new ArrayList<String>();
				GNormPlus.BioCDocobj.Annotations.add(AnnotationInPMID);
				AnnotationInPMID = new ArrayList();
				countPMID++;
				countPassage=0;
			}
			else if( (!pmid_last.equals("")) && (!pmid.equals(pmid_last)) ) // pmid change
			{
				GNormPlus.BioCDocobj.Annotations.add(AnnotationInPMID);
				AnnotationInPMID = new ArrayList();
				countPMID++;
				countPassage=0;
			}
			
			if(F == true)
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>countPMID && GNormPlus.BioCDocobj.PassageContexts.get(countPMID).size()>countPassage && GNormPlus.BioCDocobj.PassageContexts.get(countPMID).get(countPassage).length()>=last && (last-start)<1000)
				{
					String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(countPMID).get(countPassage); // Passage context
					String Mention = PassageContext.substring(start, last); 
					String Mention_nospace = Mention.replaceAll("[\\W\\-\\_]", ""); 
					if(Mention.toLowerCase().matches("(figure|tables|fig|tab|exp\\. [0-9]+).*")){}
					else if(Mention.matches("[A-Z][A-Z]s")){}
					else if(Mention.matches(".*\\|.*")){} 
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\;\\,\\'\\/\\\\].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\(].*") && !Mention.matches(".*[\\)].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\[].*") && !Mention.matches(".*[\\]].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\)].*") && !Mention.matches(".*[\\(].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\]].*") && !Mention.matches(".*[\\[].*")){}
					else if((GNormPlus.Abb2Longformtok_hash.containsKey(Mention_nospace.toLowerCase())) && (PassageContext.toLowerCase().matches(".*[\\W\\-\\-]("+GNormPlus.Abb2Longformtok_hash.get(Mention_nospace.toLowerCase())+")[\\W\\-\\-].*")))
					{
						//System.out.println(Mention_nospace.toLowerCase()+"\t"+GNormPlus.Abb2Longformtok_hash.get(Mention_nospace.toLowerCase()));
					}
					else
					{
						AnnotationInPassage.add(start+"\t"+last+"\t"+Mention+"\t"+MentionType);
					}
				}
				i--;
			}
			paragraph_num_last=paragraph_num;
			pmid_last=pmid;
		}// outputArr1
		AnnotationInPMID.add(AnnotationInPassage);
		GNormPlus.BioCDocobj.Annotations.add(AnnotationInPMID);
		
		/** outputArr2 */
		pmid_last="";
		paragraph_num_last="";
		pmid="";
		paragraph="";
		paragraph_num="";
		countPMID=0;
		countPassage=0; 
		size_Arr=outputArr2.size();
		if(locationArr.size()<outputArr2.size())
		{
			size_Arr=locationArr.size();
		}
		for(int i=0;i<size_Arr;i++)
		{
			String outputsRow[]=outputArr2.get(i).split("\\t");
			String locationRow[]=locationArr.get(i).split("\\t");
			int start=100000;
			int last=0;
			String MentionType="";
			if(locationRow.length>2)
			{
				pmid=locationRow[0];
				paragraph=locationRow[1];
				paragraph_num=locationRow[2];
			}
			
			boolean F = false; //Flag of Finding
			if(outputsRow.length>=1)
			{
				Matcher mat = pat_B.matcher(outputsRow[outputsRow.length-1]); // last column : Status
				while(mat.find() && locationRow.length==6)
				{
					MentionType=mat.group(2);	
					pmid=locationRow[0];
					int start_tmp=Integer.parseInt(locationRow[4])-1;
					int last_tmp=Integer.parseInt(locationRow[5]);
					if(start_tmp<start){start=start_tmp;}
					if(last_tmp>last){last=last_tmp;}
					i++;
					outputsRow=outputArr2.get(i).split("\\t");
					locationRow=locationArr.get(i).split("\\t");
					mat = pat_IE.matcher(outputsRow[outputsRow.length-1]);
					F = true;
				}
			}
			
			if( (!paragraph_num_last.equals("")) && (!paragraph_num.equals(paragraph_num_last)) ) // paragraph change
			{
				countPassage++;
			}
			
			if( !pmid.equals(pmid_last) && paragraph_num.equals("0") && paragraph_num_last.equals("0") ) // pmid change (special case : the article only has one paragrpah)
			{
				countPMID++;
				countPassage=0;
			}
			else if( (!pmid_last.equals("")) && (!pmid.equals(pmid_last)) ) // pmid change
			{
				countPMID++;
				countPassage=0;
			}
			
			if(F == true)
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>countPMID && GNormPlus.BioCDocobj.PassageContexts.get(countPMID).size()>countPassage && GNormPlus.BioCDocobj.PassageContexts.get(countPMID).get(countPassage).length()>=last && (last-start)<1000)
				{
					String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(countPMID).get(countPassage); // Passage context
					String Mention = PassageContext.substring(start, last); 
					String Mention_nospace = Mention.replaceAll("[\\W\\-\\_]", ""); 
					if(Mention.toLowerCase().matches("(figure|tables|fig|tab|exp\\. [0-9]+).*")){}
					else if(Mention.matches("[A-Z][A-Z]s")){}
					else if(Mention.matches(".*\\|.*")){} 
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\;\\,\\'\\/\\\\].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\(].*") && !Mention.matches(".*[\\)].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\[].*") && !Mention.matches(".*[\\]].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\)].*") && !Mention.matches(".*[\\(].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\]].*") && !Mention.matches(".*[\\[].*")){}
					else if((GNormPlus.Abb2Longformtok_hash.containsKey(Mention_nospace.toLowerCase())) && (PassageContext.toLowerCase().matches(".*[\\W\\-\\-]("+GNormPlus.Abb2Longformtok_hash.get(Mention_nospace.toLowerCase())+")[\\W\\-\\-].*")))
					{
						//System.out.println(Mention_nospace.toLowerCase()+"\t"+GNormPlus.Abb2Longformtok_hash.get(Mention_nospace.toLowerCase()));
					}
					else if(Double.parseDouble(outputArr2_score.get(i))>threshold)
					{
						boolean overlap=false;
						for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).size();j++)
						{
							String GetData[]=GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).get(j).split("\t");
							int startj=Integer.parseInt(GetData[0]);
							int lastj=Integer.parseInt(GetData[1]);
							String Mention_tmp = Mention.replaceAll("([^A-Za-z0-9@ ])", "\\\\$1");
							if(MentionType.equals("Gene") && Double.parseDouble(outputArr2_score.get(i))>threshold_GeneType && GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).get(j).matches(start+"\t"+last+"\t"+Mention_tmp+"\t(FamilyName|DomainMotif)") )
							{
								GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).set(j, start+"\t"+last+"\t"+Mention+"\t"+MentionType);
							}
							else if( (start>=startj && start<lastj) || (last>startj && last<=lastj) )
							{
								overlap=true;
							}
						}
						if(overlap == false)
						{
							GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).add(start+"\t"+last+"\t"+Mention+"\t"+MentionType);
						}
					}
				}
				i--;
			}
			
			paragraph_num_last=paragraph_num;
			pmid_last=pmid;
		}// outputArr2
		
		/** outputArr3 */
		pmid_last="";
		paragraph_num_last="";
		pmid="";
		paragraph="";
		paragraph_num="";
		countPMID=0;
		countPassage=0;
		size_Arr=outputArr3.size();
		if(locationArr.size()<outputArr3.size())
		{
			size_Arr=locationArr.size();
		}
		for(int i=0;i<size_Arr;i++)
		{
			String outputsRow[]=outputArr3.get(i).split("\\t");
			String locationRow[]=locationArr.get(i).split("\\t");
			int start=100000;
			int last=0;
			String MentionType="";
			if(locationRow.length>2)
			{
				pmid=locationRow[0];
				paragraph=locationRow[1];
				paragraph_num=locationRow[2];
			}
			
			boolean F = false; //Flag of Finding
			if(outputsRow.length>=1)
			{
				Matcher mat = pat_B.matcher(outputsRow[outputsRow.length-1]); // last column : Status
				while(mat.find() && locationRow.length==6)
				{
					MentionType=mat.group(2);	
					pmid=locationRow[0];
					paragraph_num=locationRow[2];
					int start_tmp=Integer.parseInt(locationRow[4])-1;
					int last_tmp=Integer.parseInt(locationRow[5]);
					if(start_tmp<start){start=start_tmp;}
					if(last_tmp>last){last=last_tmp;}
					i++;
					outputsRow=outputArr3.get(i).split("\\t");
					locationRow=locationArr.get(i).split("\\t");
					mat = pat_IE.matcher(outputsRow[outputsRow.length-1]);
					F = true;
				}
			}
			
			if( (!paragraph_num_last.equals("")) && (!paragraph_num.equals(paragraph_num_last)) ) // paragraph change
			{
				countPassage++;
			}
			
			if( !pmid.equals(pmid_last) && paragraph_num.equals("0") && paragraph_num_last.equals("0") ) // pmid change (special case : the article only has one paragrpah)
			{
				countPMID++;
				countPassage=0;
			}
			else if( (!pmid_last.equals("")) && (!pmid.equals(pmid_last)) ) // pmid change
			{
				countPMID++;
				countPassage=0;
			}
			
			if(F == true)
			{
				if(GNormPlus.BioCDocobj.PassageContexts.size()>countPMID && GNormPlus.BioCDocobj.PassageContexts.get(countPMID).size()>countPassage && GNormPlus.BioCDocobj.PassageContexts.get(countPMID).get(countPassage).length()>=last && (last-start)<1000)
				{
					String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(countPMID).get(countPassage); // Passage context
					String Mention = PassageContext.substring(start, last); 
					String Mention_nospace = Mention.replaceAll("[\\W\\-\\_]", ""); 
					if(Mention.toLowerCase().matches("(figure|tables|fig|tab|exp\\. [0-9]+).*")){}
					else if(Mention.matches("[A-Z][A-Z]s")){}
					else if(Mention.matches(".*\\|.*")){} 
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\;\\,\\'\\/\\\\].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\(].*") && !Mention.matches(".*[\\)].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\[].*") && !Mention.matches(".*[\\]].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\)].*") && !Mention.matches(".*[\\(].*")){}
					else if(Mention_nospace.length()<=3 && Mention.matches(".*[0-9].*") && Mention.matches(".*[\\]].*") && !Mention.matches(".*[\\[].*")){}
					else if((GNormPlus.Abb2Longformtok_hash.containsKey(Mention_nospace.toLowerCase())) && (PassageContext.toLowerCase().matches(".*[\\W\\-\\-]("+GNormPlus.Abb2Longformtok_hash.get(Mention_nospace.toLowerCase())+")[\\W\\-\\-].*")))
					{
						//System.out.println(Mention_nospace.toLowerCase()+"\t"+GNormPlus.Abb2Longformtok_hash.get(Mention_nospace.toLowerCase()));
					}
					else if(Double.parseDouble(outputArr3_score.get(i))>threshold)
					{
						boolean overlap=false;
						for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).size();j++)
						{
							String GetData[]=GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).get(j).split("\t");
							int startj=Integer.parseInt(GetData[0]);
							int lastj=Integer.parseInt(GetData[1]);
							String Mention_tmp = Mention.replaceAll("([^A-Za-z0-9@ ])", "\\\\$1");
							if(MentionType.equals("Gene") && Double.parseDouble(outputArr3_score.get(i))>threshold_GeneType && GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).get(j).matches(start+"\t"+last+"\t"+Mention_tmp+"\t(FamilyName|DomainMotif)") )
							{
								GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).set(j, start+"\t"+last+"\t"+Mention+"\t"+MentionType);
							}
							else if( (start>=startj && start<lastj) || (last>startj && last<=lastj) )
							{
								overlap=true;
							}
						}
						if(overlap == false)
						{
							GNormPlus.BioCDocobj.Annotations.get(countPMID).get(countPassage).add(start+"\t"+last+"\t"+Mention+"\t"+MentionType);
						}
					}
				}
				i--;
			}
			
			paragraph_num_last=paragraph_num;
			pmid_last=pmid;
		}// outputArr3
		
		//GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,false); //save in BioC file
	}
	
	public void PostProcessing(String Filename,String FilenameBioC) throws XMLStreamException, IOException
	{
		/** Develop Cell | FamilyName | DomainMotif lists */
		String Disease_Suffix="disease|diseases|syndrome|syndromes|tumor|tumour|deficiency|dysgenesis|atrophy|frame|dystrophy";
		String Cell_Suffix="cell|cells";
		String FamilyName_Suffix="disease|diseases|syndrome|syndromes|tumor|tumour|deficiency|dysgenesis|atrophy|frame|dystrophy|frame|factors|family|families|superfamily|superfamilies|subfamily|subfamilies|complex|genes|proteins";
		String DomainMotif_Suffix="domain|motif|domains|motifs|sequences";
		String Strain_Suffix="alpha|beta|gamma|kappa|theta|delta|[A-Ga-g0-9]";
		ArrayList<String> Translate2Family = new ArrayList<String>();
		
		for(int i=0;i<GNormPlus.BioCDocobj.Annotations.size();i++) // PMID
		{
			/** Pre-processing of  Mention type refinement */
			HashMap<String, String> Mention2Type_Hash = new HashMap<String, String>(); // for substring detection - Extract all mentions in the target PMID : MentionList
			ArrayList<String> GeneMentionPattern = new ArrayList<String>(); // pattern match to extend Gene
			HashMap<String, Integer> MentionType2Num = new HashMap<String, Integer>(); // for frequency calculation
			if(GNormPlus.BioCDocobj.PMIDs.size()>=i)
			{
				String pmid=GNormPlus.BioCDocobj.PMIDs.get(i);
				for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
				{	
					for(int k=0;k<GNormPlus.BioCDocobj.Annotations.get(i).get(j).size();k++) // Annotation
					{
						String Anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\\t");
						String start = Anno[0];
						String last = Anno[1];
						String mention = Anno[2];
						String type = Anno[3];
						Mention2Type_Hash.put(mention.toLowerCase(),type);
						if(MentionType2Num.containsKey(mention+"\t"+type))
						{
							MentionType2Num.put(mention.toLowerCase()+"\t"+type,MentionType2Num.get(mention+"\t"+type)+1);
							if( GNormPlus.PmidLF2Abb_lc_hash.containsKey(pmid+"\t"+mention.toLowerCase()) )
							{
								MentionType2Num.put(GNormPlus.PmidLF2Abb_lc_hash.get(pmid+"\t"+mention.toLowerCase())+"\t"+type,MentionType2Num.get(mention+"\t"+type)+1);
							}
							else
							{
								MentionType2Num.put(GNormPlus.PmidLF2Abb_lc_hash.get(pmid+"\t"+mention.toLowerCase())+"\t"+type,1);
							}
							if( GNormPlus.PmidAbb2LF_lc_hash.containsKey(pmid+"\t"+mention.toLowerCase()) )
							{
								MentionType2Num.put(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+mention.toLowerCase())+"\t"+type,MentionType2Num.get(mention+"\t"+type)+1);
							}
							else
							{
								MentionType2Num.put(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+mention.toLowerCase())+"\t"+type,1);
							}
						}
						else
						{
							MentionType2Num.put(mention.toLowerCase()+"\t"+type,1);
							if( GNormPlus.PmidLF2Abb_lc_hash.containsKey(pmid+"\t"+mention.toLowerCase()) )
							{
								MentionType2Num.put(GNormPlus.PmidLF2Abb_lc_hash.get(pmid+"\t"+mention.toLowerCase())+"\t"+type,1);
							}
							if( GNormPlus.PmidAbb2LF_lc_hash.containsKey(pmid+"\t"+mention.toLowerCase()) )
							{
								MentionType2Num.put(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+mention.toLowerCase())+"\t"+type,1);
							}
						}
						if(Anno[3].equals("Gene")) //Anno[3] is type
						{
							String mentmp=mention.toLowerCase();
							if(mentmp.matches(".*[0-9].*") || mentmp.matches(".*(alpha|beta|gamma|theta|zeta|delta).*"))
							{
								if(!mentmp.matches(".*\\{(alpha|beta|gamma|theta|zeta|delta)\\}.*"))
								{
									mentmp=mentmp.replaceAll("([^A-Za-z0-9\\| ])", "\\\\$1");
									mentmp=mentmp.replaceAll("[0-9]", "[0-9]");
									mentmp=mentmp.replaceAll("(alpha|beta|gamma|theta|zeta|delta)", "(alpha\\|beta\\|gamma\\|theta\\|zeta\\|delta)");
									if(!GeneMentionPattern.contains(mentmp))
									{
										GeneMentionPattern.add(mentmp);
									}
								}
							}
						}
					}
				}
							
				for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
				{
					ArrayList<Integer> RemoveList = new ArrayList<Integer>();
					for(int k=0;k<GNormPlus.BioCDocobj.Annotations.get(i).get(j).size();k++) // Annotation
					{
						String Anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\\t");
						String start = Anno[0];
						String last = Anno[1];
						String mention = Anno[2];
						String type = Anno[3];
						String mention_tmp=mention.toLowerCase().replaceAll("([^A-Za-z0-9@ ])", "\\\\$1");
						
						/**
						 * 1. Mention type refinement 
						 */
						/* 1. Check substring : update Gene -> Family name (TIF & TIF1) */
						boolean SubSt=false;
						/*
						// GDNFb -> GDNF (not work on 12682085_J_Cell_Biol_2003.xml)
						for (String men : Mention2Type_Hash.keySet())
						{
							if((!men.equals(mention.toLowerCase())) && men.matches(mention_tmp+"[\\W\\-\\_]*("+Strain_Suffix+")"))
							{
								GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\tFamilyName");
								if(GNormPlus.PmidLF2Abb_lc_hash.containsKey(GNormPlus.BioCDocobj.PMIDs.get(i)+"\t"+mention.toLowerCase()))
								{
									Translate2Family.add(GNormPlus.PmidLF2Abb_lc_hash.get(GNormPlus.BioCDocobj.PMIDs.get(i)+"\t"+mention.toLowerCase()));
								}
								else if(GNormPlus.PmidAbb2LF_lc_hash.containsKey(GNormPlus.BioCDocobj.PMIDs.get(i)+"\t"+mention.toLowerCase()))
								{
									Translate2Family.add(GNormPlus.PmidAbb2LF_lc_hash.get(GNormPlus.BioCDocobj.PMIDs.get(i)+"\t"+mention.toLowerCase()));
								}
								SubSt=true;
								break;
							}
						}
						*/
						if(SubSt == false)
						{
							int BoundaryLen=15;
							if(GNormPlus.BioCDocobj.PassageContexts.get(i).get(j).length()<Integer.parseInt(last)+15)
							{
								BoundaryLen=GNormPlus.BioCDocobj.PassageContexts.get(i).get(j).length()-Integer.parseInt(last);
							}
							String SurroundingString="";
							if(BoundaryLen<=0){}
							else if (GNormPlus.BioCDocobj.PassageContexts.get(i).get(j).length()<Integer.parseInt(last)+BoundaryLen)
							{
								BoundaryLen=GNormPlus.BioCDocobj.PassageContexts.get(i).get(j).length()-Integer.parseInt(last)-1;
								SurroundingString = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j).substring(Integer.parseInt(last),Integer.parseInt(last)+BoundaryLen).toLowerCase();
							}
							else
							{
								SurroundingString = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j).substring(Integer.parseInt(last),Integer.parseInt(last)+BoundaryLen).toLowerCase();
							}

							/* 2. Check suffix and surrounding words - Gene -> Family/Domain/Cell */
							if( mention.toLowerCase().matches(".*("+Cell_Suffix+")") ||	SurroundingString.matches("("+Cell_Suffix+")")	)
							{
								type="Cell";
								GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\t"+type);
							}
							else if( mention.toLowerCase().matches(".*("+FamilyName_Suffix+")") || SurroundingString.matches("("+FamilyName_Suffix+")")	)
							{
								type="FamilyName";
								GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\t"+type);
							}
							else if( mention.toLowerCase().matches(".*("+DomainMotif_Suffix+")")|| SurroundingString.matches("("+DomainMotif_Suffix+")")	)
							{
								type="DomainMotif";
								GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\t"+type);
							}
							else if(!type.equals("Gene"))
							{
								/* 3. Check (Family+Domain+Cell)/All rate (threshold = 0.5) - Family/Domain/Cell -> Gene */
								double Num_FDC=0;
								double Num_Gene=0;
								if(MentionType2Num.containsKey(mention.toLowerCase()+"\tFamilyName"))
								{
									Num_FDC = Num_FDC + MentionType2Num.get(mention.toLowerCase()+"\tFamilyName");
								}
								if(MentionType2Num.containsKey(mention.toLowerCase()+"\tDomainMotif"))
								{
									Num_FDC = Num_FDC + MentionType2Num.get(mention.toLowerCase()+"\tDomainMotif");
								}
								if(MentionType2Num.containsKey(mention.toLowerCase()+"\tCell"))
								{
									Num_FDC = Num_FDC + MentionType2Num.get(mention.toLowerCase()+"\tCell");
								}
								if(MentionType2Num.containsKey(mention.toLowerCase()+"\tGene"))
								{
									Num_Gene = Num_Gene + MentionType2Num.get(mention.toLowerCase()+"\tGene");
								}
								if(Num_Gene/(Num_FDC+Num_Gene)>=0.5)
								{
									GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\tGene");
								}
								
								/* 4. Extend Genes to Family/Domain mentions by pattern match - Family/Domain/Cell -> Gene */
								for(int p=0;p<GeneMentionPattern.size();p++)
								{
									if(mention.toLowerCase().matches(GeneMentionPattern.get(p)))
									{
										GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\tGene");
									}
								}
							}
							
							/*
							 *  5. Abbreviation resolution
							 *   - LF/Abb extracted : LF.type -> Abb.type 
							 *   - Abb only : Abb.type -> LF.type
							 *   - LF only : LF.type -> Abb.type
							 */
							String lc_ment=mention.toLowerCase();
							if(GNormPlus.PmidAbb2LF_lc_hash.containsKey(pmid+"\t"+lc_ment)) //the target mention is abbreviation
							{
								//Infer Abbreviation by Long form
								if(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+lc_ment).matches(".*("+Disease_Suffix+")"))
								{
									//remove the mention (Abb), because the LF is a disease 
								}
								else if(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+lc_ment).matches(".*("+Cell_Suffix+")"))
								{
									//GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Anno[0]+"\t"+Anno[1]+"\tCell");
								}
								else if(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+lc_ment).matches(".*("+FamilyName_Suffix+")") && !lc_ment.matches(".+[a-z][0-9][a-z]")) //AtRPA1a in pmid:19153602
								{
									GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\tFamilyName");
								}
								else if(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+lc_ment).matches(".*("+DomainMotif_Suffix+")"))
								{
									GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\tDomainMotif");
								}
								else
								{
									if(Mention2Type_Hash.containsKey(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+lc_ment)) 
									&& Mention2Type_Hash.get(GNormPlus.PmidAbb2LF_lc_hash.get(pmid+"\t"+lc_ment)).equals("Gene")
									&& !(type.equals("Gene"))
									) // if Long Form is recognized as a Gene, and Abb is recognized as not a Gene
									{
										GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, start+"\t"+last+"\t"+mention+"\tGene");
									}
								}
							}
						} //if(Remov == true)
					}
				}
				
				for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
				{
					for(int k=0;k<GNormPlus.BioCDocobj.Annotations.get(i).get(j).size();k++) // Annotation
					{
						String Anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\\t");
						if(Translate2Family.contains(Anno[2].toLowerCase()))
						{
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Anno[0]+"\t"+Anno[1]+"\t"+Anno[2]+"\tFamilyName");
						}
					}
				}
				
				//ArrayList<String> GeneMentionPattern = new ArrayList<String>(); // pattern match to extend Gene
				HashMap<String,String> GeneMentions = new HashMap<String,String>(); // Extending Gene mentions
				HashMap<String,String> GeneMentionLocationGNR = new HashMap<String,String>(); // Extending Gene mentions
				for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
				{
					for(int k=0;k<GNormPlus.BioCDocobj.Annotations.get(i).get(j).size();k++) // Annotation
					{
						String Anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\\t");
						int start = Integer.parseInt(Anno[0]);
						int last = Integer.parseInt(Anno[1]);
						String mention = Anno[2];
						String type = Anno[3];
						for (int s=start ;s<=last;s++)
						{
							GeneMentionLocationGNR.put(j+"\t"+s,type);
						}
						GeneMentions.put(mention.toLowerCase(),type);
					}
				}
				
				//Extend to all gene mentions
				for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
				{
					if(GNormPlus.BioCDocobj.PassageContexts.size()>i && GNormPlus.BioCDocobj.PassageContexts.get(i).size()>j)
					{
						String PassageContexts = " " + GNormPlus.BioCDocobj.PassageContexts.get(i).get(j) + " ";
						String PassageContexts_tmp = PassageContexts.toLowerCase();
						for(String gm : GeneMentions.keySet())
						{
							String type=GeneMentions.get(gm);
							if(type.equals("Gene"))
							{
								gm = gm.replaceAll("([\\W\\-\\_])", "\\\\$1");
								gm=gm.replaceAll("[0-9]", "\\[0\\-9\\]");
								gm=gm.replaceAll("(alpha|beta|gamma|theta|zeta|delta)", "(alpha\\|beta\\|gamma\\|theta\\|zeta\\|delta)");
								gm=gm.replaceAll("\\-[a-z]$", "\\-\\[a\\-z\\]");
								Pattern ptmp = Pattern.compile("^(.*[\\W\\-\\_])("+gm+")([\\W\\-\\_].*)$");
								Matcher mtmp = ptmp.matcher(PassageContexts_tmp);
								while(mtmp.find())
								{
									String pre = mtmp.group(1);
									String gmtmp = mtmp.group(2);
									String post = mtmp.group(3);
			
									int start = pre.length()-1;
									int last = start+gmtmp.length();
									if(PassageContexts.length()>last)
									{
										String mention = PassageContexts.substring(start+1,last+1);
										if(!GeneMentionLocationGNR.containsKey(j+"\t"+start) && !GeneMentionLocationGNR.containsKey(j+"\t"+last))
										{
											if(GNormPlus.BioCDocobj.Annotations.get(i).get(j).contains(start+"\t"+last+"\t"+mention+"\tFamilyName"))
											{
												GNormPlus.BioCDocobj.Annotations.get(i).get(j).remove(start+"\t"+last+"\t"+mention+"\tFamilyName");
											}
											else if(GNormPlus.BioCDocobj.Annotations.get(i).get(j).contains(start+"\t"+last+"\t"+mention+"\tDomainMotif"))
											{
												GNormPlus.BioCDocobj.Annotations.get(i).get(j).remove(start+"\t"+last+"\t"+mention+"\tDomainMotif");
											}
											GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tGene");
										}
										gmtmp = gmtmp.replaceAll(".", "X");
										PassageContexts_tmp=pre+""+gmtmp+""+post;
										mtmp = ptmp.matcher(PassageContexts_tmp);
									}
								}
							}
						}
					}
				}
				
				//Extend to all family mentions
				for(int j=0;j<GNormPlus.BioCDocobj.Annotations.get(i).size();j++) // Paragraph
				{
					if(GNormPlus.BioCDocobj.PassageContexts.size()>i && GNormPlus.BioCDocobj.PassageContexts.get(i).size()>j)
					{
						String PassageContexts = " " + GNormPlus.BioCDocobj.PassageContexts.get(i).get(j) + " ";
						String PassageContexts_tmp = PassageContexts.toLowerCase();
						for(String gm : GeneMentions.keySet())
						{
							String type=GeneMentions.get(gm);
							if(type.matches("(FamilyName|DomainMotif)"))
							{
								gm = gm.replaceAll("([\\W\\-\\_])", "\\\\$1");
								gm=gm.replaceAll("s$", "(s\\|)");
								Pattern ptmp = Pattern.compile("^(.*[\\W\\-\\_])("+gm+")([\\W\\-\\_].*)$");
								Matcher mtmp = ptmp.matcher(PassageContexts_tmp);
								while(mtmp.find())
								{
									String pre = mtmp.group(1);
									String gmtmp = mtmp.group(2);
									String post = mtmp.group(3);
									
									int start = pre.length()-1;
									int last = start+gmtmp.length();
									if(PassageContexts.length()>last)
									{
										String mention = PassageContexts.substring(start+1,last+1);
										if(!GeneMentionLocationGNR.containsKey(j+"\t"+start) && !GeneMentionLocationGNR.containsKey(j+"\t"+last))
										{
											if(!GNormPlus.BioCDocobj.Annotations.get(i).get(j).contains(start+"\t"+last+"\t"+mention+"\tGene"))
											{
												GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\t"+type);
											}
										}
										gmtmp = gmtmp.replaceAll(".", "X");
										PassageContexts_tmp=pre+""+gmtmp+""+post;
										mtmp = ptmp.matcher(PassageContexts_tmp);
									}
								}
							}
						}
					}
				}
			}
		}
		GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false,false); //save in BioC file
	}
}


